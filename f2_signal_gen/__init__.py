# f2_signal_gen class 
# Notes, this is an _Generator_ this system should not consider if the signal 
# is generated for downlink or uplink
# Assumptions:
#   Every transmitter may have multiple TX antennas
#   Every transmitter has the same number of antennas
#   Users can be in the same (Downlink) of in different (Uplink) transmitter
#   Generator does not take into account where the user signals are merged
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 06.08.2018 18:23
import sys
import numpy as np
import scipy.signal as sig
import tempfile
import subprocess
import shlex
import time
import functools as fnc
from thesdk import *
from refptr import *
import modem as mdm #Function definitions
from signal_generator_802_11n import *

#Simple buffer template
class f2_signal_gen(thesdk):
    def __init__(self,*arg): 
        self.proplist=[]
        self.Rs = 80e6                          #Sample frequency
        self.Txantennas=4                       #Number of transmitting antennas
        self.Txpower=30                         #Output power per antenna in dBm
        self.Users=2                            #Number of users
        self.Digital='False'                    #If true, the ouput is quantized to Bits
        self.Bits=10
        self.Digital_mode='2C'                  #Two's complement
        self.Disableuser=[]
        self.Disableuser= [ self.Disableuser.append(False) for i in range(self.Users) ]         #Disable data transmission for cerrtain users
        self.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14, 'BBRs':40e6 };  #Mode of the baseband signal. Let's start with sinusoids

        self.model='py';                        #can be set externally, but is not propagated
        self._filterlist=[]                     #list of interpolation filters
        self._qam_reference=[]                  #Reference symbols stream for EVM
        self._bitstream_reference=[]            #Reference bit stream for BER 
        self._Z = refptr();
        self._classfile=__file__                #needed only if rtl defined as superclass
        self.DEBUG= False
        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;
        self.init()

    def init(self):
        #adding the 802.11n generator
        self.sg802_11n=signal_generator_802_11n(self)

        if self.bbsigdict['mode']=='sinusoid':
            self.sinusoid()
        if self.bbsigdict['mode']=='ofdm_sinusoid':
            self.ofdm_sinusoid()
            self._qam_reference=self.sg802_11n._qam_reference
        if self.bbsigdict['mode']=='ofdm_random_qam':
            self.ofdm_random_qam()
            self._qam_reference=self.sg802_11n._qam_reference
            self._bitstream_reference=self.sg802_11n._bitstream_reference
        if self.bbsigdict['mode']=='ofdm_random_802_11n':
            self.ofdm_random_802_11n()
            self._qam_reference=self.sg802_11n._qam_reference
            self._bitstream_reference=self.sg802_11n._bitstream_reference
        
        if self.Digital=='True':
            digitize_argdict={'signal':self._Z.Value, 'Bits':self.Bits, 
                    'Scale':self.Txpower, 'mode':self.Digital_mode }
            if digitize_argdict['Scale']>0:
                self.print_log({'type':'I', 'msg':"Digitizer scale > 0dB. Defaulting to 0 dB"})
                digitize_argdict['Scale']=0
            self._Z.Value=digitize(**digitize_argdict)

    def run(self): #Just an alias for init to be consistent: run() executes the core function
        self.init()

    #Methods for signal generation. Add a function and add it to init()
    #controlled with bbsigdict
    def sinusoid(self):
            length=self.bbsigdict['length']
            phi=np.transpose(np.array(np.mat(self.bbsigdict['freqs'])))*np.array(range(length))*2*np.pi/(self.bbsigdict['BBRs'])
            usersig=np.transpose(np.ones((self.Txantennas,1))*np.sum(np.exp(1j*phi),0)/len(self.bbsigdict['freqs'])) #All antennas emit the same signal
            #All users have the same signal
            out=np.zeros((self.Users,usersig.shape[0],usersig.shape[1]),dtype='complex')
            for i in range(self.Users):
                out[i,:,:]=usersig
            
            out=self.interpolate_at_antenna({'signal':out})
            self._Z.Value=out 

    def ofdm_sinusoid(self):
         self.sg802_11n.ofdm_sinusoid()
         self._Z=self.sg802_11n._Z

    
    def ofdm_random_802_11n(self):
         out=self.sg802_11n.gen_random_802_11n_ofdm()
         out=self.interpolate_at_antenna({'signal':out})
         self.print_log({'type':'D', 'msg':out.shape})
         self.print_log({'type':'D', 'msg':"Test"})
         test=out[0,320+16:320+80,0]
         test.shape=(-1,1)
         self.print_log({'type':'D', 'msg':test.shape})
         test=np.fft.fft(test,axis=0)/64
         self.print_log({'type':'D', 'msg':test[Freqmap]})
         self._Z.Value=out
         #self.print_log({'type':'D', 'msg':self._Z.Value[0,320+16:320+80,0]})
         test=self._Z.Value[0,320+16:320+80,0]
         test.shape=(-1,1)
         self.print_log({'type':'D', 'msg':test.shape})
         test=np.fft.fft(test,axis=0)/64
         self.print_log({'type':'D', 'msg':test[Freqmap]})

    def ofdm_random_qam(self):
         self.sg802_11n.ofdm_random_qam()
         self._Z=self.sg802_11n._Z

    def set_transmit_power(self):
         t=[]
         for user in range(self._Z.Value.shape[0]):
             for antenna in range(self._Z.Value.shape[2]):
                 if not self.sg802_11n.Disableuser[user]:
                     t=np.r_['0',t, self._Z.Value[user,:,antenna]]
         Vrmscurrent=np.std(t)

         Vrms=np.sqrt(1e-3*50*10**(self.Txpower/10))
         for user in range(self._Z.Value.shape[0]):
             for antenna in range(self._Z.Value.shape[2]):
                 msg="Setting transmit Rms signal amplitude to from %f to %f Volts corresponding to %f dBm transmit power to 50 ohms" %(float(Vrmscurrent), float(Vrms), float(self.Txpower))
                 self.print_log({'type':'I', 'msg': msg}) 
                 self._Z.Value[user,:,antenna]=self._Z.Value[user,:,antenna]/Vrmscurrent*Vrms

    def interpolate_at_antenna(self,argdict={'signal':[]}):
        ratio=self.Rs/self.bbsigdict['BBRs']
        signal=argdict['signal']
        #Currently fixeed interpolation. check the function definitions for details
        factors=factor({'n':ratio})
        msg="Interpolation factors at antenna are %s" %(factors)
        self.print_log({'type':'I', 'msg': msg}) 
        filterlist=self.generate_interpolation_filterlist({'interp_factor':ratio})
        msg="Signal length is now %i" %(signal.shape[1])
        self.print_log({'type':'I', 'msg': msg}) 
        #This is to enable growth of the signal length that better mimics the hardware
        #sig.resample_poly is more effective, but does not allow growth.
        for user in range(signal.shape[0]):
            for antenna in range(signal.shape[2]):
                t=signal[user,:,antenna]
                for i in range(factors.shape[0]):
                #signali=sig.resample_poly(signal, fact, 1, axis=1, window=fircoeffs)
                    #signali=sig.resample_poly(signal[user,:,antenna], fact, 1, window=i)
                    t2=np.zeros((int(t.shape[0]*factors[i])),dtype='complex')
                    t2[0::int(factors[i])]=t
                    t=sig.convolve(t2, filterlist[i],mode='full')
                if user==0 and antenna==0:
                    signali=np.zeros((signal.shape[0],t.shape[0],signal.shape[2]),dtype='complex')
                    signali[user,:,antenna]=t
                else:
                    signali[user,:,antenna]=t
        msg="Signal length is now %i" %(signali.shape[1])
        self.print_log({'type':'I', 'msg': msg}) 
        self._filterlist=filterlist
        return signali

    def generate_interpolation_filterlist(self,argdict={'interp_factor':1}):
        #Use argument dictionary. Makes modifications easier.
        interp_factor=argdict['interp_factor']
        
        attenuation=70 #Desired attenuation in decibels
        factors=factor({'n':interp_factor})
        #self.print_log({'type':'D', 'msg':factors})
        fsample=1
        BW=0.45
        numtaps=65     # TAps for the first filterThis should be somehow verified
        #Harris rule. This is to control stability of Rmez
        #numtaps= int(np.ceil(attenuation*fsample*factors[0]/(fsample/2-BW)))    # TAps for the first filterThis should be somehow verified
        desired=np.array([ 1, 10**(-attenuation/10)] )
        #check the mask specs from standard 
        #mask=np.array([ 1, 10**(-28/10)    ] )
        filterlist=list()
        if interp_factor >1:
            for i in factors:
                fact=i
                #self.print_log({'type':'D', 'msg':fsample})
                if  fsample/(0.5)<= 8: #FIR is needed
                    msg= "BW to sample rate ratio is now %s" %(fsample/0.5)
                    self.print_log({'type': 'I', 'msg':msg })
                    msg="Interpolation by %i" %(fact)
                    self.print_log({'type': 'I', 'msg':msg })
                    bands=np.array([0, BW, (fsample*fact/2-BW), fact*fsample/2])
                    filterlist.append(sig.remez(numtaps, bands, desired, Hz=fact*fsample))
                    fsample=fsample*fact #increase the sample frequency
                    numtaps=np.amax([3, int(np.floor(numtaps/fact)) + int((np.floor(numtaps/fact)%2-1))]) 
                else:
                    self.print_log({'type':'I', 'msg':"BW to sample rate ratio is now %s" %(fsample/0.5)})
                    fact=fnc.reduce(lambda x,y:x*y,factors)/fsample
                    self.print_log({'type':'I', 'msg':"Interpolation with 3-stage CIC-filter by %i" %(fact)})
                    fircoeffs=np.ones(int(fact))/(fact) #do the rest of the interpolation with 3-stage CIC-filter
                    fircoeffs=fnc.reduce(lambda x,y: np.convolve(x,y),list([fircoeffs, fircoeffs, fircoeffs]))
                    filterlist.append(fircoeffs)
                    #self.print_log({'type':'D', 'msg':filterlist})
                    fsample=fsample*fact #increase the sample frequency
                    self.print_log({'type':'I', 'msg':"BW to sample rate ratio is now %s" %(fsample/0.5)})
                    break
        else:
            self.print_log({'type':'I', 'msg':"Interpolation ratio is 1. Generated unit coefficient"})
            filterlist.append([1.0]) #Ensure correct operation in unexpected situations.
        return filterlist
    

#Funtion definitions
def digitize(**kwargs):
    signal=kwargs.get('signal')
    bits=kwargs.get('Bits',10)
    #Scale factor in decibels. By default, signal is scaled to cover the full dynamic range
    #0 means full scale.
    scale=10**(kwargs.get('Scale',0)/20)
    mode=kwargs.get('Mode','2C') #2C =Two's complement, BO is binary offset
    max=np.amax(np.abs(np.r_['1', np.real(signal), np.imag(signal)]))
    #Default is two's complement, i.e. negative numbers remain negative
    digitized=np.round(signal/max*scale*(2**(bits-1)-1))
    if mode=='BO':
       #Not rescaled in order to retain symmetry relatively to the midpoint
       digitized=digitized+2**(bits-1)-1
    return digitized
        

def factor(argdict={'n':1}):
    #This is a function to calculate factors of an integer as in Matlab
    # "Everything in Matlab is available in Python" =False.
    reminder=argdict['n']
    factors=np.array([])
    while reminder >1:
        notfound=True
        for i in range(2,int(np.sqrt(reminder)+2)):
            if reminder % i == 0:
                factors=np.r_[ factors, i]
                reminder=reminder/i
                notfound=False
                break
        if notfound==True:
                factors=np.r_[ factors, reminder]
                reminder=0
    return factors


if __name__=="__main__":
    import scipy as sci
    import numpy as np
    import matplotlib.pyplot as plt
    from thesdk import *
    from  f2_signal_gen import *
    t=f2_signal_gen()
    t.Rs=8*11*20e6
    t.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14, 'BBRs': 20e6 }; #Mode of the baseband signal. Let's start with sinusoids
    t.Txantennas=1                       #Number of transmitting antennas
    t.Txpower=0                          #Output power per antenna in dBm
    t.Users=4                            #Number of users
    t.Digital='True'
    t.DEBUG='True'
    t.init()
    #t.set_transmit_power()
    #self.print_log({'type':'D', 'msg':np.std(t._Z.Value,axis=1)})
    #self.print_log({'type':'D', 'msg':t._Z.Value})
    #self.print_log({'type':'D', 'msg':t._Z.Value.shape})
    #n=t._Z.Value/np.std(t._Z.Value,axis=1)
    t.print_log({'type':'D', 'msg':np.max(t._Z.Value)})
    t.print_log({'type':'D', 'msg':t._Z.Value.shape})
    #self.print_log({'type':'D', 'msg':filt})
    #self.print_log({'type':'D', 'msg':filt.shape})
    #tf=factor(8)
    #self.print_log({'type':'D', 'msg':tf})
    #tf=factor(80)
    #self.print_log({'type':'D', 'msg':tf})
    filt=t._filterlist
    for i in filt:
        w, h = sig.freqz(i)
        plt.plot(w/(2*np.pi), 20*np.log10(abs(h)))
    plt.show()

    
