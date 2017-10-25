# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 24.10.2017 17:02
import sys
sys.path.append ('/home/projects/fader/TheSDK/Entities/refptr/py')
sys.path.append ('/home/projects/fader/TheSDK/Entities/thesdk/py')
sys.path.append ('/home/projects/fader/TheSDK/Entities/modem/py')
import numpy as np
import scipy.signal as sig
import tempfile
import subprocess
import shlex
import time
import functools as fnc
import modem as mdm #Function definitions

#Classes are distinguished by the class
from refptr import *
from thesdk import *
#from rtl import *

#Simple buffer template
class f2_signal_gen(thesdk):

    def __init__(self,*arg): 
        self.proplist = [ 'Txantennas', 'Txpower', 'Users', 'Rs', 'bbsigdict', 'ofdmdict' ];  #Properties that can be propagated from parent
        self.Rs = 80e6                          #Sample frequency
        self.Txantennas=4                       #Number of transmitting antennas
        self.Txpower=30                         #Output power per antenna in dBm
        self.Users=2                            #Number of users
        self.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14, 'BBRs':20e6 }; #Mode of the baseband signal. Let's start with sinusoids
        self.ofdmdict={ 'framelen':64,'data_loc': np.r_[1:11+1, 13:25+1, 
        27:39+1, 41:53+1, 55:64+1]-1, 'pilot_loc' : np.r_[-21, -7, 7, 21] + 32, 'CPlen':16}

        self.model='py';                        #can be set externally, but is not propagated
        self._filterlist=[]                     #list of interpolation filters
        self._Z = refptr();
        self._classfile=__file__ #needed only if rtl defined as superclass

        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;
        self.init()

    def init(self):
        if self.bbsigdict['mode']=='sinusoid':
            self.sinusoid()
        if self.bbsigdict['mode']=='ofdm_sinusoid':
            self.ofdm_sinusoid()
        if self.bbsigdict['mode']=='ofdm_random_qam':
            self.ofdm_random_qam()

    def run(self): #Just an alias for init to be consistent: run() always executes the core function
        self.init()


    #Methods for signal generation. Add a function and add it to init()
    #controlled with bbsigdict
    def sinusoid(self):
            length=self.bbsigdict['length']
            phi=np.transpose(np.array(np.mat(self.bbsigdict['freqs'])))*np.array(range(length))*2*2*np.pi/(self.Rs)
            usersig=np.transpose(np.ones((self.Txantennas,1))*np.sum(np.exp(1j*phi),0)/len(self.bbsigdict['freqs'])) #All antennas emit the same signal
            #All users have the same signal
            out=np.zeros((self.Users,usersig.shape[0],usersig.shape[1]),dtype='complex')
            for i in range(self.Users):
                out[i,:,:]=usersig
            
            out=self.interpolate_at_antenna({'signal':out})
            self._Z.Value=out 

    def ofdm_sinusoid(self):
            ofdmdict=self.ofdmdict
            length=np.floor(self.bbsigdict['length']/ofdmdict['framelen'])
            signalindexes=np.round(np.array(self.bbsigdict['freqs'])/self.Rs*ofdmdict['framelen']).astype(int)
            frame=np.zeros((length,ofdmdict['framelen']))
            frame[:,signalindexes]=1
            datasymbols=frame[:,ofdmdict['data_loc']]
            pilotsymbols=frame[:,ofdmdict['pilot_loc']]
            usersig=np.ones((self.Txantennas,1))*mdm.ofdmMod(ofdmdict,datasymbols,pilotsymbols).T
            
            out=np.zeros((self.Users,usersig.shape[0],usersig.shape[1]),dtype='complex')
            for i in range(self.Users):
                out[i,:,:]=usersig

            self._Z.Value=out 



    
    def ofdm_random_qam(self):
            #Local vars just to clear to code
            ofdmdict=self.ofdmdict
            bbsigdict=self.bbsigdict
            framelen=ofdmdict['framelen']
            length=bbsigdict['length']
            CPlen=ofdmdict['CPlen']
            QAM=bbsigdict['QAM']
            BBRs=bbsigdict['BBRs']
            #The length is approx this many frames
            frames=np.floor(length/(framelen+CPlen))
            bitspersymbol=np.log2(QAM).astype(int)
            
            #generate random bitstreams per antenna
            bitstream=np.random.randint(2,size=(self.Txantennas,frames*bitspersymbol*framelen))
            #Init the qam signal, frame and out
            qamsignal=np.zeros((self.Txantennas,frames*framelen),dtype='complex')
            frame=np.zeros((frames,framelen),dtype='complex')
            #usersig=np.zeros((self.Txantennas,frames*(framelen+CPlen)),dtype='complex')

            for i in range(self.Txantennas):
                wordstream, qamsignal[i]= mdm.qamModulateBitStream(bitstream[i], QAM)
                qamsignal[i]=qamsignal[i].reshape((1,qamsignal.shape[1]))
                frame= qamsignal[i].reshape((-1,framelen))
                datasymbols=frame[:,ofdmdict['data_loc']]
                pilotsymbols=frame[:,ofdmdict['pilot_loc']]
                interpolated=mdm.ofdmMod(ofdmdict,datasymbols,pilotsymbols)
                interpolated=interpolated.reshape((-1,(framelen+CPlen)))
                interpolated=self.interpolate_at_ofdm({'signal':interpolated})
                win=window({'Tr':100e-9, 'length':interpolated.shape[1], 'fs':80e6, 'duration': (framelen+CPlen)*self.Rs/BBRs})
                interpolated=interpolated*win               
                print(interpolated.shape)
                if i==0:
                    usersig=np.zeros((self.Txantennas,interpolated.shape[0]*interpolated.shape[1]),dtype='complex')
                    usersig[i,:]=interpolated.reshape(1,-1)
                else:
                    usersig[i,:]=interpolated.reshape(1,-1)
                
            usersig=usersig.T #usersig.shape[0] is time
            print(usersig.shape)
            out=np.zeros((self.Users,usersig.shape[0],usersig.shape[1]))
            for i in range(self.Users):
                print(usersig.shape)
                out[i,:,:]=usersig

            self._Z.Value=out 

    def set_transmit_power(self):
         for user in range(self._Z.Value.shape[0]):
             for antenna in range(self._Z.Value.shape[2]):
                 self._Z.Value[user,:,antenna]=self._Z.Value[user,:,antenna]/np.std(self._Z.Value[user,:,antenna])*np.sqrt(10**(self.Txpower/10)*1e-3*50) #Rms voltage normalized to 50 ohm
                 #print(np.std(self._Z.Value[user,:,antenna]))


    def interpolate_at_antenna(self,argdict={'signal':[]}):
        ratio=self.Rs/self.bbsigdict['BBRs']
        signal=argdict['signal']
        #Currently fixeed interpolation. check the fucntion definitions for details
        factors=factor(ratio)
        filterlist=generate_interpolation_filterlist({'interp_factor':ratio})
        print("Signal length is now %i" %(signal.shape[1]))
        #This is to enable growth of the signal length that better mimics the hardware
        #sig.resample_poly is more effective, but does not allow growth.
        for user in range(signal.shape[0]):
            for antenna in range(signal.shape[2]):
                t=signal[user,:,antenna]
                for i in range(factors.shape[0]):
                #signali=sig.resample_poly(signal, fact, 1, axis=1, window=fircoeffs)
                    #signali=sig.resample_poly(signal[user,:,antenna], fact, 1, window=i)
                    t2=np.zeros((t.shape[0]*factors[i]))
                    t2[0::factors[i]]=t
                    t=sig.convolve(t2, filterlist[i],mode='full')
                if user==0 and antenna==0:
                    signali=np.zeros((signal.shape[0],t.shape[0],signal.shape[2]),dtype='complex')
                    signali[user,:,antenna]=t
                else:
                    signali[user,:,antenna]=t
        print("Signal length is now %i" %(signali.shape[1]))
        self._filterlist=filterlist
        return signali

    def interpolate_at_ofdm(self,argdict={'signal':[]}):
        ratio=self.Rs/self.bbsigdict['BBRs']
        signal=argdict['signal']
        #Currently fixeed interpolation. check the fucntion definitions for details
        factors=factor(ratio)
        filterlist=generate_interpolation_filterlist({'interp_factor':ratio})
        print("Signal length is now %i" %(signal.shape[1]))
        #This is to enable growth of the signal length that better mimics the hardware
        #sig.resample_poly is more effective, but does not allow growth.
        for symbol in range(signal.shape[0]):
            t=signal[symbol,:]
            for i in range(factors.shape[0]):
                #signali=sig.resample_poly(signal[user,:,antenna], fact, 1, window=i)
                t2=np.zeros((t.shape[0]*factors[i]),dtype='complex')
                t2[0::factors[i]]=t
                t=sig.convolve(t2, filterlist[i],mode='full')
            if symbol==0:
                signali=np.zeros((signal.shape[0],t.shape[0]),dtype='complex')
                signali[symbol,:]=t
            else:
                signali[symbol,:]=t
        print("Signal length is now %i" %(signali.shape[1]))
        self._filterlist=filterlist
        print(filterlist)
        return signali

#Window to taper OFDM symbols
def window(argdict={'Tr':100e-9, 'length':478, 'fs':80e6, 'duration':240 }):
    #Tr is the window rise/fall time 100ns is from the specification
    Tr=argdict['Tr']
    fs=argdict['fs']
    length=argdict['length']
    duration=argdict['duration']
    T=argdict['duration']/fs
    prefix=(length-duration)/2
    t=(np.arange(length)-prefix)/fs # Time: T is the duration of the symbol without tails
    window=np.zeros(length)
    for i in range(len(t)):
     if t[i]>= -Tr/2 and t[i]<=Tr/2:
         window[i]=np.sin(np.pi/2*(0.5+t[i]/Tr))**2
     elif t[i]>Tr/2 and t[i]<T-Tr/2:
         window[i]=1
     elif t[i]>= T-Tr/2 and t[i] <= T+Tr/2:
         window[i]=np.sin(np.pi/2*(0.5-(t[i]-T)/Tr))**2
    return window

#Funtion definitions
def generate_interpolation_filterlist(argdict={'interp_factor':1}):
    #Use argument dictionary. Makes modifications easier.
    interp_factor=argdict['interp_factor']

    attenuation=70 #Desired attenuation in decibels
    numtaps=65     # TAps for the first filterThis should be somehow verified
    factors=factor(interp_factor)
    print(factors)
    fsample=1
    BW=0.45
    desired=np.array([ 1, 10**(-attenuation/10)] )
    filterlist=list()
    if interp_factor >1:
        for i in factors:
            fact=i
            print(fsample)
            if  fsample/(0.5)<= 8: #FIR is needed
                print("BW to sample rate ratio is now %s" %(fsample/0.5))
                print("Interpolation by %i" %(fact))
                bands=np.array([0, BW, (fsample*fact/2-BW), fact*fsample/2])
                filterlist.append(sig.remez(numtaps, bands, desired, Hz=fact*fsample))
                fsample=fsample*fact #increase the sample frequency
                numtaps=np.amax([3, int(np.floor(numtaps/fact)) + int((np.floor(numtaps/fact)%2-1))]) 
            else:
                print("BW to sample rate ratio is now %s" %(fsample/0.5))
                fact=fnc.reduce(lambda x,y:x*y,factors)/fsample
                print("Interpolation with 3-stage CIC-filter by %i" %(fact))
                fircoeffs=np.ones(int(fact))/(fact) #do the rest of the interpolation with 3-stage CIC-filter
                fircoeffs=fnc.reduce(lambda x,y: np.convolve(x,y),list([fircoeffs, fircoeffs, fircoeffs]))
                filterlist.append(fircoeffs)
                print(filterlist)
                fsample=fsample*fact #increase the sample frequency
                print("BW to sample rate ratio is now %s" %(fsample/0.5))
                break
    else:
        print("Interpolation ratio is 1. Generated unit coefficient")
        filterlist.append([1.0]) #Ensure correct operation in unexpected situations.
    return filterlist



def factor(n):
    #This is a function to calculate factors of an integer as in Matlab
    # "Everything in Matlab is available in Python" =False.
    reminder=n
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
    from  f2_signal_gen import *
    t=f2_signal_gen()
    t.Rs=8*11*20e6
    t.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14, 'BBRs': 20e6 }; #Mode of the baseband signal. Let's start with sinusoids
    t.init()
    t.set_transmit_power()
    print(np.std(t._Z.Value,axis=1))
    #print(t._Z.Value)
    #print(t._Z.Value.shape)
    #n=t._Z.Value/np.std(t._Z.Value,axis=1)
    print(t._Z.Value)
    #print(filt)
    #print(filt.shape)
    #tf=factor(8)
    #print(tf)
    #tf=factor(80)
    #print(tf)
    filt=t._filterlist
    for i in filt:
        w, h = sig.freqz(i)
        plt.plot(w/(2*np.pi), 20*np.log10(abs(h)))
    plt.show()


