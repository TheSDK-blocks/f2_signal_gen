# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 20.10.2017 13:35
import sys
sys.path.append ('/home/projects/fader/TheSDK/Entities/refptr/py')
sys.path.append ('/home/projects/fader/TheSDK/Entities/thesdk/py')
sys.path.append ('/home/projects/fader/TheSDK/Entities/modem/py')
import numpy as np
import tempfile
import subprocess
import shlex
import time

import modem as mdm #Function definitions

#Classes are distinguished by the class
from refptr import *
from thesdk import *
#from rtl import *

#Simple buffer template
class f2_signal_gen(thesdk):
    def __init__(self,*arg): 
        self.proplist = [ 'Txantennas', 'Txpower', 'Users', 'Rs', 'bbsigdict', 'ofdmdict' ];  #Properties that can be propagated from parent
        self.Rs = 100e6                         #Sampling frequency
        self.Txantennas=4                       #Number of transmitting antennas
        self.Txpower=30                         #Output power per antenna in dBm
        self.Users=2                            #Number of users
        self.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14 }; #Mode of the baseband signal. Let's start with sinusoids
        self.ofdmdict={ 'framelen':64,'data_loc': np.r_[1:11+1, 13:25+1, 
        27:39+1, 41:53+1, 55:64+1]-1, 'pilot_loc' : np.r_[-21, -7, 7, 21] + 32, 'CPlen':16}
        self.model='py';                               #can be set externally, but is not propagated
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

            self._Z.Value=out 

    def ofdm_sinusoid(self):
            ofdmdict=self.ofdmdict
            length=np.floor(self.bbsigdict['length']/ofdmdict['framelen'])
            signalindexes=np.round(np.array(self.bbsigdict['freqs'])/self.Rs*ofdmdict['framelen']).astype(int)
            frame=np.zeros((length,ofdmdict['framelen']))
            frame[:,signalindexes]=1
            datasymbols=frame[:,ofdmdict['data_loc']]
            pilotsymbols=frame[:,ofdmdict['pilot_loc']]
            out=np.ones((self.Txantennas,1))*mdm.ofdmMod(ofdmdict,datasymbols,pilotsymbols).T
            
            out=np.zeros(self.Users,usersig.shape[0],usersig,shape[1])
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
            #The length is approx this many frames
            frames=np.floor(length/(framelen+CPlen))
            bitspersymbol=np.log2(QAM).astype(int)
            
            #generate random bitstreams per antenna
            bitstream=np.random.randint(2,size=(self.Txantennas,frames*bitspersymbol*framelen))
            #Init the qam signal, frame and out
            qamsignal=np.zeros((self.Txantennas,frames*framelen),dtype='complex')
            frame=np.zeros((frames,framelen))
            usersig=np.zeros((self.Txantennas,frames*(framelen+CPlen)),dtype='complex')
            for i in range(bitstream.shape[0]):
                wordstream, qamsignal[i]= mdm.qamModulateBitStream(bitstream[i], QAM)
                qamsignal[i]=qamsignal[i].reshape((1,qamsignal.shape[1]))
                frame= qamsignal[i].reshape((-1,framelen))
                datasymbols=frame[:,ofdmdict['data_loc']]
                pilotsymbols=frame[:,ofdmdict['pilot_loc']]
                usersig[i]=mdm.ofdmMod(ofdmdict,datasymbols,pilotsymbols)

            out=np.zeros(self.Users,usersig.shape[0],usersig,shape[1])
            for i in range(self.Users):
                out[i,:,:]=usersig

            self._Z.Value=out/np.st 

    def set_transmit_power(self):
         for user in range(self._Z.Value.shape[0]):
             for antenna in range(self._Z.Value.shape[2]):
                 self._Z.Value[user,:,antenna]=self._Z.Value[user,:,antenna]/np.std(self._Z.Value[user,:,antenna])*np.sqrt(10**(self.Txpower/10)*1e-3*50) #Rms voltage normalized to 50 ohm
                 #print(np.std(self._Z.Value[user,:,antenna]))


if __name__=="__main__":

    import numpy as np
    from  f2_signal_gen import *
    t=f2_signal_gen()
    t.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14 }; #Mode of the baseband signal. Let's start with sinusoids
    t.init()
    t.set_transmit_power()
    print(np.std(t._Z.Value,axis=1))
    #print(t._Z.Value)
    #print(t._Z.Value.shape)
    #n=t._Z.Value/np.std(t._Z.Value,axis=1)
    print(t._Z.Value)


