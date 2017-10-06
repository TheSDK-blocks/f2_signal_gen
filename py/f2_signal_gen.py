# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 05.10.2017 17:23
import numpy as np
import tempfile
import subprocess
import shlex
import time

import modem as mdm #Function definitions

#Classes are distinguished by the class
from refptr import *
from thesdk import *
from rtl import *

#Simple buffer template
class f2_signal_gen(thesdk):
    def __init__(self,*arg): 
        self.proplist = [ 'M', 'K', 'Rs', 'bbsigdict', 'ofdmdict' ];  #Properties that can be propagated from parent
        self.Rs = 100e6                                #Sampling frequency
        self.M=4                                       #Number of antennas
        self.K=2                                       #Number of users
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

    def run(self): #Just an alias for init to be consistent: run() always executes the core function
        self.init()


    #Methods for signal generation. Add a function and add it to init()
    #controlled with bbsigdict
    def sinusoid(self):
            length=self.bbsigdict['length']
            phi=np.transpose(np.array(np.mat(self.bbsigdict['freqs'])))*np.array(range(length))*2*2*np.pi/(self.Rs)
            out=np.sum(np.exp(1j*phi),0)/len(self.bbsigdict['freqs'])
            self._Z.Value=np.ones((self.M,1))*out #All antennas emit the same signal

    def ofdm_sinusoid(self):
            ofdmdict=self.ofdmdict
            length=np.floor(self.bbsigdict['length']/ofdmdict['framelen'])
            signalindexes=np.round(np.array(self.bbsigdict['freqs'])/self.Rs*ofdmdict['framelen']).astype(int)
            frame=np.zeros((length,ofdmdict['framelen']))
            frame[:,signalindexes]=1
            datasymbols=frame[:,ofdmdict['data_loc']]
            pilotsymbols=frame[:,ofdmdict['pilot_loc']]
            out=np.ones((self.M,1))*mdm.ofdmMod(ofdmdict,datasymbols,pilotsymbols)
            self._Z.Value=np.ones((self.M,1))*out #All antennas emit the same signal



