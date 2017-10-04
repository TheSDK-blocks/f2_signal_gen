# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 04.10.2017 09:58
import numpy as np
import tempfile
import subprocess
import shlex
import time

from refptr import *
from thesdk import *
from rtl import *

#Simple buffer template
class f2_signal_gen(thesdk):
    def __init__(self,*arg): 
        self.proplist = [ 'M' 'K' 'Rs' 'bbsigdict' ];  #Properties that can be propagated from parent
        self.Rs = 100e6                                #Sampling frequency
        self.M=4                                       #Number of antennas
        self.K=2                                       #Number of users
        self.bbsigdict={ 'mode':'sinusoid', 'freqs':[11.0e6 , 13e6, 17e6 ], 'length':2**14 };                     #Mode of the baseband signal. Let's start with sinusoids
        self.model='py';                               #can be set externally, but is not propagated
        self._Z = refptr();
        self._classfile=__file__
        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;
        self.init()

    def init(self):
        if self.bbsigdict['mode']=='sinusoid':
            self.sinusoid()

    def run(self): #Just an alias for init to be consistent: run() always executes the core function
        self.init()

    def sinusoid(self):
            length=self.bbsigdict['length']
            phi=np.transpose(np.array(np.mat(self.bbsigdict['freqs'])))*np.array(range(length))*2*2*np.pi/(self.Rs)
            out=np.sum(np.exp(1j*phi),0)/len(self.bbsigdict['freqs'])
            self._Z.Value=np.ones((self.M,1))*out #All antennas emit the same signal


