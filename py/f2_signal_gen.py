# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 03.10.2017 16:44
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
        self.bbsigdict={ 'mode':'sinusoid', 'freqs':[1e6 , 10e6, 15e6 ], 'lenght':2**14 };                     #Mode of the baseband signal. Let's start with sinusoids
        self.model='py';                               #can be set externally, but is not propagated
        self._Z = refptr();
        self._classfile=__file__
        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;
        self.init()

    def init(self):
        if self.bbigdict['mode']=='sinusoid'
            length=self.bbsigdict['length']
            phi=(1/np.transpose(np.array(np.mat(bbsigdict['freqs']))))*np.array(range(length))*np.pi/(self.Rs)
            out=np.exp(1j*phi)
            self._Z.Value=out

    def run(self):
        self.init()

