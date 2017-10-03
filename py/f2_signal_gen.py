# f2_signal_gen class 
# Last modification by Marko Kosunen, marko.kosunen@aalto.fi, 03.10.2017 11:42
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
        self.proplist = [ 'M' 'K' 'Rs' 'bbsigmode' ];  #Properties that can be propagated from parent
        self.Rs = 1                                    #Sampling frequency
        self.M=4                                       #Number of antennas
        self.K=2                                       #Number of users
        self.bbsigmode='sinusoid';                     #Mode of the basebans signal. Let's start with sinusoids
        self.model='py';                               #can be set externally, but is not propagated
        self._Z = refptr();
        self._classfile=__file__
        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;

    def init(self):
        pass
    def run(self):
        if mode=='sinusoid'
            phase=
            out=np.array(self.iptr_A.Value)
            if par:
                queue.put(out)
            self._Z.Value=out
        else: 
          try:
              os.remove(self._infile)
          except:
              pass

