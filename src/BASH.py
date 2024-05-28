import os
import subprocess as sub

def BASH(__in, __out=False, __overwrite=False):

    if type(__in) != list: __in = __in.split()

    if __out:
        if __overwrite: sub.call(__in, stdout=open(__out, "a"))
        else: sub.call(__in, stdout=open(__out, "w"))
    else:
        sub.call(__in)
