import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.signal.windows as sw
import scipy.signal as ss

def is_coded(ztx):
    """
    quick way to determine if this is a coded long pulse
    """
    ztx=ztx/n.max(n.abs(ztx))
    thresh=len(ztx)*0.8
    tsum=n.abs(n.sum(ztx))
#    print(tsum)
    if  tsum > thresh:
        return(False)
    else:
        return(True)
