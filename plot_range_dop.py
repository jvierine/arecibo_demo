import numpy as n
import h5py
import matplotlib.pyplot as plt

def plot_ambiguity(fname="clp_rda.h5",title="Coded long pulse",xlim=[-1,1]):
    h=h5py.File(fname,"r")
    print(h.keys())
    ranges=h["range"][()]
    dops=h["doppler"][()]/1e6

    plt.imshow(10.0*n.log10(h["psf"][()].T),vmax=0,vmin=-40,extent=[n.min(dops),n.max(dops),n.min(ranges),n.max(ranges)],aspect="auto")
    plt.xlim(xlim)
    cb=plt.colorbar()
    cb.set_label("dB")
    plt.xlabel("Doppler (MHz)")
    plt.ylabel("Range (km)")
    plt.title("%s Range-Doppler ambiguity"%(title))
    plt.show()
    h.close()


plot_ambiguity(fname="clp_rda.h5",title="Coded long pulse",xlim=[-1,1])
plot_ambiguity(fname="lp_rda.h5",title="Uncoded long pulse",xlim=[-0.1,0.1])
