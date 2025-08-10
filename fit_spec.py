import h5py
import matplotlib.pyplot as plt
import numpy as n
import scipy.optimize as so
import il_interp as il

ils=il.ilint(fname="ion_line_interpolate_32_16_430_00.h5")

hrd=h5py.File("lp_rda.h5","r")
psf=hrd["psf"][()]
dop=hrd["doppler"][()]
hrd.close()

rpsf=n.sum(psf,axis=1)
rpsf=rpsf/n.max(rpsf)
fidx=n.where(n.abs(dop)<100e3)[0]
rpsf=rpsf[fidx]
#plt.plot(freq,rpsf)
#plt.show()

h=h5py.File("isr_spec_meas.h5","r")
freq=h["doppler_freq"][()]
spec=h["spec"][()]

S=ils.getspec(ne=n.array([1e11]),
              te=n.array([2000]),
              ti=n.array([1500]),
              mol_frac=n.array([0]),
              vi=n.array([0]),
              acf=False,
              interpf=True
              )


plt.plot(freq,spec/n.max(spec))
plt.plot(freq,rpsf)
model=S[0](freq)
plt.plot(freq,model/n.max(model))
plt.show()
h.close()
