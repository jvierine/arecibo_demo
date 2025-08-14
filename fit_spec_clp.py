import h5py
import matplotlib.pyplot as plt
import numpy as n
import scipy.optimize as so
import il_interp as il

def forward_model(magic_const,te,ti,vi,psf,freqs,ils):
    """
    simple ambiguity function that assumes zero range gradient in plasma-parameters, 
    only taking into account doppler broadening
    """
    n_freq=len(freqs)
    specf=ils.getspec(ne=n.array([1e11]),
                      te=n.array([te]),ti=n.array([ti]),vi=n.array([vi]),acf=False,interpf=True)[0]
    real_spec=specf(freqs)
    conv_spec=n.zeros(n_freq)

    zidx=n.where(freqs==0)[0][0]
    # this can be sped up with fft!
    for i in range(n_freq):
        for j in range(n_freq):
            if ((i-j+zidx) > 0) and ((i-j+zidx)<n_freq):
                conv_spec[i]+=psf[i-j+zidx]*real_spec[j]
    return(magic_const*conv_spec)

ils=il.ilint(fname="ion_line_interpolate_32_16_430_00.h5")

hrd=h5py.File("clp_rda.h5","r")
psf=hrd["psf"][()]
dop=hrd["doppler"][()]
hrd.close()

rpsf=n.sum(psf,axis=1)

fidx=n.where(n.abs(dop)<100e3)[0]
rpsf=rpsf[fidx]
rpsf=rpsf/n.sum(rpsf)

#        ho=h5py.File("isr_spec_clp_meas.h5","w")
h=h5py.File("isr_spec_clp_meas.h5","r")
freq=h["doppler_freq"][()][fidx]
print(len(freq))
print(len(dop))
spec=h["spec"][()]
range_gates=h["range_gates"][()]




S=ils.getspec(ne=n.array([1e11]),
              te=n.array([2000]),
              ti=n.array([1500]),
              mol_frac=n.array([0]),
              vi=n.array([0]),
              acf=False,
              interpf=True
              )

ri=n.argmin(n.abs(range_gates-300))
print(ri)
print(range_gates[ri])
spec=n.sum(spec[(ri-20):(ri+20),fidx],axis=0)


meas=spec/n.max(spec)
plt.plot(freq,meas)
plt.plot(freq,forward_model(25,4200,1400,0,rpsf,freq,ils))
plt.show()

fidx2=n.where(n.abs(freq)<15e3)[0]
def ss(x):
    #    model=forward_model(1.5,4200,1400,0,rpsf,freq,ils)
    model=forward_model(x[0],x[1],x[2],x[3],rpsf,freq,ils)
    
    return(n.sum(n.abs(model[fidx2]-meas[fidx2])**2))

x=so.fmin(ss,[25,4200,2000,0])

print(x)

model=forward_model(x[0],x[1],x[2],x[3],rpsf,freq,ils)




plt.plot(freq,meas,"x")
plt.plot(freq,model)
plt.title(r"$T_e=%1.0f$ K $T_i=%1.0f$ K $v_i=%1.0f$ m/s"%(x[1],x[2],x[3]))
plt.xlim([-20e3,20e3])
plt.xlabel("Doppler (kHz)")
plt.show()
#model=S[0](freq)
#plt.plot(freq,model/n.max(model))
#plt.show()
h.close()
