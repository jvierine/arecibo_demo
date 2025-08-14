import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.signal.windows as sw
import scipy.signal as ss
import h5py

import arecibo_lib as al
import scipy.fft as sfft

def idx_array(r0=100,r1=600,dr=1,txlen=11100):
    drg=dr*150
    n_rg=int((r1-r0)*1e3/drg)
    drg_samp=25*dr
    rg0=(r0*1e3/150)*25
    rgs=n.arange(n_rg)*drg_samp + rg0
    idx_arr=n.zeros([n_rg,txlen],dtype=int)
    txidx=n.arange(txlen,dtype=int)
    
    for ri in range(n_rg):
        idx_arr[ri,:] = txidx + rgs[ri]
    range_gates=rgs*0.15/25
    return(idx_arr,range_gates)

# python3 -m venv /home/j/arecibo
# source /home/j/arecibo/bin/activate

d=drf.DigitalRFReader("/mnt/data/juha/arecibo/2016.05.11")
c=d.get_channels()
#print(c)
b=d.get_bounds("ch")
print("samples since 1970")
#print(b)
sample_rate = 25000000

# 10 ms ipp
ipp=10000*25

# coded long pulse start
i0=14765+b[0] + 1068*ipp

# uncoded long pulse start
#i0=14765+b[0] 

n_ipp = 10000

pulse_length = 11100
fft_length = pulse_length*2
doppler_freq = n.fft.fftshift(n.fft.fftfreq(fft_length,d=1/sample_rate))
idxmat,range_gates=idx_array(txlen=pulse_length)

n_rg=len(range_gates)
S=n.zeros([n_rg,fft_length],dtype=n.float32)

for i in range(n_ipp):
    print("%d/%d"%(i,n_ipp))
    # read raw voltage for IPP i
    z=d.read_vector_c81d(i0+i*ipp,10000*25,"ch")
    # conjugated transmit pulse for decoding
    z_tx = n.conj(n.copy(z[0:pulse_length]))
    
    # decode and build matrix with echoes of each range gate on each row
    Z=z[idxmat]*z_tx[None,:]
    S+=n.abs(n.fft.fftshift(sfft.fft(Z,fft_length,axis=1,workers=48),axes=1))**2.0
    
    if i%100 == 99:
        # median filter is approximately the self noise and background noise
        Sl = n.copy(S)
        for ri in range(n_rg):
            print("%d/%d"%(ri,n_rg))
            Sl[ri,:]=S[ri,:]-ss.medfilt(S[ri,:],31)
        vlow,vhigh=n.percentile(Sl,[5,99])
        #        plt.imshow(10*n.log10(S[::-1,:]),aspect="auto")
        plt.imshow(Sl[::-1,:],aspect="auto",vmin=0,vmax=vhigh,extent=[n.min(doppler_freq)/1e6,n.max(doppler_freq)/1e6,n.min(range_gates),n.max(range_gates)])#10*n.log10(S[::-1,:]),aspect="auto")
        plt.xlabel("Doppler shift (MHz)")
        plt.ylabel("Range (km)")
        plt.colorbar() 
        plt.show()        

        ho=h5py.File("isr_spec_clp_meas.h5","w")
        ho["spec"]=S
        ho["range_gates"]=range_gates
        ho["doppler_freq"]=doppler_freq
        ho.close()
