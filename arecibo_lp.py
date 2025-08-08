import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.signal.windows as sw
import scipy.signal as ss

import arecibo_lib as al    



# python3 -m venv /home/j/arecibo
# source /home/j/arecibo/bin/activate

d=drf.DigitalRFReader("/mnt/data/juha/arecibo_snippet")
c=d.get_channels()
#print(c)
b=d.get_bounds("ch")
print("samples since 1970")
#print(b)
sample_rate = 25000000

# 10 ms ipp
ipp=10000*25

# coded long pulse start
#i0=14765+b[0] + 1068*ipp

# uncoded long pulse start
i0=14765+b[0] 

n_ipp = 10000

range_gate = 800e3 # km
range_gates=n.linspace(80e3,1400e3,num=250)
n_rg = len(range_gates)
range_gate_idx = int((2*range_gate/sc.c)*sample_rate)

range_gate_idxs = n.array(((2*range_gates/sc.c)*sample_rate),dtype=int)

pulse_length = 11100

window=sw.hann(pulse_length)
fft_length = pulse_length*2
spectrum=n.zeros(fft_length)
spectrums=n.zeros([2,fft_length,n_rg])

doppler_freq = n.fft.fftshift(n.fft.fftfreq(fft_length,d=1/sample_rate))

means=[]

for i in range(n_ipp):

    # read raw voltage for IPP i
    z=d.read_vector_c81d(i0+i*ipp,10000*25,"ch")
    # remove DC
    dc=n.mean(z[150000:190000])
#    z=z-dc
    #plt.plot(z.real)
    #plt.plot(z.imag)
    #plt.show()
    
    z_tx = n.conj(n.copy(z[0:pulse_length]))
    
    # ground clutter removal
    z[0:(5000+pulse_length)]=0
    means.append(dc)
    print(n.mean(means))
    # normalize
    z_tx=z_tx/n.max(n.abs(z_tx))
    
    coded_pulse=al.is_coded(z_tx)
    print("pulse %d coded %s"%(i,coded_pulse))
    if coded_pulse:
        continue
    
    for ri in range(n_rg):
        echo = z[(range_gate_idxs[ri]):(range_gate_idxs[ri]+pulse_length)]

        spectrums[0,:,ri] += n.abs(n.fft.fftshift(n.fft.fft(z_tx*echo,fft_length)))**2.0            
        

    if i%100 == 99:
        fidx=n.where(n.abs(doppler_freq/1e3)<100)[0]
        spectrums2=n.copy(spectrums)
#        for ri in range(n_rg):
 #           spectrums2[0,:,ri]=spectrums2[0,:,ri]-spectrums2[0,:,-1]
        plt.pcolormesh(doppler_freq[::100]/1e3,range_gates/1e3,10.0*n.log10(n.abs(spectrums2[0,::100,:].T)))
        plt.colorbar()
        plt.show()

        plt.pcolormesh(doppler_freq[fidx]/1e3,range_gates/1e3,10.0*n.log10(n.abs(spectrums2[0,fidx,:].T)))
        plt.colorbar()
        plt.show()

