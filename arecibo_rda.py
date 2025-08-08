import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.signal.windows as sw
import scipy.signal as ss

import range_doppler_ambiguity as rda

import arecibo_lib as al        
    
from mpi4py import MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

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

n_ipp=int(n.floor( (b[1]-i0)/ipp))


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

r_clp=rda.range_dop_amb(pulse_length,2,0,sample_rate,delta_rg=25)
r_lp=rda.range_dop_amb(pulse_length,2,0,sample_rate,delta_rg=25)



for i in range(n_ipp):
    # read raw voltage for IPP i
    z=d.read_vector_c81d(i0+i*ipp,10000*25,"ch")
    z_tx = n.copy(z[0:pulse_length])
    coded_pulse=al.is_coded(z_tx)
    print("%d/%d coded %d"%(i,n_ipp,coded_pulse))
    print(coded_pulse)
    plt.plot(z_tx.real)
    plt.plot(z_tx.imag)
    plt.show()
    if coded_pulse:
        if rank==1:
            r_clp.add_tx(z_tx)
            r_clp.save("clp_rda.h5")
    else:
        if rank==0:
            r_lp.add_tx(z_tx)
            r_lp.save("lp_rda.h5")        
