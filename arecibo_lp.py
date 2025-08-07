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
    thresh=len(ztx)*0.8
    if n.abs(n.sum(ztx)) > thresh:
        return(False)
    else:
        return(True)

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

range_gate = 800e3 # km
range_gates=n.linspace(80e3,1000e3,num=250)
n_rg = len(range_gates)
range_gate_idx = int((2*range_gate/sc.c)*sample_rate)

range_gate_idxs = n.array(((2*range_gates/sc.c)*sample_rate),dtype=int)

pulse_length = 11100

window=sw.hann(pulse_length)
fft_length = pulse_length*2
spectrum=n.zeros(fft_length)
spectrums=n.zeros([2,fft_length,n_rg])

doppler_freq = n.fft.fftshift(n.fft.fftfreq(fft_length,d=1/sample_rate))

for i in range(n_ipp):

    # read raw voltage for IPP i
    z=d.read_vector_c81d(i0+i*ipp,10000*25,"ch")
    #plt.plot(z.real)
    #plt.plot(z.imag)
    #plt.show()
    
    z_tx = n.conj(n.copy(z[0:pulse_length]))
    # ground clutter removal
    z[0:(5000+pulse_length)]=0
    # normalize
    z_tx=z_tx/n.max(n.abs(z_tx))
    
    coded_pulse=is_coded(z_tx)
    print("pulse %d coded %s"%(i,coded_pulse))
    
    for ri in range(n_rg):
        # average 13 neighbouring range gates. 2 microsecond code bit length
        # means that we need 2 microsecond steps to get independent samples of the ionospheric scatter
        for ra in range(-6,6):
            echo = z[(range_gate_idxs[ri]+ra*50):(range_gate_idxs[ri]+ra*50+pulse_length)]

            if coded_pulse:
                spectrums[1,:,ri] += n.abs(n.fft.fftshift(n.fft.fft(z_tx*echo,fft_length)))**2.0
            else:
                spectrums[0,:,ri] += n.abs(n.fft.fftshift(n.fft.fft(z_tx*echo,fft_length)))**2.0            
        

    if i%100 == 99:
        if True:
            fidx=n.where(n.abs(doppler_freq/1e3)<100)[0]
            if False:
                plt.pcolormesh(doppler_freq[::100]/1e3,range_gates/1e3,10.0*n.log10(spectrums[0,::100,:].T))
                plt.colorbar()
                plt.show()

                plt.pcolormesh(doppler_freq[fidx]/1e3,range_gates/1e3,10.0*n.log10(spectrums[0,fidx,:].T))
                plt.colorbar()
                plt.show()

            pspec=n.copy(spectrums[1,fidx,:])
            for ri in range(n_rg):
                pspec[:,ri]=(pspec[:,ri]-n.median(pspec[:,ri]))*(range_gates[ri]/1e3)**2.0
                    
            plt.pcolormesh(doppler_freq[fidx]/1e3,range_gates/1e3,pspec[:,:].T)
            plt.colorbar()
            plt.show()
        if False:
            # plot
            plt.subplot(211)
            plt.plot(z.real)
            plt.plot(z.imag)
            plt.axvline(range_gate_idx)
            plt.axvline(range_gate_idx+pulse_length)
            plt.xlabel("Time (samples)")
            plt.ylabel("Complex voltage")

            plt.subplot(212)
            plt.plot(doppler_freq/1e3,10.0*n.log10(spectrum))
            plt.xlabel("Doppler frequency (kHz)")
            plt.ylabel("Received power (dB)")
            plt.title(r"$N_{\mathrm{ipp}}=%d$"%(i+1))
            plt.tight_layout()
            plt.show()
