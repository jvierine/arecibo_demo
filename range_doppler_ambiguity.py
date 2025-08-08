import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.signal.windows as sw
import scipy.signal as ss
import h5py

class range_dop_amb:
    def __init__(self,txlen,freq_oversampling_factor,missing_samples=0,sample_rate=25e6,delta_rg=25):
        self.txlen=txlen
        self.freq_oversampling_factor=freq_oversampling_factor
        self.missing_samples=missing_samples
        self.sample_rate=sample_rate
        self.delta_rg=delta_rg
        self.echo=n.zeros(3*txlen+1,dtype=n.complex64)

        self.fftlen=txlen*freq_oversampling_factor
        
        self.rgs=n.arange(-txlen,txlen,delta_rg,dtype=int)
        self.n_rg=len(self.rgs)
        self.echo = n.zeros(3*txlen+1,dtype=n.complex64)

        self.P=n.zeros([self.fftlen,self.n_rg])
        self.n_tx=0
        #echo[0:missing_samples]=0
        
        self.fvec=n.fft.fftshift(n.fft.fftfreq(self.fftlen,d=1/self.sample_rate))
        self.rvec=(self.rgs/self.sample_rate)*sc.c/2/1e3
        self.P=n.zeros([self.fftlen,self.n_rg],dtype=n.float64)
        
    def add_tx(self,z_tx):
        """
        range-Doppler ambiguity function
        """
        z_tx=z_tx/n.sqrt(n.sum(n.abs(z_tx)**2.0))
        self.echo[(self.txlen-1):(self.txlen+self.txlen-1)]=z_tx
    
        for i in range(len(self.rgs)):
            self.P[:,i]+=n.abs(n.fft.fftshift(n.fft.fft(self.echo[(self.rgs[i]+self.txlen):(self.rgs[i]+self.txlen+self.txlen)]*n.conj(z_tx),self.fftlen)))**2.0
        self.n_tx+=1


    def get_psf(self):
        return(self.P/self.n_tx)

    def save(self,fname="clp.h5"):
        ho=h5py.File(fname,"w")
        ho["psf"]=n.array(self.P/self.n_tx,dtype=n.float32)
        ho["doppler"]=self.fvec
        ho["range"]=self.rvec
        ho.close()
        
#    fidx=n.where(n.abs(fvec)<100e3)[0]

#    if plot:
 #       plt.pcolormesh(fvec[fidx]/1e3,rvec,10.0*n.log10(P[fidx,:].T),vmax=0,vmin=-30)
  #      cb=plt.colorbar()
   #     cb.set_label("dB")
    #    plt.xlabel("Doppler shift (kHz)")
     #   plt.ylabel("Range offset (km)")
      #  plt.show()

#    return(P[fidx,:])    
