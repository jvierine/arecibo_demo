import glob
import numpy as n
import h5py
import digital_rf as drf
import os
#
# convert old style arecibo digital rf data to modern digital_rf
# this problem only exists because arecibo data was taken with an early
# prototype of digital_rf. Not a big problem luckily!
# 
# pip install digital_rf
# pip install numpy==1.26.4
# 

def digital_rf_hdf5(dirname="/media/j/c6295cbc-61f6-48df-b91e-8c9a3d520bea/arecibo/2016.05.11/ch",outdir="/data2/arecibo/2016.05.11/ch"):

    os.system("mkdir -p %s"%(outdir))
    fl=glob.glob("%s/*/rf*.h5"%(dirname))
    fl.sort()
    h=h5py.File(fl[0],"r")
    i0=h["rf_data_index"][()][0,0]
    sr=25000000
    sample_rate_numerator=25000000
    sample_rate_denominator=1
    dtype_str="i2"
    sub_cadence_secs=3600
    file_cadence_millisecs=1000
    compression_level=0
    checksum=False
    is_complex=True
    is_continuous=True
    num_subchannels=1
    marching_periods=True
    uuid="FAKE UUID"
    vector_length=25000000
    start_global_index=i0
    
    dwo = drf.DigitalRFWriter(
        outdir,
        dtype_str,
        sub_cadence_secs,
        file_cadence_millisecs,
        start_global_index,
        sample_rate_numerator,
        sample_rate_denominator,
        uuid,
        compression_level,
        checksum,
        is_complex,
        num_subchannels,
        is_continuous,
        marching_periods,
    )
    raw_volt=h["rf_data"][()]
    dwo.rf_write(raw_volt)
    h.close()
    for f in fl:
        h=h5py.File(f,"r")
        i1=h["rf_data_index"][()][0,0]
        raw_volt=h["rf_data"][()]
        dwo.rf_write_blocks(raw_volt,n.array([i1-i0+25000000]),n.array([0]))
        h.close()

digital_rf_hdf5()
    
        
