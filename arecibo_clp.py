import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
# python3 -m venv /home/j/arecibo
# source /home/j/arecibo/bin/activate

#drf.recreate_properties_file("/media/j/6a5b1924-878b-4ba6-84eb-61abe98217b2/arecibo/2016.05.13/ch")

d=drf.DigitalRFReader("/data2/arecibo/2016.05.11")
c=d.get_channels()
print(c)
b=d.get_bounds("ch")
print("samples since 1970")
print(b)

i0=14765+b[0]

n_ipp = 1000
# 10 ms ipp
ipp=10000*25
for i in range(n_ipp):
    
    z=d.read_vector_c81d(i0+i*ipp*100,10000*25,"ch")
    plt.plot(z.real)
    plt.plot(z.imag)
    plt.xlabel("Time (samples)")
    plt.ylabel("Complex voltage")
    plt.show()
