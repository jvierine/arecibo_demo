# Ionospheric Collective Thomson Scatter Demo

The demo comes with 20 seconds of Arecibo Observatory 430 MHz raw voltage data sampled at 25 MHz sample-rate. 10 seconds of uncoded long pulse and 10 seconds of coded long pulse data is included. This data allows observing a plasma-line profile as well as the ion-line profile. The coded long pulse provides high range and Doppler resolution, and the uncoded long pulse provides high SNR for top-side observations.

In order to follow the demo, install digital_rf

> pip install digital_rf

Download 20 seconds of 25 MHz Arecibo raw voltage

Navigate to https://zenodo.org/records/16780141 and download arecibo_snippet.tar.gz

Untar the data:

> tar xvfz arecibo_snippter.tar.gz

Test that you can read the data:

<code>
> ipython3 
> import digital_rf as drf
> # path to data is to directory where you untarred arecibo_snippet.tar.gz
> d=drf.DigitalRFReader("./arecibo_snippet")
> print(d.get_channels())
["ch"]
> b=d.get_bounds("ch")
> # print first and last sample index (samples since 1970)
> print(b)
(36574175850000000, 36574176499999999)
> # read one 10 ms ipp of daa (25 MHz sample-rate) 
> z=d.read_vector_c81d(b[0],250000,"ch")
> # plot it
> import matplotlib.pyplot as plt
> plt.plot(z.real)
> plt.plot(z.imag)
> plt.show()
</code>
You should get something like this:

<img width="584" height="421" alt="Screenshot 2025-08-13 at 10 27 51" src="https://github.com/user-attachments/assets/4854a6bd-ca5e-4af2-adee-5c69fb8f52e0" />


Vierinen, J. (2025). 20 seconds of Arecibo 430 MHz carriage house ionospheric radar raw voltage data sampled at 25 MHz sample-rate [Data set]. Zenodo. https://doi.org/10.5281/zenodo.16780141
