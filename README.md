# Ionospheric Collective Thomson Scatter Demo

The demo comes with 20 seconds of Arecibo Observatory 430 MHz raw voltage data sampled at 25 MHz sample-rate. 10 seconds of uncoded long pulse and 10 seconds of coded long pulse data is included. This data allows observing a plasma-line profile as well as the ion-line profile. The coded long pulse provides high range and Doppler resolution, and the uncoded long pulse provides high SNR for top-side observations.

In order to follow the demo, install digital_rf

> pip install digital_rf

Download 20 seconds of 25 MHz Arecibo raw voltage

Navigate to https://zenodo.org/records/16780141 and download arecibo_snippet.tar.gz

Untar the data:

> tar xvfz arecibo_snippter.tar.gz

Test that you can read the data:

> ipython3 
> import digital_rf as drf
# path to data is to directory where you untarred arecibo_snippet..tar.gz 
> d=drf.DigitalRFReader("


Vierinen, J. (2025). 20 seconds of Arecibo 430 MHz carriage house ionospheric radar raw voltage data sampled at 25 MHz sample-rate [Data set]. Zenodo. https://doi.org/10.5281/zenodo.16780141
