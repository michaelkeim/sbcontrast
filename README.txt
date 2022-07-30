Authors
-------
Michael Alan Keim, Pieter van Dokkum, and Jiaxuan Li.

Email
-----
michael [dot] keim [at] yale [dot] edu

Description
-----------
The sbcontrast code is a method to obtain surface brightness limits that
accurately reflect an images depth at a given spatial scale. Rather than
relying on the naive Poisson expectation, this code will estimate limitations
presented by large scale variations. A full description of this method is given
in Keim et al. 2022 (https://arxiv.org/pdf/2109.09778.pdf).

Usage
-----
sbcontrast may be run from the command line. Enter 'sbcontrast -h' for help. It
requires only an image fits file to run, however users should also specify a
fits file containing masks for genuine sources (via -masks), otherwise no masks
will be applied and the reported limit will not be as deep as the true value.
Note that if the header does not contain the pixel size and photometric
zeropoint as expected, the user will need to provide these (via -pix and -zp).
Moreover, sbcontrast defaults to a 1sigma, 60 arcsec scale, though other limits
may also be specified (via -N and -s). Additional commands (-binmap and
-fluctmap) to save the binned image map and the final fluctuation map used to
determine the limit are also available.

sbcontrast may also be run within python scripts, using numpy arrays rather
than fits files. Simply import the sblimit function from sbcontrast.

Note that the specified scale may be rounded so that the binning factors are
integers. Note also that the accuracy of this limit is limited by data
reduction - for scales exceeding that of background subtraction, genuine
features above the calculated limit may have been removed. Moreover, at the
single pixel scale sbcontrast may diverge from the true rms due to
correlations introduced by re-sampling.

Citation
--------
We ask that users utilizing this tool please cite Keim et al. 2022
(https://arxiv.org/pdf/2109.09778.pdf).
