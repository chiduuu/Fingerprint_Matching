# Fingerprint_Matching

A fingerprint recognition software using Discrete Wavelet Transformation (DWT). 

Fingerprint based identification is one of the most important biometric technologies which has drawn a substantial amount of attention recently. Fingerprints are believed to be unique across individuals and across fingers of same individual. These observations have led to the increased use of automatic fingerprint based recognition technique.

Minutiae based methods are most widely used generate spurious data that makes recognition of fingerprint difficult especially in noisy fingerprint. Lot of pre-processing is required like eliminate noise, reduce the redundant data and validating the minutiae points etc. All these issues add to the tie complexity thus making un-suitable for real time applications.

This project demonstrates an alternative method that is based on image transforms, that generates transform coefficients (A set of 36 values) and these represent the content of fingerprint image. Finally fingerprints are compared in terms of transform coefficients using distance measure between feature vectors. This technique is also useful in devices which have small amount of memory like smart cards etc.,
