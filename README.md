# A small library for chopper calculations in C, intended for McStas

McStas can be used to simulate time-of-flight neutron scattering instrument,
but was originally develop for constant wavelength instruments.

This small library is intended to provide extra utility to McStas
time-of-flight instruments.

To start, the only utility is intended for use by spectrometers 
which have a train of choppers in their primary spectrometer to define
a wavelength band that can reach the sample.
For direct-geometry spectrometers the band(s) are very sharp,
and for indirect-geometry spectrometers the band(s) are very broad.
In either case, the possible wavelengths that can pass through the
primary spectrometer is reduced by the chopper train.
If the chopper parameters are known, as they must be for the simulation
to progress, they can be used to limit which wavelengths McStas simulates
to only those which *could* make it to the sample position.
The first utility performs the necessary chopper acceptance intersection
calculations and identifies the *envelope* of possible wavelengths.

Embedded in the library code is a utility which returns the *list* of
possible wavelength ranges.
This list could be used as input to a new source which selects
from multiple wavelength bands; or a semi-automatic `GROUP` of 
sources with a pre-source random selection between the bands/sources.
The latter solution forces registration of a new `USERVAR` in 
the particle structure, and a possibly-large group size;
both of which are undesirable.
