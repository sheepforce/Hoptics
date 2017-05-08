# HOptics -- *calculation of optical data from thin film spectra*
-------------

## Build 
Hoptics is built and installed by a Makefile and typing

> make all && make install

### Dependencies
In order to build Hoptics you will need the Glasgow Haskell Compiler `ghc`, the GNU Scientific Library, BLAS and some Haskell packages. They can be installed by

> aptitude install libgsl-dev libopenblas-dev ghc cabal-install
> cabal install text hmatrix hmatrix-gsl either-unwrap attoparsec

-------------
## Usage
Call HOptics with the file name of the spectrum as argument `hoptics /path/to/spectrum.dat`. You may now choose how to proceed.

 - HOptics
	 - derive index of refraction from spectrum
		 - calculate index of refraction by Kramers-Kronig relation
	 - mix two sets of spectra
		 - Maxwell Garnet
		 - Bruggeman


