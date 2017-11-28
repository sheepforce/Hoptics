# HOptics -- *calculation of optical data from thin film spectra*
## Build
Hoptics is built and installed by Stack

		git clone https://github.com/sheepforce/Hoptics.git
		cd Hoptics
		stack setup
		stack install


### Dependencies
In order to build Hoptics you will need the Glasgow Haskell Compiler `ghc`, the GNU Scientific Library, BLAS and some Haskell packages. They can be installed by

	aptitude install libgsl-dev libopenblas-dev ghc cabal-install
	cabal install text hmatrix hmatrix-gsl either-unwrap attoparsec


## Usage
HOptics can be used to calculate the index of refraction from thin film transmission spectra or to get an averaged index of refraction for inclusion/matrix system, where the two components has been measured separately. Thin film spectra have to be a two column ascii file with whitespace as column separators.


### Index of Refraction
For calculating the real and imaginary part of the index of refraction call hoptics with the transmission-spectrum as the argument `hoptics /path/to/spectrum.dat`. You will see the main menu.

	Hoptics
	∟ Main Menu

	(1)  derive index of refraction from spectrum
	(2)  mix two sets of indices of refraction

press `1` to get to the menu for analysing the spectra. You will see the current settings for calculation, which can be changed by entering the corresponding number an [enter].

	HOptics
	∟ Main Menu
	  ∟ analyse spectrum

	(-1) return to main menu
	(0)  start computation
	(1)  path to spectrum                                  "spectrum.dat"
	(2)  thickness of the sample [nm]                      100.0
	(3)  spectral range for the calculation                (45000.0,590000.0)
	(4)  security distance around poles and boundaries     500.0
	(5)  seed value for real part of index of refraction   1.0
	(6)  dimension on x axis                               Wavenumber
	(7)  integration method for Kramers Kronig             Akima

 - (1) change the file name of the spectrum or give it here, if it was not given a argument
 - (2) thickness of the sample in nano metres
 - (3) for the integration in the Kramers-Kronig relation these are the integration bounds in m⁻¹
 - (4) which range to ommit (in m⁻¹) in the integration around the poles, if you get GSL errors, try increasing this value
 - (5) seed value for the real part of the index of refraction
 - (6) is the dimension on the x-axis a wavelength (assumed to be given in nm) or a wavenumber (assumed to be given in cm⁻¹)
 - (7) integration method to use for Kramers-Kronig relation. `Akima` is very good choice and fast. `Naive` is slow but works for problematic cases. Accuracy of `Naive` is only accurate enough for well resolved spectra

After the settings are made, you can hit `0`  for starting the calculation. Be patient, the integration has to be done many times. Afterwards you will find the files `${prefix}_trans.dat ${prefix}_alpha.dat ${prefix}_k.dat ${prefix}_n0.dat` in your directory. The contain the transmission spectrum in m⁻¹, the extinction coefficient in m⁻¹, the imaginary part of the index of refraction and the real part of the index of refraction.

### Mixing of Indices of Refraction
From the main menu press `2` to enter the menu for getting the index of refraction of a mixture of two components with know indices of refraction. Again, you can change the settings by pressing the according number

	Hoptics
	∟ Main Menu
	  ∟ mix spectra

	(-1) return to main menu
	(0)  start computation
	(1)  prefix of spectrum of component 1 (inclusion)       "spectrum1"
	(2)  prefix of spectrum of component 2 (matrix)          "spectrum2"
	(3)  volume fraction of component 1                      0.5
	(4)  magnetic permittivity of the mixture                1.0

 - (1) prefix name of the spectra of the inclusion. If you have the files `sample1_n0.dat sample1_k.dat` your prefix is `sample1`
 - (2) prefix name of the spectra of the matrix
 - (3) the volume fraction of the inclusion. Should typically be <0.5 as the justification for value >0.5 is difficult
 - (4) relative magnettic permittivity of the mixture, assuming a diamagnetic sample

After pressing `0` the computation starts and you will get two files labeled `${prefix1}+${prefix2}_k_MaxwellGarnet.dat ${prefix1}+${prefix2}_n0_MaxwellGarnet.dat` containing the resulting imaginary and real part of the index of refraction for the mixture.
