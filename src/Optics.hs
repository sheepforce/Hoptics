module Optics where

import Numeric.GSL

-- speed of light
c = 299792458.0e0 :: Double

-- wavelength on x to frequency on x, assuming meter as input
wavelength2frequency :: [(Double,Double)] -> [(Double,Double)]
wavelength2frequency spectrum = map conversion spectrum
    where
        conversion (lambda,y) = (c/lambda,y)

wavelength2wavenumber :: [(Double,Double)] -> [(Double,Double)]
wavelength2wavenumber spectrum = map conversion spectrum
    where
        conversion (lambda,y) = (1.0 / lambda,y)

-- wavenumber on x to frequency on x, assuming inverse meter as input
wavenumber2frequency :: [(Double,Double)] -> [(Double,Double)]
wavenumber2frequency spectrum = map conversion spectrum
    where
        conversion (nu,y) = (nu*c,y)

wavenumber2wavelength :: [(Double,Double)] -> [(Double,Double)]
wavenumber2wavelength spectrum = map conversion spectrum
    where
        conversion (nu,y) = (1.0 / nu,y)

frequency2wavenumber :: [(Double,Double)] -> [(Double,Double)]
frequency2wavenumber spectrum = map conversion spectrum
    where
        conversion (freq,y) = (freq / c,y)

frequency2wavelength :: [(Double,Double)] -> [(Double,Double)]
frequency2wavelength spectrum = map conversion spectrum
    where
        conversion (freq,y) = (c / freq,y)

-- absorption coefficient alpha from transmission spectrum
absorption_coeff :: Double -> [(Double,Double)] -> [(Double,Double)]
absorption_coeff thickness trans_spectrum = map get_alpha trans_spectrum
    where
        get_alpha (nu,trans) = (nu,1.0/thickness * log (1.0 / trans))

-- imaginary part of the index of refraction calculated from the alpha spectrum
indOfRef_k :: [(Double,Double)] -> [(Double,Double)]
indOfRef_k alpha_spectrum = map get_k alpha_spectrum
    where
        get_k (nu,alpha) = (nu,alpha / (4.0 * pi * nu))


-- calculate the real part of the index of refraction by Kramers Kronig relation from
-- the alpha spectra. Needs a seed value
indOfRef_n :: (Double,Double) -> Double -> Double -> [(Double,Double)] -> [(Double,Double)]
indOfRef_n (low_limit,high_limit) n_seed safety_distance alpha_spectrum = [(nu_indOfRef, n_seed + 1.0 / (2.0 * pi**2.0) * (alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum)) 
                                                          | nu_indOfRef <- (map fst alpha_spectrum), nu_indOfRef >= low_limit, nu_indOfRef <= high_limit ]
    where
        -- integrate the quotient term for a specific index of refraction with the corresponding
        -- series of point of the quotient
        alpha_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
        --alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,high_limit)
        alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = low_integral + high_integral        
            where
                alpha_quotient_spectrum = alpha_quotients nu_indOfRef alpha_spectrum
                
                low_integral = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,nu_indOfRef - 0.0)
                high_integral = evaluateIntegral Akima alpha_quotient_spectrum (nu_indOfRef + 0.0,high_limit)
                
                -- for a specific wavenumber generate the series of points, that need to 
                -- be integrated. This gives the term
                --
                -- alpha(nu') / (nu'² - nu²)
                --
                -- as a series the full spectral range of the alpha spectrum (nu') at the specific
                -- wavenumber nu
                alpha_quotients :: Double -> [(Double,Double)] -> [(Double,Double)]
                --alpha_quotients nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (snd a) <= 10.0]
                alpha_quotients nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (fst a) <= (nu_indOfRef - safety_distance)] ++ [a | a <- alpha_quotients_raw, (fst a) >= (nu_indOfRef + safety_distance)]
                    where
                        alpha_quotients_raw = map get_quotient alpha_spectrum
                        get_quotient (nu,alpha) = (nu,alpha / (nu**2.0 - nu_indOfRef**2.0))


{-
-- works
alpha_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,high_limit)
    where
        alpha_quotient_spectrum = alpha_quotients nu_indOfRef alpha_spectrum
-}
-- works
{-
alpha_quotients :: Double -> [(Double,Double)] -> [(Double,Double)]
alpha_quotients nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (snd a) <= 10.0]
    where
        alpha_quotients_raw = map get_quotient alpha_spectrum
        get_quotient (nu,alpha) = (nu,alpha / (nu**2.0 - nu_indOfRef**2.0))

-- take a alpha spectrum and create a list of wavenumbers with corresponding quotient spectra
alpha_quotients_spectrum :: [(Double,Double)] -> [(Double,[(Double,Double)])]
alpha_quotients_spectrum alpha_spectrum = [(nu,alpha_quotients nu alpha_spectrum) | nu <- (map fst alpha_spectrum)]

alpha_quot_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
alpha_quot_integral (low_limit,high_limit) nu_indOfRef alpha_quotient_spectrum = low_integral + high_integral
    where
        low_integral = evaluateIntegral Akima alpha_quotient_spectrum (low_limit, nu_indOfRef - 10.0)
        high_integral = evaluateIntegral Akima alpha_quotient_spectrum (nu_indOfRef + 10.0, high_limit)

--indOfRef_n :: (Double,Double) -> Double -> [(Double,[(Double,Double)])] -> [(Double,Double)]
indOfRef_n :: (Double,Double) -> Double -> [(Double,Double)] -> [(Double,Double)]
indOfRef_n (low_limit,high_limit) n_seed alpha_spectrum = zip wavenumbers indicesOfRef
    where
        wavenumbers = map fst alpha_quotients_spectra
        indicesOfRef = [n_seed + 1.0 / (2.0 * pi^2) * (alpha_quot_integral (low_limit,high_limit) (wavenumbers!!spec_ind) ((map snd alpha_quotients_spectra)!!spec_ind)) | spec_ind <- [0 .. ((length wavenumbers) - 1)],  (wavenumbers!!spec_ind) >= low_limit, (wavenumbers!!spec_ind) <= high_limit]
        
        alpha_quotients_spectra = alpha_quotients_spectrum alpha_spectrum
        
        alpha_quotients_spectrum :: [(Double,Double)] -> [(Double,[(Double,Double)])]
        alpha_quotients_spectrum alpha_spectrum = [(nu,alpha_quotients nu alpha_spectrum) | nu <- (map fst alpha_spectrum)]

        alpha_quotients :: Double -> [(Double,Double)] -> [(Double,Double)]
        alpha_quotients nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (snd a) <= 10.0]
            where
                alpha_quotients_raw = map get_quotient alpha_spectrum
                get_quotient (nu,alpha) = (nu,alpha / (nu**2.0 - nu_indOfRef**2.0))
        
        alpha_quot_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
        alpha_quot_integral (low_limit,high_limit) nu_indOfRef alpha_quotient_spectrum = low_integral + high_integral
            where
                low_integral = evaluateIntegral Akima alpha_quotient_spectrum (low_limit, nu_indOfRef - 1000.0)
                high_integral = evaluateIntegral Akima alpha_quotient_spectrum (nu_indOfRef + 1000.0, high_limit)
-}
