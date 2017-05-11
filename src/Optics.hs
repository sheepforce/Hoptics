module Optics where

import Numeric.GSL
--import Data.Complex

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

nkSep2complexIndOfRef :: [(Double,Double)] -> [(Double,Double)] -> [(Double,Complex Double)]
nkSep2complexIndOfRef n0_spectrum k_spectrum = zip xVals indOfRef
    where
        xVals = map fst n0_spectrum
        nVals = map snd n0_spectrum
        kVals = [evaluate Akima k_spectrum points | points <- xVals]
        indOfRef = zipWith (:+) nVals kVals

n2epsilon :: Double -> [(Double,Complex Double)] -> [(Double,Complex Double)]
n2epsilon magnetic_permittivity n_spectrum = zip xVals epsilonVals
    where
        xVals = map fst n_spectrum
        nVals = map snd n_spectrum
        n_length = length n_spectrum
        epsilonVals = [((nVals !! ind)^2) `commul` (1 / magnetic_permittivity) | ind <- [0 .. (n_length - 1)]]

epsilon2n :: Double -> [(Double,Complex Double)] -> [(Double,Complex Double)]
epsilon2n magnetic_permittivity epsilon_spectrum = zip xVals nVals
    where
        xVals = map fst epsilon_spectrum
        epsilonVals = map snd epsilon_spectrum
        epsilon_length = length epsilon_spectrum
        nVals = [sqrt ((epsilonVals !! ind)) `commul` magnetic_permittivity | ind <- [0 .. (epsilon_length - 1)]]

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
indOfRef_n :: InterpolationMethod -> (Double,Double) -> Double -> Double -> [(Double,Double)] -> [(Double,Double)]
indOfRef_n int_method (low_limit,high_limit) n_seed safety_distance alpha_spectrum = [(nu_indOfRef, n_seed + 1.0 / (2.0 * pi**2.0) * (alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum)) 
                                                          | nu_indOfRef <- (map fst alpha_spectrum), nu_indOfRef >= low_limit, nu_indOfRef <= high_limit ]
    where
        -- integrate the quotient term for a specific index of refraction with the corresponding
        -- series of point of the quotient
        alpha_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
        --alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,high_limit)
        alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = low_integral + high_integral        
            where
                alpha_quotient_spectrum_low = alpha_quotients_low nu_indOfRef alpha_spectrum
                alpha_quotient_spectrum_high = alpha_quotients_high nu_indOfRef alpha_spectrum
                
                low_integral = evaluateIntegral int_method alpha_quotient_spectrum_low (low_low_limit,low_high_limit)
                high_integral = evaluateIntegral int_method alpha_quotient_spectrum_high (high_low_limit,high_high_limit)
                
                low_low_limit = last $ filter (< low_limit) $ map fst alpha_quotient_spectrum_low
                low_high_limit = fst $ last alpha_quotient_spectrum_low
                high_low_limit = fst $ head alpha_quotient_spectrum_high
                high_high_limit = head $ filter (> high_limit) $ map fst alpha_quotient_spectrum_high
                
                
                -- for a specific wavenumber generate the series of points, that need to 
                -- be integrated. This gives the term
                --
                -- alpha(nu') / (nu'² - nu²)
                --
                -- as a series the full spectral range of the alpha spectrum (nu') at the specific
                -- wavenumber nu
                alpha_quotients_low :: Double -> [(Double,Double)] -> [(Double,Double)]
                --alpha_quotients nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (snd a) <= 10.0]
                alpha_quotients_low nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (fst a) <= (nu_indOfRef - safety_distance), (snd a) <= 10.0] 
                    where
                        alpha_quotients_raw = map get_quotient alpha_spectrum
                        get_quotient (nu,alpha) = (nu,alpha / (nu**2.0 - nu_indOfRef**2.0))
                alpha_quotients_high :: Double -> [(Double,Double)] -> [(Double,Double)]
                alpha_quotients_high nu_indOfRef alpha_spectrum = [a | a <- alpha_quotients_raw, (fst a) >= (nu_indOfRef + safety_distance), (snd a) <= 10.0]
                    where
                        alpha_quotients_raw = map get_quotient alpha_spectrum
                        get_quotient (nu,alpha) = (nu,alpha / (nu**2.0 - nu_indOfRef**2.0))

-- calculate the real part of the index of refraction by Kramers Kronig relation from
-- the alpha spectra. Needs a seed value
indOfRef_n' :: (Double,Double) -> Double -> Double -> [(Double,Double)] -> [(Double,Double)]
indOfRef_n' (low_limit,high_limit) n_seed safety_distance alpha_spectrum = [(nu_indOfRef, n_seed + 1.0 / (2.0 * pi**2.0) * (alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum)) 
                                                          | nu_indOfRef <- (map fst alpha_spectrum), nu_indOfRef >= low_limit, nu_indOfRef <= high_limit ]
    where
        -- integrate the quotient term for a specific index of refraction with the corresponding
        -- series of point of the quotient
        alpha_integral :: (Double,Double) -> Double -> [(Double,Double)] -> Double
        --alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,high_limit)
        alpha_integral (low_limit,high_limit) nu_indOfRef alpha_spectrum = integral
            where
                alpha_quotient_spectrum = alpha_quotients nu_indOfRef alpha_spectrum
                
                -- low_integral = evaluateIntegral Akima alpha_quotient_spectrum (low_limit,nu_indOfRef - 0.0)
                -- high_integral = evaluateIntegral Akima alpha_quotient_spectrum (nu_indOfRef + 0.0,high_limit)
                
                integral = naiveIntegrate alpha_quotient_spectrum (low_limit, high_limit)
                
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
             [             (   e_i - e_m   ) ]
             [     3 * f * ( ------------- ) ]
             [             ( e_i + 2 * e_m ) ]
e_av = e_m * [ 1 + ------------------------- ]
             [             (   e_i - e_m   ) ]
             [     1 - f * ( ------------- ) ]
             [             ( e_i + 2 * e_m ) ]

e_i - e_m = a
e_i + 2 * e_m = b

e_av -> permittivity of mixture
e_m -> permittivity of matrix
e_i -> permittivity of inclusions
f -> volume fraction of inclusions             
-}
maxwellGarnet :: Double -> [(Double,Complex Double)] -> [(Double,Complex Double)] -> [(Double,Complex Double)]
maxwellGarnet volume_fraction epsilon1_spectrum epsilon2_spectrum = zip xVals epsilonAverVals
    where
        xVals = map fst epsilon1_spectrum
        e_iVals = map snd epsilon1_spectrum
        e_mVals = map snd epsilon2_spectrum 
        ll = length epsilon1_spectrum
        a = zipWith (-) e_iVals e_mVals
        b = zipWith (+) e_iVals e_m2Vals
        f = volume_fraction
        e_m2Vals = map (*2) e_mVals
        epsilonAverVals = [(e_mVals!!ind) * (1 + (((a !! ind) / (b !! ind)) `commul` f * 3) / (1 - ((a !! ind) / (b !! ind)) `commul` f)) | ind <- [0 .. (ll -1)]]


commul :: Complex Double -> Double -> Complex Double
a `commul` b = ((realPart a) * b) :+ ((imagPart a) * b)

naiveIntegrate :: [(Double,Double)] -> (Double,Double) -> Double
naiveIntegrate gridFunc (low_limit,high_limit) = sum parts
    where
        xVals = map fst gridFunc
        yVals = map snd gridFunc
        gridFuncInRange = [(xVals !! ind, yVals !! ind) | ind <- [0 .. (length gridFunc -1)], (xVals !! ind) >= low_limit && (xVals !! ind) <= high_limit]
        xValsIR = map fst gridFuncInRange
        yValsIR = map snd gridFuncInRange
        xDistances_middle = [0.5 * abs ((xValsIR !! ind) - (xValsIR !! (ind + 1)))  + 0.5 * abs ((xValsIR !! ind) - (xValsIR !! (ind + 1))) | ind <- [1 .. (length xValsIR - 2)]]
        xDistances_start = abs $ (xValsIR !! 0) - (xValsIR !! 1)
        xDistances_end = abs $ (last xValsIR) - (xValsIR !! (length xValsIR - 2))
        xDistances = [xDistances_start] ++ xDistances_middle ++ [xDistances_end]
        parts = zipWith (*) yValsIR xDistances
        
        
