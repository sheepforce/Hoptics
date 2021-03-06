{- 
hoptics - command line tool for deriving optical constants from transmission
spectra and combining optical constants

freq   -> frequency in Hz
omega  -> frequency in Hz
nu     -> wavenumber in 1/m
lambda -> wavelength in m
alpha  -> exctinction coefficient in 1/m

n = n0 + ik
n  -> complex index of refraction
n0 -> real part
k  -> imaginary part

alpha = 1/d * ln(I0/I)
alpha -> absorption_coeff
d     -> sample thickness
I/I0  -> transmission

k = alpha / (4 pi nu)
nu   -> wavenumber

n0(nu) = n0_seed + 1/(2 pi^2) * \int \limits_0^\infty (alpha(nu')) / (nu'^2 - nu^2) d nu'
nu  -> wavenumber at which the index of refraction is calculated
nu' -> integration variable

Maxwell Garnet formalism in Hoptics
-}

import qualified Data.Text.IO as TextIO
import qualified Data.Text as Text
import Data.Complex
import Numeric.GSL
import Data.Either
import Text.Read
import Data.Attoparsec.Text
import Data.Either.Unwrap (fromRight,fromLeft)
import Text.Printf
import System.Environment
import System.IO
import Control.Exception
-- modules specific for hoptics
import Parsers         -- parsers for text data
import Optics          -- calculation functions from spectra

x2meter = 1.0e-9      -- conversion factor from wavelength as in spectrum to wavelength in meter
x2invmeter = 1.0e2    -- conversion for x as wavenumber inverse meter

-- data for unit on x-axis
data XData = Wavenumber | Wavelength deriving (Show,Eq)

-- choose a integration method for Kramers-Kronig relation.
-- Naive is implemented in Optics, the others are provided by the GNU Scientific library
data IntegrationMethod = Naive | Linear | Polynomial | CSpline | Akima deriving (Show,Eq)

-- main module, only showing the version
main = do
    putStrLn "         ********************"
    putStrLn "         *** HOptics v0.1 ***"
    putStrLn "         ********************"
    putStrLn ""
    
    -- call the IO routine for the main menu
    mainMenu

mainMenu = do
    -- set defaults for other menus
    -- spectrumMenu
    thickness <- return 100.0                           -- slab thickness in nano metre
    n_seed <- return 1.0                                -- seed value for real part of the index of refraction (Kramers Kronig)
    spectralRange <- return (45000.0,590000.0)          -- only applies to integration at Kramers-Kronig, which to which wavenumber
    security_distance <- return 500.0                   -- around the poles in Kramers-Kronig relation, how large is the gap where values are ommited for integration
    spectrum_path_raw <- getArgs                        -- get the arguments supplied with the call of hoptics executable
    spectrum_path <- do
        if ((length spectrum_path_raw) < 1)             -- and test if any were supplied
           then do                                      -- if no arguments were found
               return $ "spectrum.dat"                  -- set dummy file name
           else do
               return $ head spectrum_path_raw          -- if arguments were found, use them as guess for the input spectrum
    unitOnX <- return Wavenumber                        -- set wavenumbers as default unit on x-axis
    integration_method <- return Main.Akima             -- use robust akima interpolation/integration for Kramers-Kronig
    -- mixingMenu
    spectrum1_prefix <- do                              -- guess for prefix of spectrum 1
        if ((length spectrum_path_raw) < 1)             -- if no arguments are found
           then do                                      -- use a dummy file name
               return $ "spectrum.dat"
           else do                                      -- if arguments are found, use them as a guess
               return $ spectrum_path_raw !! 0
    spectrum2_prefix <- do                              -- same as above
        if ((length spectrum_path_raw) < 2)
           then do
               return $ "spectrum.dat"
           else do
               return $ spectrum_path_raw !! 1
    volume_fraction <- return 0.5                       -- assume a 1/1 mixture of inclusion and matrix
    magnetic_permittivity <- return 1.0                 -- assume a non magnetic material
    
    -- show infos of the current menu
    putStrLn "\nHoptics"
    putStrLn "∟ Main Menu"
    putStrLn ""
    
    putStrLn "(1)  derive index of refraction from spectrum"
    putStrLn "(2)  mix two sets of indices of refraction"
    
    -- get user selection for this menus
    main_menu_input <- getLine
    -- check if user input is valid and react accordingly
    case main_menu_input of
         "1" -> do
             -- call the menu for analysing measured spectra
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
         "2" -> do
             -- call the menu for mixing two already analysed spectra
             mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         _   -> do
             -- if no correct answer was given, go back to the main menu
             putStrLn "enter correct number"
             mainMenu

-- the menu and routines for analysing measured spectra
spectrumMenu :: String -> Double -> (Double,Double) -> Double -> Double -> XData -> IntegrationMethod -> IO ()
spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method = do
    -- show options for this menu
    putStrLn ""
    putStrLn "Hoptics"
    putStrLn "∟ Main Menu"
    putStrLn "  ∟ analyse spectrum"
    putStrLn ""
    putStrLn   "(-1) return to main menu"
    putStrLn $ "(0)  start computation"
    putStrLn $ "(1)  path to spectrum                                  " ++ (show spectrum_path)
    putStrLn $ "(2)  thickness of the sample [nm]                      " ++ (show thickness)
    putStrLn $ "(3)  spectral range for the calculation                " ++ (show spectralRange)
    putStrLn $ "(4)  security distance around poles and boundaries     " ++ (show security_distance)
    putStrLn $ "(5)  seed value for real part of index of refraction   " ++ (show n_seed)
    putStrLn $ "(6)  dimension on x axis                               " ++ (show unitOnX)
    putStrLn $ "(7)  integration method for Kramers Kronig             " ++ (show integration_method)
    
    -- read in user selection for menu and apply changes or do the computation
    spectrumMenu_input <- getLine
    case spectrumMenu_input of
         "-1" -> do
             mainMenu
         "0" -> do
             -- this will start with the computation
             -- try reading in the spectrum from the input file.
             -- if it fails, catch the exception, show the user what went wrong
             -- and let him try again.
             spectrum_raw_content <- Control.Exception.try (TextIO.readFile spectrum_path) :: IO (Either SomeException Text.Text)
             case spectrum_raw_content of
                  Left exception -> do
                      putStrLn ("can't open file " ++ (show exception))
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  Right spectrum_raw_content -> do
                      -- parse succesfully read in spectrum
                      let spectrum_xUnkwn = fromRight $ parseOnly parse_spectrum spectrum_raw_content
                      -- convert to wavenumber and order for increasing wavenumbers
                      trans_spectrum_unordered <- if (unitOnX == Wavelength)
                                                     then do
                                                         let trans_spectrum_inMeter = map (scaleX x2meter) spectrum_xUnkwn
                                                         return $ wavelength2wavenumber trans_spectrum_inMeter
                                                     else do
                                                         let trans_spectrum_inInvMeter = map (scaleX x2invmeter) spectrum_xUnkwn
                                                         return $ trans_spectrum_inInvMeter
                      let trans_spectrum = order_spectrum trans_spectrum_unordered
                          -- calculate the absorption coefficient for every wavelength
                          alpha_spectrum = absorption_coeff (thickness * 1.0e-9) trans_spectrum
                          -- the imaginary part of the index of refraction
                          k_spectrum = indOfRef_k alpha_spectrum
                          -- and the real part of the index of refraction
                          n0_spectrum
                              | integration_method == Naive = indOfRef_n' spectralRange n_seed security_distance alpha_spectrum
                              | integration_method == Main.Linear = indOfRef_n Numeric.GSL.Linear spectralRange n_seed security_distance alpha_spectrum
                              | integration_method == Main.Polynomial = indOfRef_n Numeric.GSL.Polynomial spectralRange n_seed security_distance alpha_spectrum
                              | integration_method == Main.CSpline = indOfRef_n Numeric.GSL.CSpline spectralRange n_seed security_distance alpha_spectrum
                              | integration_method == Main.Akima = indOfRef_n Numeric.GSL.Akima spectralRange n_seed security_distance alpha_spectrum
                              | otherwise = indOfRef_n Numeric.GSL.Akima spectralRange n_seed security_distance alpha_spectrum
                      
                      -- start outputting them, calculation happens here
                      putStrLn "\ncalculating index of refraction... (can take some time)"  
                      -- parse the file name and find its prefix and suffix. use the prefix as prefix for the output
                      let spectrum_name = fromRight $ parseOnly parse_filename (Text.pack spectrum_path)
                          spectrum_basename = fst spectrum_name
                      
                      -- open handles for the output spectra
                      trans_handle <- openFile (spectrum_basename ++ "_trans.dat") WriteMode
                      alpha_handle <- openFile (spectrum_basename ++ "_alpha.dat") WriteMode
                      k_handle <- openFile (spectrum_basename ++ "_k.dat") WriteMode
                      n0_handle <- openFile (spectrum_basename ++ "_n0.dat") WriteMode
                      -- monadically map a formatting output routine to the spectra
                      mapM_ (print_specpoint trans_handle) trans_spectrum
                      mapM_ (print_specpoint alpha_handle) alpha_spectrum
                      mapM_ (print_specpoint k_handle) k_spectrum
                      mapM_ (print_specpoint n0_handle) n0_spectrum
                      -- close the handles
                      hClose trans_handle
                      hClose alpha_handle
                      hClose k_handle
                      hClose n0_handle
                      
                      putStrLn "finished calculation"
                      mainMenu
         -- changes to the settings are processed from here
         "1" -> do
             putStrLn ""
             putStrLn "enter file name of the spectrum"
             spectrum_path <- getLine
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
         "2" -> do
             putStrLn ""
             putStrLn "enter thickness of the sample [nm]"
             thickness_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) thickness_raw) of
                  Just x -> do
                      let thickness = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method 
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
         "3" -> do
             putStrLn ""
             putStrLn "spectral range for the computation (low,high)"
             spectralRange_raw <- getLine
             case ((readMaybe :: String -> Maybe (Double,Double)) spectralRange_raw) of
                  Just x -> do
                      let spectralRange = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method 
                  Nothing -> do
                      putStrLn "can not read this numbers, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method 
         "4" -> do
             putStrLn ""
             putStrLn "security distance around poles and boundaries in the dimension on x"
             security_distance_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) security_distance_raw) of
                  Just x -> do
                      let security_distance = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method 
         "5" -> do
            putStrLn ""
            putStrLn "seed value fro the index of refraction"
            n_seed_raw <- getLine
            case ((readMaybe :: String -> Maybe Double) n_seed_raw) of
                 Just x -> do
                     let n_seed = x
                     spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                 Nothing -> do
                     putStrLn "can not read this number, try again"
                     spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
         "6" -> do
             putStrLn ""
             putStrLn "wavenumber [cm⁻¹] (wn) or wavelength [nm] (wl) on x axis?"
             unitOnX_raw <- getLine
             case unitOnX_raw of
                  "wl" -> do
                      let unitOnX = Wavelength
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  "wn" -> do
                      let unitOnX = Wavenumber
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  _ -> do
                      let unitOnX = Wavenumber
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
         "7" -> do
             putStrLn ""
             putStrLn "integration method for real part"
             putStrLn "Naive | Linear | Polynomial | CSpline | Akima"
             putStrLn "Naive recommended for UVVis, Akima otherwise"
             integration_method_raw <- getLine
             case integration_method_raw of
                  "Naive" -> do
                      let integration_method = Main.Naive
                      let security_distance = 0.0
                      putStrLn "security distance does not apply to this method"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  "Linear" -> do
                      let integration_method = Main.Linear
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  "Polynomial" -> do
                      let integration_method = Main.Polynomial
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  "CSpline" -> do
                      let integration_method = Main.CSpline
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  "Akima" -> do
                      let integration_method = Main.Akima
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                  _ -> do
                      putStrLn "not a valid choice, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method
                      
                  
         _ -> do
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX integration_method

-- menu and routines for mixing dielectric permittivities by Maxwell-Garnet
mixingMenu :: String -> String -> Double -> Double -> IO ()
mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity = do
    -- show options for this menu
    putStrLn ""
    putStrLn "Hoptics"
    putStrLn "∟ Main Menu"
    putStrLn "  ∟ mix spectra"
    putStrLn ""
    putStrLn   "(-1) return to main menu"
    putStrLn $ "(0)  start computation"
    putStrLn $ "(1)  prefix of spectrum of component 1 (inclusion)       " ++ (show spectrum1_prefix)
    putStrLn $ "(2)  prefix of spectrum of component 2 (matrix)          " ++ (show spectrum2_prefix)
    putStrLn $ "(3)  volume fraction of component 1                      " ++ (show volume_fraction)
    putStrLn $ "(4)  magnetic permittivity of the mixture                " ++ (show magnetic_permittivity)
    -- if a strange input is found, echo a warning
    if (volume_fraction > 0.5)
       then do
           putStrLn "WARNING: volume fraction of the inclusion is higher than 0.5. BE SURE YOU WANT TO DO THIS"
       else do
           return ()
    
    -- get choice of the user
    mixingMenu_input <- getLine
    -- check for its validity and act accordingly
    case mixingMenu_input of
         "-1" -> do
             mainMenu
         "0" -> do
             -- try reading 2 spectra, each consisting of the real part (n0) and the complex part (k), catching possible exceptions
             spectrum1_n0_raw <- Control.Exception.try (TextIO.readFile (spectrum1_prefix ++ "_n0.dat")) :: IO (Either SomeException Text.Text)
             spectrum2_n0_raw <- Control.Exception.try (TextIO.readFile (spectrum2_prefix ++ "_n0.dat")) :: IO (Either SomeException Text.Text)
             spectrum1_k_raw <- Control.Exception.try (TextIO.readFile (spectrum1_prefix ++ "_k.dat")) :: IO (Either SomeException Text.Text)
             spectrum2_k_raw <- Control.Exception.try (TextIO.readFile (spectrum2_prefix ++ "_k.dat")) :: IO (Either SomeException Text.Text)
             if (isLeft spectrum1_n0_raw || isLeft spectrum2_n0_raw || isLeft spectrum1_k_raw || isLeft spectrum2_k_raw) 
                -- exception handling
                then do
                    putStrLn ("can't open some files, look at the errors and try again")
                    if (isLeft spectrum1_n0_raw)
                       then do
                           print (fromLeft spectrum1_n0_raw)
                       else do
                           return ()
                    if (isLeft spectrum2_n0_raw)
                       then do
                           print (fromLeft spectrum2_n0_raw)
                       else do
                           return ()
                    if (isLeft spectrum1_k_raw)
                       then do
                           print (fromLeft spectrum1_k_raw)
                       else do
                           return ()
                    if (isLeft spectrum2_k_raw)
                       then do
                           print (fromLeft spectrum2_k_raw)
                       else do
                           return ()
                    mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
                -- doing the actual calculation if everything else went fine
                else do
                    -- parse the spectra
                    let spectrum1_n0 = fromRight $ parseOnly parse_spectrum (fromRight spectrum1_n0_raw)
                        spectrum2_n0 = fromRight $ parseOnly parse_spectrum (fromRight spectrum2_n0_raw)
                        spectrum1_k_= fromRight $ parseOnly parse_spectrum (fromRight spectrum1_k_raw)
                        spectrum2_k = fromRight $ parseOnly parse_spectrum (fromRight spectrum2_k_raw)
                        -- combine the two pars (n0 and k) read for each from file to the complex index of refraction
                        spectrum1_indOfRef = nkSep2complexIndOfRef spectrum1_n0 spectrum1_k_
                        spectrum2_indOfRef = nkSep2complexIndOfRef spectrum2_n0 spectrum2_k
                        -- convert index of refraction to dielectric permittivity
                        spectrum1_permittivity = n2epsilon magnetic_permittivity spectrum1_indOfRef
                        spectrum2_permittivity = n2epsilon magnetic_permittivity spectrum2_indOfRef
                        -- apply Maxwell Garnet formalism to the two sets of permittivities
                        spectrum_permittivity_comb_MaxwellGarnet = maxwellGarnet volume_fraction spectrum1_permittivity spectrum2_permittivity 
                        -- convert back to index of refraction
                        spectrum_n_comb_MaxwellGarnet = epsilon2n magnetic_permittivity spectrum_permittivity_comb_MaxwellGarnet
                        -- get the real part from the index of refraction
                        spectrum_n0_comb_MaxwellGarnet = zip (map fst spectrum_n_comb_MaxwellGarnet) (map realPart $ map snd spectrum_n_comb_MaxwellGarnet)
                        -- get the complex part from the index of refraction
                        spectrum_k_comb_MaxwellGarnet = zip (map fst spectrum_n_comb_MaxwellGarnet) (map imagPart $ map snd spectrum_n_comb_MaxwellGarnet)
                    
                    -- start actual calculation
                    putStrLn ""
                    putStrLn "mixing spectra now by Maxwell Garnet formalism"
                    -- open a output handle
                    maxwellGarnetHandle_n0 <- openFile (spectrum1_prefix ++ "+" ++ spectrum2_prefix ++ "_n0_MaxwellGarnet.dat") WriteMode
                    maxwellGarnetHandle_k <- openFile (spectrum1_prefix ++  "+" ++ spectrum2_prefix ++ "_k_MaxwellGarnet.dat") WriteMode
                    -- apply monadic print
                    mapM_ (print_specpoint maxwellGarnetHandle_n0) spectrum_n0_comb_MaxwellGarnet
                    mapM_ (print_specpoint maxwellGarnetHandle_k) spectrum_k_comb_MaxwellGarnet
                    -- close the handles again
                    hClose maxwellGarnetHandle_k
                    hClose maxwellGarnetHandle_n0 
                    
                    putStr "finished spectra mixing"
                    
                    mainMenu
         -- allow changes to the settings
         "1" -> do
             putStrLn ""
             putStrLn "enter prefix of spectrum 1 (\"spectrum\" for \"spectrum_n.dat\" and \"spectrum_k.dat\")"
             spectrum1_prefix <- getLine
             mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         "2" -> do
             putStrLn ""
             putStrLn "enter prefix of spectrum 2 (\"spectrum\" for \"spectrum_n.dat\" and \"spectrum_k.dat\")"
             spectrum2_prefix <- getLine
             mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         "3" -> do
             putStrLn ""
             putStrLn "volume fraction of component 1"
             volume_fraction_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) volume_fraction_raw) of
                  Just x -> do
                      let volume_fraction = x
                      mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         "4" -> do
             putStrLn ""
             putStrLn "relative magnetic permittivity"
             magnetic_permittivity_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) magnetic_permittivity_raw) of
                  Just x -> do
                      let magnetic_permittivity = x
                      mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         _ -> do
             mixingMenu spectrum1_prefix spectrum2_prefix volume_fraction magnetic_permittivity
         
-- scaling only the x-value of a point
scaleX :: (Num a) => a -> (a,a) -> (a,a)
scaleX a (x,y) = (a*x,y)

-- scaling only the y-value of a point
scaleY :: (Num a) => a -> (a,a) -> (a,a)
scaleY a (x,y) = (x,a*y)

-- make sure, that the spectrum is in increasing order, if it was in decreasing order before
order_spectrum :: [(Double,Double)] -> [(Double,Double)]
order_spectrum spectrum
    | (fst $ head spectrum) < (fst $ last spectrum) = spectrum
    | (fst $ head spectrum) > (fst $ last spectrum) = reverse spectrum

-- formatted printing of spectrum point
print_specpoint :: Handle -> (Double,Double) -> IO()
print_specpoint handle (x,y) = do
    hPrintf handle "%+12.6e      " x
    hPrintf handle "%+12.6e      \n" y
