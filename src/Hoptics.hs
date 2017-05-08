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
-}

import qualified Data.Text.IO as TextIO
import qualified Data.Text as Text
import Text.Read
import Data.Attoparsec.Text
import Data.Either.Unwrap (fromRight)
import Text.Printf
import System.Environment
import System.IO
--import Control.Applicative
import Control.Exception

-- modules specific for hoptics
import Parsers         -- parsers for text data
import Optics          -- calculation functions from spectra

x2meter = 1.0e-9      -- conversion factor from wavelength as in spectrum to wavelength in meter
x2invmeter = 1.0e2    -- conversion for x as wavenumber inverse meter

data XData = Wavenumber | Wavelength deriving (Show,Eq)

main = do
    putStrLn "         ********************"
    putStrLn "         *** HOptics v0.1 ***"
    putStrLn "         ********************\n"
    
    mainMenu

mainMenu = do
    -- set defaults
    -- spectrumMenu
    thickness <- return 100.0
    n_seed <- return 1.0
    spectralRange <- return (45000.0,590000.0)
    security_distance <- return 500.0
    spectrum_path_raw <- getArgs
    spectrum_path <- return $ head spectrum_path_raw
    unitOnX <- return Wavenumber
    
    
    putStrLn "\nHoptics"
    putStrLn "∟ Main Menu"
    putStrLn ""
    
    putStrLn "(1)  derive index of refraction from spectrum"
    putStrLn "(2)  mix two sets of indices of refraction"
    
    main_menu_input <- getLine
    case main_menu_input of
         "1" -> do
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
         "2" -> do
             putStrLn "not implemented yet"
         _   -> do
             putStrLn "enter correct number"
             mainMenu

spectrumMenu :: String -> Double -> (Double,Double) -> Double -> Double -> XData -> IO ()
spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX = do
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
    
    -- read in user selection for menu and apply changes or do the computation
    spectrumMenu_input <- getLine
    case spectrumMenu_input of
         "-1" -> do
             mainMenu
         "0" -> do
             spectrum_raw_content <- Control.Exception.try (TextIO.readFile spectrum_path) :: IO (Either SomeException Text.Text)
             case spectrum_raw_content of
                  Left exception -> do
                      putStrLn ("can't open file " ++ (show exception))
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
                  Right spectrum_raw_content -> do
                      let spectrum_xUnkwn = fromRight $ parseOnly parse_spectrum spectrum_raw_content
                      trans_spectrum_unordered <- if (unitOnX == Wavelength)
                                                     then do
                                                         let trans_spectrum_inMeter = map (scaleX x2meter) spectrum_xUnkwn
                                                         return $ wavelength2wavenumber trans_spectrum_inMeter
                                                     else do
                                                         let trans_spectrum_inInvMeter = map (scaleX x2invmeter) spectrum_xUnkwn
                                                         return $ trans_spectrum_inInvMeter
                      let trans_spectrum = order_spectrum trans_spectrum_unordered
                          alpha_spectrum = absorption_coeff (thickness * 1.0e-9) trans_spectrum
                          k_spectrum = indOfRef_k alpha_spectrum
                          n_spectrum = indOfRef_n spectralRange n_seed security_distance alpha_spectrum
                      
                      putStrLn "\ncalculating index of refraction... (can take some time)"  
                      
                      let spectrum_name = fromRight $ parseOnly parse_filename (Text.pack spectrum_path)
                          spectrum_basename = fst spectrum_name
                      
                      
                      trans_handle <- openFile ((spectrum_basename ++ "_trans.dat")) WriteMode
                      alpha_handle <- openFile ((spectrum_basename ++ "_alpha.dat")) WriteMode
                      k_handle <- openFile ((spectrum_basename ++ "_k.dat")) WriteMode
                      n_handle <- openFile ((spectrum_basename ++ "_n.dat")) WriteMode
                      
                      mapM_ (print_specpoint trans_handle) trans_spectrum
                      mapM_ (print_specpoint alpha_handle) alpha_spectrum
                      mapM_ (print_specpoint k_handle) k_spectrum
                      mapM_ (print_specpoint n_handle) n_spectrum

                      hClose trans_handle
                      hClose alpha_handle
                      hClose k_handle
                      hClose n_handle
                      
                      putStrLn "finished calculation"
                      mainMenu
         "1" -> do
             putStrLn ""
             putStrLn "enter file name of the spectrum"
             spectrum_path <- getLine
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
         "2" -> do
             putStrLn ""
             putStrLn "enter thickness of the sample [nm]"
             thickness_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) thickness_raw) of
                  Just x -> do
                      let thickness = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX 
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
         "3" -> do
             putStrLn ""
             putStrLn "spectral range for the computation (low,high)"
             spectralRange_raw <- getLine
             case ((readMaybe :: String -> Maybe (Double,Double)) spectralRange_raw) of
                  Just x -> do
                      let spectralRange = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX 
                  Nothing -> do
                      putStrLn "can not read this numbers, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX 
         "4" -> do
             putStrLn ""
             putStrLn "security distance around poles and boundaries in the dimension on x"
             security_distance_raw <- getLine
             case ((readMaybe :: String -> Maybe Double) security_distance_raw) of
                  Just x -> do
                      let security_distance = x
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
                  Nothing -> do
                      putStrLn "can not read this number, try again"
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX 
         "5" -> do
            putStrLn ""
            putStrLn "seed value fro the index of refraction"
            n_seed_raw <- getLine
            case ((readMaybe :: String -> Maybe Double) n_seed_raw) of
                 Just x -> do
                     let n_seed = x
                     spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
                 Nothing -> do
                     putStrLn "can not read this number, try again"
                     spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
         "6" -> do
             putStrLn ""
             putStrLn "wavenumber [cm⁻¹] (wn) or wavelength [nm] (wl) on x axis?"
             unitOnX_raw <- getLine
             case unitOnX_raw of
                  "wl" -> do
                      let unitOnX = Wavelength
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
                  "wn" -> do
                      let unitOnX = Wavenumber
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
                  _ -> do
                      let unitOnX = Wavenumber
                      spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX
         _ -> do
             spectrumMenu spectrum_path thickness spectralRange security_distance n_seed unitOnX


scaleX :: (Num a) => a -> (a,a) -> (a,a)
scaleX a (x,y) = (a*x,y)

scaleY :: (Num a) => a -> (a,a) -> (a,a)
scaleY a (x,y) = (x,a*y)

order_spectrum :: [(Double,Double)] -> [(Double,Double)]
order_spectrum spectrum
    | (fst $ head spectrum) < (fst $ last spectrum) = spectrum
    | (fst $ head spectrum) > (fst $ last spectrum) = reverse spectrum

print_specpoint :: Handle -> (Double,Double) -> IO()
print_specpoint handle (x,y) = do
    hPrintf handle "%+12.6e      " x
    hPrintf handle "%+12.6e      \n" y
