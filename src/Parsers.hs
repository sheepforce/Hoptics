{-
implements parsers for hoptics
-}

module Parsers
( parse_spectrum
, parse_filename
) where
import qualified Data.Text as Text
import Data.Attoparsec.Text


-- takes simple two column ascii, separated with spaces or tabs
parse_spectrum :: Parser [(Double,Double)]
parse_spectrum = do
    many' (many' space <* endOfLine)
    spectrum <- many1 $ parse_spectrum_line
    endOfInput
    return $ spectrum

-- parse a single line of a two column ascii spectrum
parse_spectrum_line :: Parser (Double,Double)
parse_spectrum_line = do
    skipSpace
    x <- double
    skipSpace
    y <- double
    skipSpace
    return $ (x,y)

parse_filename :: Parser (String,String)
parse_filename = do
    name <- manyTill' anyChar (char '.')
    suffix <- many1 anyChar 
    return $ (name,suffix)
