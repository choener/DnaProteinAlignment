
module Main where

import BioInf.Alignment.DnaProtein



main = do
  dna <- getLine
  pro <- getLine
  print $ dnaProtein dna pro
