
module Main where

import BioInf.Alignment.DnaProtein



main = do
  dna <- getLine
  pro <- getLine
  let (s,bs) = dnaProtein dna pro
  print s
  mapM_ (\(o,u) -> putStrLn o >> putStrLn u) $ take 1 bs
