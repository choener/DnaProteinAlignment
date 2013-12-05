{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

module Main where

import Data.Conduit.Lazy
import qualified Data.ByteString.Char8 as B
import Control.Monad
import Data.Conduit
import Data.Conduit.List (consume)
import Data.Conduit.Binary (sourceFile)
import System.Console.CmdArgs
import Text.Printf
import Control.Applicative
import Bio.Core.Sequence (Offset (..))
import Data.List.Split (chunksOf)

import Biobase.Fasta
import Biobase.Fasta.Import
import Biobase.SubstMatrix
import Biobase.SubstMatrix.Import
import Biobase.Primary

import BioInf.Alignment.DnaProtein



data Option = Option
  { dna           :: String
  , dnawindow     :: Int
  , protein       :: String
  , proteinwindow :: Int
  , blastMatrix   :: String
  }
  deriving (Show,Data,Typeable)

option = Option
  { dna           = def
  , dnawindow     = 2000
  , protein       = def
  , proteinwindow = 500
  , blastMatrix   = def
  }

main = do
  Option{..} <- cmdArgs option
  ds <- if null dna     then return [] else runResourceT $ sourceFile dna     $= parseFastaWindows dnawindow     $$ consume
  ps <- if null protein then return [] else runResourceT $ sourceFile protein $= parseFastaWindows proteinwindow $$ consume
  when (null blastMatrix) $ error "require Blast matrix"
  mat <- fromFile blastMatrix
  let n3m = mkNuc3SubstMatrix mat
  let n2m = mkNuc2SubstMatrix max id mat
  let n1m = mkNuc2SubstMatrix max id mat
  forM_ ds $ \d -> forM_ ps $ \p -> do
    printf "DNA: %s @%d   |||   Protein: %s @%d\n"
            (B.unpack $ _identifier d)
            (unOff $ _offset d)
            (B.unpack $ _identifier p)
            (unOff $ _offset p)
    printf "DNA length: %d Protein length: %d\n" (B.length $ _fasta d) (B.length $ _fasta p)
    let (s,bs) = dnaProtein (mkPrimary $ _fasta d) (B.unpack $ _fasta p)
    printf "Score: %d\n" s
    if null bs then putStrLn "NO ALIGNMENT?" else do
      let os = chunksOf 100 . fst $ head bs
          us = chunksOf 100 . snd $ head bs
      zipWithM_ (\o u -> putStrLn o >> putStrLn u) os us
  {-
  dna <- getLine
  pro <- getLine
  let (s,bs) = dnaProtein dna pro
  print s
  mapM_ (\(o,u) -> putStrLn o >> putStrLn u) $ take 1 bs
  -}
