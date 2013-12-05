{-# LANGUAGE BangPatterns #-}
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
import Control.Parallel.Strategies
import Data.DList (toList)

import Biobase.Fasta
import Biobase.Fasta.Import
import Biobase.SubstMatrix
import Biobase.SubstMatrix.Import
import Biobase.Primary

import BioInf.Alignment.DnaProtein



data Option = Option
  { dna           :: String
--  , dnawindow     :: Int
  , protein       :: String
--  , proteinwindow :: Int
  , windowMult    :: Int
  , blastMatrix   :: String
  , insertAA      :: Int
  , deleteAA      :: Int
  , rf1S          :: Int
  , rf1delS       :: Int
  , rf2S          :: Int
  , rf2delS       :: Int
  , minScore      :: Int
  , parallelism   :: Int
  }
  deriving (Show,Data,Typeable)

option = Option
  { dna           = def
--  , dnawindow     = 2000
  , protein       = def
--  , proteinwindow = 10000000
  , windowMult    =   3
  , blastMatrix   = def
  , insertAA      = -50
  , deleteAA      = -50
  , rf1S          = -50
  , rf1delS       = -50
  , rf2S          = -50
  , rf2delS       = -50
  , minScore      = -999999
  , parallelism   = 16
  }

main = do
  Option{..} <- cmdArgs option
  ps <- if null protein then return [] else runResourceT $ sourceFile protein $= parseFastaWindows 999999999 $$ consume
  forM_ ps $ \p -> do
    let lenP = fromIntegral $ B.length $ _fasta p
    ds <- if null dna     then return [] else runResourceT $ sourceFile dna     $= parseFastaWindows (lenP * windowMult)  $$ consume
    when (null blastMatrix) $ error "require Blast matrix"
    mat <- fromFile blastMatrix
    let !n3m = mkNuc3SubstMatrix mat
    let !n2m = mkNuc2SubstMatrix max id mat
    let !n1m = mkNuc1SubstMatrix max id mat
    let xs = [ (d,inpD,offD,p, dnaProtein n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS inpD (B.unpack $ _fasta p))
             | d <- ds
             , let inpD = _past d `B.append` _fasta d
             , let offD = (unOff $ _offset d) - (fromIntegral . B.length $ _past d)
             ]
    forM_ (xs `using` parBuffer parallelism (evalTuple5 r0 r0 r0 r0 (evalTuple2 rdeepseq r0))) $ \(d,inpD,offD,p,(s,bs)) -> do
      when (s>=minScore) $ do
        printf "DNA: %s @ %d   |||   Protein: %s @ %d\n"
                (B.unpack $ _identifier d)
                offD
                (B.unpack $ _identifier p)
                (unOff $ _offset p)
        printf "DNA length: %d Protein length: %d\n" (B.length inpD) (B.length $ _fasta p)
        printf "Score: %d\n" s
        if null bs then putStrLn "NO ALIGNMENT?" else do
          let os = chunksOf 100 . toList . fst $ head bs
              us = chunksOf 100 . toList . snd $ head bs
          zipWithM_ (\o u -> putStrLn o >> putStrLn u) os us

