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
  { dna           = def &= help "DNA fasta file to read"
--  , dnawindow     = 2000
  , protein       = def &= help "Protein fasta file to read"
--  , proteinwindow = 10000000
  , windowMult    =   3 &= help "window of k nucleotides for each amino acid (actually 2k, as we use sliding windows)"
  , blastMatrix   = def &= help "Blast matrix (PAM / BLOSUM) to use"
  , insertAA      = -50 &= help "cost for inserting an amino acid (indel)"
  , deleteAA      = -50 &= help "cost for deleting an amino acid (indel)"
  , rf1S          = -50 &= help "cost for aligning only two nucleotides with an AA and frame shifting by 1"
  , rf1delS       = -50 &= help "cost for deleting two nucleotides and frame shifting by 1"
  , rf2S          = -50 &= help "cost for aligning only one nucleotide with an AA and frame shifting by 2"
  , rf2delS       = -50 &= help "cost for deleting a nucleotide and frame shifting by 1"
  , minScore      = -999999 &= help "display only scores above this threshold"
  , parallelism   = 16 &= help "maximum parallelism (should be set 2-4x or more of the number of CPUs"
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
        let ll = if null bs then 0 else length . takeWhile (=='.') . toList . snd $ head bs
        printf "DNA: %s @ %d   |||   Protein: %s @ %d\n"
                (B.unpack $ _identifier d)
                (offD + fromIntegral ll)
                (B.unpack $ _identifier p)
                (unOff $ _offset p)
        printf "DNA length: %d Protein length: %d\n" (B.length inpD) (B.length $ _fasta p)
        printf "Score: %d\n" s
        if null bs then putStrLn "NO ALIGNMENT?" else do
          let tt = length . takeWhile (/='.') . drop ll . toList . snd $ head bs
              os = chunksOf 100 . take tt . drop ll . toList . fst $ head bs
              us = chunksOf 100 . take tt . drop ll . toList . snd $ head bs
          zipWithM_ (\o u -> putStrLn o >> putStrLn u) os us

