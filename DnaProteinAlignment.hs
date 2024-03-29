{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE ScopedTypeVariables #-}
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
import qualified Text.PrettyPrint.ANSI.Leijen as PP

import Biobase.Fasta
import Biobase.Fasta.Import
import Biobase.SubstMatrix
import Biobase.SubstMatrix.Import
import Biobase.Primary

import BioInf.Alignment.DnaProtein
import BioInf.Alignment.DnaProtein.Pretty



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
  , minadjScore   :: Double
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
  , insertAA      = -10 &= help "cost for inserting an amino acid (indel)"
  , deleteAA      = -15 &= help "cost for deleting an amino acid (indel)"
  , rf1S          = -20 &= help "cost for aligning only two nucleotides with an AA and frame shifting by 1"
  , rf1delS       = -45 &= help "cost for deleting two nucleotides and frame shifting by 1"
  , rf2S          = -60 &= help "cost for aligning only one nucleotide with an AA and frame shifting by 2"
  , rf2delS       = -75 &= help "cost for deleting a nucleotide and frame shifting by 1"
  , minScore      = -999999 &= help "display only scores above this threshold"
  , minadjScore   = -999999 &= help "minimal (score / protein length)"
  , parallelism   = 16 &= help "maximum parallelism (should be set 2-4x or more of the number of CPUs"
  }

data Direction
  = Forward
  | Backward
  deriving (Eq,Ord,Show)

main = do
  o@Option{..} <- cmdArgs option
  ps <- if null protein then return [] else runResourceT $ sourceFile protein $= parseFastaWindows 999999999 $$ consume
  forM_ ps $ \p -> do
    let lenP = fromIntegral $ B.length $ _fasta p
    ds <- if null dna     then return [] else runResourceT $ sourceFile dna     $= parseFastaWindows (lenP * windowMult)  $$ consume
    when (null blastMatrix) $ error "require Blast matrix"
    mat <- fromFile blastMatrix
    let !n3m = mkNuc3SubstMatrix mat
    let !n2m = mkNuc2SubstMatrix max id mat
    let !n1m = mkNuc1SubstMatrix max id mat
    let xsF = [ (d,inpD,offD,p,Forward, dnaProtein n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS inpD (B.unpack $ _fasta p))
              | d <- ds
              , let inpD = _past d `B.append` _fasta d
              , let offD = (unOff $ _offset d) - (fromIntegral . B.length $ _past d)
              ]
    let xsB = [ (d,inpD,offD,p,Backward,  dnaProtein n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS inpD (B.unpack $ _fasta p))
              | d <- ds
              , let inpD = B.map compl . B.reverse $ _past d `B.append` _fasta d
              , let offD = (unOff $ _offset d) - (fromIntegral . B.length $ _past d)
              ]
    forM_ ((xsF++xsB) `using` parBuffer parallelism (evalTuple6 r0 r0 r0 r0 r0 (evalTuple2 rdeepseq r0)))
           $ \(d,inpD,offD,p,dir,(s,bs)) -> do
      let sa :: Double = fromIntegral s / fromIntegral (B.length $ _fasta p) :: Double
      when (s>=minScore && sa>=minadjScore) $ do
        let llF = if null bs then 0 else length . takeWhile isLOC           . toList . head $ bs
        let llR = if null bs then 0 else length . takeWhile isLOC . reverse . toList . head $ bs
        let frs1s = filter (\case (FRS [_,_] [_] _) -> True ; z->False) . toList . head $ bs
        let frs2s = filter (\case (FRS [_]   [_] _) -> True ; z->False) . toList . head $ bs
        printf "DNA: %s %s @ %s %d\nProtein: %s %s @ %d\n"
                (B.unpack $ _identifier d)
                (B.unpack $ _description d)
                (show dir)
                (offD + if dir==Forward then (fromIntegral llF) else (fromIntegral llR))
                (B.unpack $ _identifier p)
                (B.unpack $ _description p)
                (unOff $ _offset p)
        printf "DNA length: %d Protein length: %d\n" (B.length inpD) (B.length $ _fasta p)
        printf "1 Nt shifts: %d   |||   2 Nt shifts: %d\n"
                (length frs1s)
                (length frs2s)
        printf "Score: %d   Length-adjusted: %.2f\n\n" s sa
        if null bs then putStrLn "NO ALIGNMENT?" else do
          let tt = length . takeWhile isPPS . drop llF . toList . head $ bs
              cs = chunksOf 30 . take tt . drop llF . toList . head $ bs
              llIdx = if dir==Forward then llF else llR
          foldM_ (\pos pps -> let adv = advancePos pos pps
                              in  PP.putDoc (pps2doc pos adv pps) >> return (advancePos pos pps)
                 ) (fromIntegral offD + llIdx + 1, 1) cs
        putStrLn ""

compl :: Char -> Char
compl 'A' = 'T'
compl 'T' = 'A'
compl 'C' = 'G'
compl 'G' = 'C'
compl 'a' = 't'
compl 't' = 'a'
compl 'c' = 'g'
compl 'g' = 'c'
compl  x  =  x

advancePos :: (Int,Int) -> [PPS] -> (Int,Int)
advancePos = foldl go where
  go (l,r) (PPS cs as _) = (l + length cs, r + length as)
  go (l,r) (FRS cs as _) = (l + length cs, r + length as)
  go (l,r) (LOC _     _) = (l + 1        , r            )

pps2doc :: (Int,Int) -> (Int,Int) -> [PPS] -> PP.Doc
pps2doc (pl,pr) (apl,apr) xs
  =      PP.text (printf "%8d " pl) PP.<> us PP.<> (PP.text (printf " %8d" (apl-1)))
  PP.<$> PP.text (printf "%8d " pr) PP.<> os PP.<> (PP.text (printf " %8d" (apr-1)))
  PP.<$> PP.empty PP.<$> PP.empty
  where
  us = PP.hcat $ map upper xs
  os = PP.hcat $ map lower xs
  upper (PPS cs as  k) =           colorize k . PP.text . take 3 $ map fromNuc cs ++ repeat '-'
  upper (FRS cs as  k) = PP.underline . PP.bold . colorize k . PP.text . take 3 $ map fromNuc cs ++ repeat '-'
  upper (LOC c      k) =           colorize k . PP.text          $ [fromNuc c]
  lower (PPS cs []  k) =           colorize k . PP.text          $ "-  "
  lower (PPS cs [a] k) =           colorize k . PP.text . take 3 $ [a]            ++ repeat ' '
  lower (FRS cs [a] k) = PP.bold . colorize k . PP.text . take 3 $ [a]            ++ repeat ' '
  lower (LOC c      k) =           colorize k . PP.text          $ "."
  colorize k
    | k>5       = PP.cyan
    | k>0       = PP.blue
    | k<(-5)    = PP.red
    | k<0       = PP.yellow
    | otherwise = id

