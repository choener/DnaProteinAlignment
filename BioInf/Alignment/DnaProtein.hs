{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeOperators #-}

{-# OPTIONS_GHC -fno-liberate-case #-}

module BioInf.Alignment.DnaProtein (dnaProtein) where

import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream.Monadic as M
import qualified Data.Vector.Fusion.Stream.Monadic
import Data.Vector.Fusion.Util
import System.IO.Unsafe
import qualified Data.Vector.Unboxed as V
import Control.Monad
import qualified Data.Array.IArray as A

import Biobase.Primary
import Data.PrimitiveArray as PA hiding (map)
import Data.PrimitiveArray.Zero as PA hiding (map)
import Data.Array.Repa.Index.Points
import ADP.Fusion.Multi
import ADP.Fusion.Multi.Classes
import ADP.Fusion.Multi.Empty
import ADP.Fusion
import ADP.Fusion.Chr
import ADP.Fusion.Table
import ADP.Fusion.Empty
import ADP.Fusion.None
import Biobase.SubstMatrix

import BioInf.Alignment.DnaProtein.Score
import BioInf.Alignment.DnaProtein.Pretty

-- |

dnaProtein n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS dna' protein' = ( _RP ! (Z:.pointL 0 nD:.pointL 0 nP), bt ) where
  ts@(_F0P, _F1P, _F2P, _LP, _RP) = unsafePerformIO $ fillTables n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS dna protein
  nD = V.length dna
  nP = V.length protein
  dna = mkPrimary dna'
  protein = V.fromList protein'
  bt = backtrack n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS dna protein ts
{-# NOINLINE dnaProtein #-}

