{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}

{-# OPTIONS_GHC -fno-liberate-case #-}

module BioInf.Alignment.DnaProtein.Common where

import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream.Monadic as M
import qualified Data.Vector.Fusion.Stream.Monadic
import Data.Vector.Fusion.Util
import System.IO.Unsafe
import qualified Data.Vector.Unboxed as V
import Control.Monad
import qualified Data.Array.IArray as A
import Data.DList (DList)

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



data SigDnaPro {-Monad-} m {-NT-} nt hResT {-T-} a c e = SigDnaPro
  {delamino :: nt -> (Z:.():.a) -> nt
  ,lcldel :: nt -> nt
  ,nilnil :: (Z:.e:.e) -> nt
  ,rf1amino :: nt -> (Z:.c:.a) -> (Z:.c:.()) -> nt
  ,rf1del :: nt -> (Z:.c:.()) -> (Z:.c:.()) -> nt
  ,rf2amino :: nt -> (Z:.c:.a) -> nt
  ,rf2del :: nt -> (Z:.c:.()) -> nt
  ,stayamino :: nt -> (Z:.c:.a) -> (Z:.c:.()) -> (Z:.c:.()) -> nt
  ,staydel :: nt -> (Z:.c:.()) -> (Z:.c:.()) -> (Z:.c:.()) -> nt
  ,eatdel :: nt -> (Z:.c:.()) -> nt
  ,h :: Data.Vector.Fusion.Stream.Monadic.Stream m nt -> m hResT}

grammarDnaPro SigDnaPro{..} {-NT-} _F0P _F1P _F2P _LP _RP {-T-} a c e =
  ((_F0P
   ,delamino <<< _F0P % (T:!None:!a)
     ||| lcldel <<< _LP
     ||| nilnil <<< (T:!e:!e)
     ||| rf1amino <<< _F1P % (T:!c:!a) % (T:!c:!None)
     ||| rf1del <<< _F1P % (T:!c:!None) % (T:!c:!None)
     ||| rf2amino <<< _F2P % (T:!c:!a)
     ||| rf2del <<< _F2P % (T:!c:!None)
     ||| stayamino <<< _F0P % (T:!c:!a) % (T:!c:!None) % (T:!c:!None)
     ||| staydel <<< _F0P % (T:!c:!None) % (T:!c:!None) % (T:!c:!None) ... h)
  ,(_F1P
   ,delamino <<< _F1P % (T:!None:!a)
     ||| lcldel <<< _LP
     ||| nilnil <<< (T:!e:!e)
     ||| rf1amino <<< _F2P % (T:!c:!a) % (T:!c:!None)
     ||| rf1del <<< _F2P % (T:!c:!None) % (T:!c:!None)
     ||| rf2amino <<< _F0P % (T:!c:!a)
     ||| rf2del <<< _F0P % (T:!c:!None)
     ||| stayamino <<< _F1P % (T:!c:!a) % (T:!c:!None) % (T:!c:!None)
     ||| staydel <<< _F1P % (T:!c:!None) % (T:!c:!None) % (T:!c:!None) ... h)
  ,(_F2P
   ,delamino <<< _F2P % (T:!None:!a)
     ||| lcldel <<< _LP
     ||| nilnil <<< (T:!e:!e)
     ||| rf1amino <<< _F0P % (T:!c:!a) % (T:!c:!None)
     ||| rf1del <<< _F0P % (T:!c:!None) % (T:!c:!None)
     ||| rf2amino <<< _F1P % (T:!c:!a)
     ||| rf2del <<< _F1P % (T:!c:!None)
     ||| stayamino <<< _F2P % (T:!c:!a) % (T:!c:!None) % (T:!c:!None)
     ||| staydel <<< _F2P % (T:!c:!None) % (T:!c:!None) % (T:!c:!None) ... h)
  ,(_LP,eatdel <<< _LP % (T:!c:!None) ||| nilnil <<< (T:!e:!e) ... h)
  ,(_RP,eatdel <<< _RP % (T:!c:!None) ||| lcldel <<< _F0P ||| lcldel <<< _F1P ||| lcldel <<< _F2P ||| nilnil <<< (T:!e:!e) ... h))
{-# INLINE grammarDnaPro #-}

type PPT = PA.Unboxed (Z:.PointL:.PointL) Int
type BtPPT = DefBtTbl Id (Z:.PointL:.PointL) Int (DList Char,DList Char)
type FunT = (Z:.PointL:.PointL) -> Id (M.Stream Id (DList Char,DList Char))
type Tbl = PA.Unboxed (Z:.PointL:.PointL) Int


