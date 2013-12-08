{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}

{-# OPTIONS_GHC -fno-liberate-case #-}

module BioInf.Alignment.DnaProtein.Score (fillTables,algScore) where

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

import BioInf.Alignment.DnaProtein.Common





-- |

algScore
  :: Monad m
  => Nuc3SubstMatrix
  -> Nuc2SubstMatrix
  -> Nuc1SubstMatrix
  -> Int -> Int -> Int -> Int -> Int -> Int
  -> SigDnaPro m Int Int Char Nuc ()
algScore n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS = SigDnaPro
  { lcldel    = \z                                     -> z
  , nilnil    = \z                                     -> 0
  , delamino  = \z (Z:.():.a)                          -> z + insertAA
  , rf1amino  = \z (Z:.c1:.a)  (Z:.c2:.())             -> z + n2m A.! (c1,c2,a) + rf1S
  , rf1del    = \z (Z:.c1:.()) (Z:.c2:.())             -> z + rf1delS
  , rf2amino  = \z (Z:.c:.a)                           -> z + n1m A.! (c,a) + rf2S
  , rf2del    = \z (Z:.c:.())                          -> z + rf2delS
  , stayamino = \z (Z:.c1:.a)  (Z:.c2:.()) (Z:.c3:.()) -> z + n3m A.! (c1,c2,c3,a)
  , staydel   = \z (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> z + deleteAA
  , eatdel    = \z (Z:.c:.())                          -> z
  , h         = M.foldl' max (-999999)
  }
{-# INLINE algScore #-}



-- |

fillTables
  :: Nuc3SubstMatrix
  -> Nuc2SubstMatrix
  -> Nuc1SubstMatrix
  -> Int -> Int -> Int -> Int -> Int -> Int
  -> V.Vector Nuc -> V.Vector Char
  -> IO ( Tbl, Tbl, Tbl, Tbl, Tbl )
fillTables !n3m !n2m !n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS dna protein = do
  let nD = V.length dna
  let nP = V.length protein
  f0p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) (-999999)
  f1p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) (-999999)
  f2p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) (-999999)
  lp'  <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) (-999999)
  rp'  <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) (-999999)
  fillFive $ grammarDnaPro
               (algScore n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS)
               (mTbl (Z:.EmptyT:.EmptyT) f0p')
               (mTbl (Z:.EmptyT:.EmptyT) f1p')
               (mTbl (Z:.EmptyT:.EmptyT) f2p')
               (mTbl (Z:.EmptyT:.EmptyT) lp')
               (mTbl (Z:.EmptyT:.EmptyT) rp')
               (chr protein) (chr dna)
               Empty
  f0p <- freeze f0p'
  f1p <- freeze f1p'
  f2p <- freeze f2p'
  lp  <- freeze lp'
  rp  <- freeze rp'
  return (f0p,f1p,f2p,lp,rp)
{-# INLINE fillTables #-}

-- type TF = ( MutArr IO Tbl, (Z:.PointL:.PointL) -> IO Int )

-- fillFive :: TF -> IO () -- ( TF, TF, TF, TF, TF) -> IO ()
fillFive ( (MTbl _ tf0p, f_f0p)
         , (MTbl _ tf1p, f_f1p)
         , (MTbl _ tf2p, f_f2p)
         , (MTbl _ tlp , f_lp )
         , (MTbl _ trp , f_rp )
         ) = do
  let (_,Z:.PointL(0:.nD):.PointL(0:.nP)) = boundsM tf0p
  forM_ [0 .. nD] $ \k -> forM_ [0 .. nP] $ \l -> do
    let i = (Z:.pointL 0 k:.pointL 0 l)
    f_f0p i >>= writeM tf0p i
    f_f1p i >>= writeM tf1p i
    f_f2p i >>= writeM tf2p i
    f_lp  i >>= writeM tlp  i
    f_rp  i >>= writeM trp  i
{-# INLINE fillFive #-}

