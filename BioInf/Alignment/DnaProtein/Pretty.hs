{-# LANGUAGE TypeOperators #-}

{-# OPTIONS_GHC -fno-liberate-case -O0 #-}

module BioInf.Alignment.DnaProtein.Pretty (backtrack,PPS(..),isLOC,isPPS) where

import           Control.Monad
import           Data.Array.Repa.Index
import           Data.Vector.Fusion.Util
import qualified Data.Array.IArray as A
import qualified Data.Vector.Fusion.Stream.Monadic
import qualified Data.Vector.Fusion.Stream.Monadic as M
import qualified Data.Vector.Unboxed as V
import           System.IO.Unsafe
import qualified Data.Array.IArray as A

import Biobase.Primary
import Data.PrimitiveArray as PA hiding (map,fromList)
import Data.PrimitiveArray.Zero as PA hiding (map,fromList)
import Data.Array.Repa.Index.Points
import ADP.Fusion.Multi
import ADP.Fusion.Multi.Classes
import ADP.Fusion.Multi.Empty
import ADP.Fusion hiding (empty)
import ADP.Fusion.Chr
import ADP.Fusion.Table
import ADP.Fusion.Empty hiding (empty)
import ADP.Fusion.None
import Biobase.SubstMatrix
import Data.DList (empty,snoc,DList,append,fromList)

import BioInf.Alignment.DnaProtein.Common
import BioInf.Alignment.DnaProtein.Score (algScore)



algPretty :: Monad m => SigDnaPro m (DList Char,DList Char) (M.Stream m (DList Char,DList Char)) Char Nuc ()
algPretty = SigDnaPro
  { lcldel = id
  , nilnil = const (empty,empty)
  , delamino = \(x,y) (Z:.():.a) -> ( x `appendL` "---"
                                    , y `snocA` a
                                    )
  , rf1amino = \(x,y) (Z:.c1:.a) (Z:.c2:.()) -> ( appendNL x [c1, c2] "-"
                                                , y `snocA` a
                                                )
  , rf1del   = \(x,y) (Z:.c1:.()) (Z:.c2:.()) -> ( appendNL x [c1,c2] "-"
                                                 , y `snocA` '-'
                                                 )
  , rf2amino = \(x,y) (Z:.c:.a) -> ( appendNL x [c] "--"
                                   , y `snocA` a
                                   )
  , rf2del   = \(x,y) (Z:.c:.()) -> ( appendNL x [c] "--"
                                    , y `snocA` '-'
                                    )
  , stayamino = \(x,y) (Z:.c1:.a) (Z:.c2:.()) (Z:.c3:.()) -> ( appendNL x [c1,c2,c3] ""
                                                             , y `snocA` a
                                                             )
  , staydel = \(x,y) (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> ( appendNL x [c1,c2,c3] ""
                                                            , y `snocA` '-'
                                                            )
  , eatdel = \(x,y) (Z:.c:.()) -> ( x `snoc` fromNuc c
                                  , y `snoc` '.'
                                  )
  , h = return . id
  } where
      appendNL l ns cs = (l `append` (fromList $ map fromNuc ns)) `append` (fromList cs)
      appendL l r = l `append` fromList r
      snocA   y a = y `append` fromList [a,' ',' ']

data PPS
  = PPS ![Nuc] ![Char] !Int
  | FRS ![Nuc] ![Char] !Int
  | LOC !Nuc           !Int
  deriving (Eq,Ord,Show)

isLOC (LOC _   _) = True
isLOC _           = False

isPPS (LOC _   _) = False
isPPS _           = True

algPPscore :: Monad m
  => Nuc3SubstMatrix
  -> Nuc2SubstMatrix
  -> Nuc1SubstMatrix
  -> Int -> Int -> Int -> Int -> Int -> Int
  -> SigDnaPro m (DList PPS) (M.Stream m (DList PPS)) Char Nuc ()
algPPscore n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS = SigDnaPro
  { lcldel = id
  , nilnil = const empty
  , delamino = \x (Z:.():.a) -> x `snoc` PPS [] [a] 0
  , rf1amino = \x (Z:.c1:.a) (Z:.c2:.()) -> x `snoc` FRS [c1,c2] [a] (n2m A.! (c1,c2,a))
  , rf1del   = \x (Z:.c1:.()) (Z:.c2:.()) -> x `snoc` PPS [c1,c2] [] 0
  , rf2amino = \x (Z:.c:.a) -> x `snoc` FRS [c] [a] (n1m A.! (c,a))
  , rf2del   = \x (Z:.c:.()) -> x `snoc` PPS [c] [] 0
  , stayamino = \x (Z:.c1:.a) (Z:.c2:.()) (Z:.c3:.()) -> x `snoc` PPS [c1,c2,c3] [a] (n3m A.! (c1,c2,c3,a))
  , staydel = \x (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> x `snoc` PPS [c1,c2,c3] [] 0
  , eatdel = \x (Z:.c:.()) -> x `snoc` LOC c 0
  , h = return . id
  } where

-- |

(<**) f g = SigDnaPro
  {delamino = \(a,aN) b -> (_Fdelamino a b, aN >>= return . M.map (\a -> _Gdelamino a b))
  ,lcldel = \(a,aN) -> (_Flcldel a, aN >>= return . M.map (\a -> _Glcldel a))
  ,nilnil = \a -> (_Fnilnil a, return $ M.singleton (_Gnilnil a))
  ,rf1amino = \(a,aN) b c -> (_Frf1amino a b c, aN >>= return . M.map (\a -> _Grf1amino a b c))
  ,rf1del = \(a,aN) b c -> (_Frf1del a b c, aN >>= return . M.map (\a -> _Grf1del a b c))
  ,rf2amino = \(a,aN) b -> (_Frf2amino a b, aN >>= return . M.map (\a -> _Grf2amino a b))
  ,rf2del = \(a,aN) b -> (_Frf2del a b, aN >>= return . M.map (\a -> _Grf2del a b))
  ,stayamino = \(a,aN) b c d -> (_Fstayamino a b c d, aN >>= return . M.map (\a -> _Gstayamino a b c d))
  ,staydel = \(a,aN) b c d -> (_Fstaydel a b c d, aN >>= return . M.map (\a -> _Gstaydel a b c d))
  ,eatdel = \(a,aN) b -> (_Featdel a b, aN >>= return . M.map (\a -> _Geatdel a b))
  ,h = \xs -> do
    hfs <- _Fh . Data.Vector.Fusion.Stream.Monadic.map fst $ xs
    let phfs = Data.Vector.Fusion.Stream.Monadic.concatMapM snd
             . Data.Vector.Fusion.Stream.Monadic.filter ((hfs==) . fst) $ xs
    _Gh phfs}
  where
    _Fdelamino = delamino f
    _Flcldel = lcldel f
    _Fnilnil = nilnil f
    _Frf1amino = rf1amino f
    _Frf1del = rf1del f
    _Frf2amino = rf2amino f
    _Frf2del = rf2del f
    _Fstayamino = stayamino f
    _Fstaydel = staydel f
    _Featdel = eatdel f
    _Fh = h f
    _Gdelamino = delamino g
    _Glcldel = lcldel g
    _Gnilnil = nilnil g
    _Grf1amino = rf1amino g
    _Grf1del = rf1del g
    _Grf2amino = rf2amino g
    _Grf2del = rf2del g
    _Gstayamino = stayamino g
    _Gstaydel = staydel g
    _Geatdel = eatdel g
    _Gh = h g

backtrack
  :: Nuc3SubstMatrix
  -> Nuc2SubstMatrix
  -> Nuc1SubstMatrix
  -> Int -> Int -> Int -> Int -> Int -> Int
  -> V.Vector Nuc
  -> V.Vector Char
  -> ( PPT, PPT, PPT, PPT, PPT )
  -> [DList PPS] -- [(DList Doc,DList Doc)]
backtrack n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS dna protein (f0p, f1p, f2p, lp, rp) =
  unId . M.toList . unId . f_rp $ Z:.pointL 0 nD:.pointL 0 nP where
  nD = V.length dna
  nP = V.length protein
  ( (_,f_f0p), (_,f_f1p), (_,f_f2p), (_,f_lp ), (_,f_rp ) ) = grammarDnaPro
        (algScore n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS
        <**
        algPPscore n3m n2m n1m insertAA deleteAA rf1S rf1delS rf2S rf2delS
        )
        (btTbl (Z:.EmptyT:.EmptyT) f0p (f_f0p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) f1p (f_f1p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) f2p (f_f2p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) lp  (f_lp  :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) rp  (f_rp  :: FunT) :: BtPPT)
        (chr protein)
        (chr dna)
        Empty
{-# NOINLINE backtrack #-}

type BtPPT = DefBtTbl Id (Z:.PointL:.PointL) Int (DList PPS) -- (DList Doc,DList Doc)
type FunT = (Z:.PointL:.PointL) -> Id (M.Stream Id (DList PPS)) -- (DList Doc,DList Doc))

