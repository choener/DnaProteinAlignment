{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeOperators #-}

module BioInf.Alignment.DnaProtein where

import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream.Monadic as M
import qualified Data.Vector.Fusion.Stream.Monadic
import Data.Vector.Fusion.Util
import System.IO.Unsafe
import qualified Data.Vector.Unboxed as V
import Control.Monad

import Data.PrimitiveArray as PA
import Data.PrimitiveArray.Zero as PA
import Data.Array.Repa.Index.Points
import ADP.Fusion.Multi
import ADP.Fusion.Multi.Classes
import ADP.Fusion.Multi.Empty
import ADP.Fusion
import ADP.Fusion.Chr
import ADP.Fusion.Table
import ADP.Fusion.Empty
import ADP.Fusion.None

import Debug.Trace


{-
data Monad m => SigDnaPro {-Monad-} m {-NT-} nt r {-T-} a c e = SigDnaPro
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
  ,h :: M.Stream m nt -> m r}

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
-}


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


-- |

algScore :: Monad m => SigDnaPro m Int Int Char Char ()
algScore = SigDnaPro
  { lcldel    = id
  , nilnil    = const 0
  , delamino  = \z (Z:.():.a)                          -> z -1
  , rf1amino  = \z (Z:.c1:.a)  (Z:.c2:.())             -> z -1
  , rf1del    = \z (Z:.c1:.()) (Z:.c2:.())             -> z -99999
  , rf2amino  = \z (Z:.c:.a)                           -> z -2
  , rf2del    = \z (Z:.c:.())                          -> z -2
  , stayamino = \z (Z:.c1:.a)  (Z:.c2:.()) (Z:.c3:.()) -> z + if c1==a && c2==a && c3==a then 5 else -1
  , staydel   = \z (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> z -1
  , eatdel    = \z (Z:.c:.())                          -> z
  , h         = M.foldl' max (-999999)
  }
{-# INLINE algScore #-}

algPretty :: Monad m => SigDnaPro m (String,String) (M.Stream m (String,String)) Char Char ()
algPretty = SigDnaPro
  { lcldel = id
  , nilnil = const ([],[])
  , delamino = \(x,y) (Z:.():.a) -> (x++" ",y++[a])
  , rf1amino = \(x,y) (Z:.c1:.a) (Z:.c2:.()) -> (x++[c1,c2],y++[a,' '])
  , rf1del   = \(x,y) (Z:.c1:.()) (Z:.c2:.()) -> (x++[c1,c2],y++"  ")
  , rf2amino = \(x,y) (Z:.c:.a) -> (x++[c],y++[a])
  , rf2del   = \(x,y) (Z:.c:.()) -> (x++[c],y++" ")
  , stayamino = \(x,y) (Z:.c1:.a) (Z:.c2:.()) (Z:.c3:.()) -> (x++[c1,c2,c3],y++[a,' ',' '])
  , staydel = \(x,y) (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> (x++[c1,c2,c3],y++"   ")
  , eatdel = \(x,y) (Z:.c:.()) -> (x++[c],y++" ")
  , h = return . id
  }

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


-- |

dnaProtein dna' protein' = ( _RP ! (Z:.pointL 0 nD:.pointL 0 nP), bt ) where
  ts@(_F0P, _F1P, _F2P, _LP, _RP) = unsafePerformIO $ fillTables dna protein
  nD = V.length dna
  nP = V.length protein
  dna = V.fromList dna'
  protein = V.fromList protein'
  bt = backtrack dna protein ts

-- |

type PPT = PA.Unboxed (Z:.PointL:.PointL) Int
type BtPPT = DefBtTbl Id (Z:.PointL:.PointL) Int (String,String)
type FunT = (Z:.PointL:.PointL) -> Id (M.Stream Id (String,String))

backtrack
  :: V.Vector Char
  -> V.Vector Char
  -> ( PPT, PPT, PPT, PPT, PPT )
  -> [(String,String)]
backtrack dna protein (f0p, f1p, f2p, lp, rp) = unId . M.toList . unId . f_rp $ Z:.pointL 0 nD:.pointL 0 nP where
  nD = V.length dna
  nP = V.length protein
  ( (_,f_f0p), (_,f_f1p), (_,f_f2p), (_,f_lp ), (_,f_rp ) ) = grammarDnaPro
        (algScore <** algPretty)
        (btTbl (Z:.EmptyT:.EmptyT) f0p (f_f0p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) f1p (f_f1p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) f2p (f_f2p :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) lp  (f_lp  :: FunT) :: BtPPT)
        (btTbl (Z:.EmptyT:.EmptyT) rp  (f_rp  :: FunT) :: BtPPT)
        (chr protein)
        (chr dna)
        Empty

-- |

type Tbl = PA.Unboxed (Z:.PointL:.PointL) Int

fillTables :: V.Vector Char -> V.Vector Char
  -> IO ( Tbl, Tbl, Tbl, Tbl, Tbl )
fillTables dna protein = do
  let nD = V.length dna
  let nP = V.length protein
  f0p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) 0
  f1p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) 0
  f2p' <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) 0
  lp'  <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) 0
  rp'  <- newWithM (Z:.pointL 0 0:.pointL 0 0) (Z:.pointL 0 nD:.pointL 0 nP) 0
  fillFive $ grammarDnaPro
               algScore
               (mTbl (Z:.EmptyT:.EmptyT) f0p')
               (mTbl (Z:.EmptyT:.EmptyT) f1p')
               (mTbl (Z:.EmptyT:.EmptyT) f2p')
               (mTbl (Z:.EmptyT:.EmptyT) lp')
               (mTbl (Z:.EmptyT:.EmptyT) rp')
               (chr protein :: GChr Char Char) (chr dna :: GChr Char Char)
               Empty
  f0p <- freeze f0p'
  f1p <- freeze f1p'
  f2p <- freeze f2p'
  lp  <- freeze lp'
  rp  <- freeze rp'
  return (f0p,f1p,f2p,lp,rp)

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

