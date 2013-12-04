{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeOperators #-}

module BioInf.Alignment.DnaProtein where

import Data.Array.Repa.Index
import qualified Data.Vector.Fusion.Stream.Monadic as M
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

data Monad m => SigDnaPro {-Monad-} m {-NT-} nt {-T-} a c e = SigDnaPro
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
  ,h :: M.Stream m nt -> m nt}

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

-- |

algScore :: SigDnaPro IO Int Char Char ()
algScore = SigDnaPro
  { lcldel    = id
  , nilnil    = const 0
  , delamino  = \z (Z:.():.a)                          -> z + error "delamino"
  , rf1amino  = \z (Z:.c1:.a)  (Z:.c2:.())             -> z + error "rf1amino"
  , rf1del    = \z (Z:.c1:.()) (Z:.c2:.())             -> z + error "rf1del"
  , rf2amino  = \z (Z:.c:.a)                           -> z + error "rf2amino"
  , rf2del    = \z (Z:.c:.())                          -> z + error "rf2del"
  , stayamino = \z (Z:.c1:.a)  (Z:.c2:.()) (Z:.c3:.()) -> z + error "stayamino"
  , staydel   = \z (Z:.c1:.()) (Z:.c2:.()) (Z:.c3:.()) -> z + error "staydel"
  , eatdel    = \z (Z:.c:.())                          -> z + error "eatdel"
  , h         = M.foldl' max (-999999)
  }
{-# INLINE algScore #-}

-- |

dnaProtein dna protein = ( _RP ! (Z:.pointL 0 nD:.pointL 0 nP), [] ) where
  (_F0P, _F1P, _F2P, _LP, _RP) = unsafePerformIO $ fillTables dna protein
  nD = V.length dna
  nP = V.length protein

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
               (chr dna :: GChr Char Char) (chr protein :: GChr Char Char)
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
    let i = (Z:.pointL 0 nD:.pointL 0 nP)
    f_f0p i >>= writeM tf0p i
    f_f1p i >>= writeM tf1p i
    f_f2p i >>= writeM tf2p i
    f_lp  i >>= writeM tlp  i
    f_rp  i >>= writeM trp  i

