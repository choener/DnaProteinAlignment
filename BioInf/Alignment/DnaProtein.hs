
module BioInf.Alignment.DnaProtein where

import ADP.Fusion
import ADP.Fusion.Multi

data SigDnaPro {-Monad-} m {-NT-} nt {-T-} c a {-E-} eps e = SigDnaPro
  {delamino :: nt -> (Z:.eps:.a) -> nt
  ,nilnil :: (Z:.e:.e) -> nt
  ,rf1amino :: nt -> (Z:.c:.a) -> (Z:.c:.eps) -> nt
  ,rf1del :: nt -> (Z:.c:.eps) -> (Z:.c:.eps) -> nt
  ,rf2amino :: nt -> (Z:.c:.a) -> nt
  ,rf2del :: nt -> (Z:.c:.eps) -> nt
  ,stayamino :: nt -> (Z:.c:.a) -> (Z:.c:.eps) -> (Z:.c:.eps) -> nt
  ,staydel :: nt -> (Z:.c:.eps) -> (Z:.c:.eps) -> (Z:.c:.eps) -> nt
  ,h :: M.Stream m nt -> nt}
grammarDnaPro SigDnaPro{..} {-NT-} _F0P _F1P _F2P {-T-} a c  {-E-} eps e =
  ((_F0P
   ,delamino <<< _F0P % (Z:.eps:.a)
     ||| nilnil <<< (Z:.e:.e)
     ||| rf1amino <<< _F1P % (Z:.c:.a) % (Z:.c:.eps)
     ||| rf1del <<< _F1P % (Z:.c:.eps) % (Z:.c:.eps)
     ||| rf2amino <<< _F2P % (Z:.c:.a)
     ||| rf2del <<< _F2P % (Z:.c:.eps)
     ||| stayamino <<< _F0P % (Z:.c:.a) % (Z:.c:.eps) % (Z:.c:.eps)
     ||| staydel <<< _F0P % (Z:.c:.eps) % (Z:.c:.eps) % (Z:.c:.eps) ... h)
  ,(_F1P
   ,delamino <<< _F1P % (Z:.eps:.a)
     ||| nilnil <<< (Z:.e:.e)
     ||| rf1amino <<< _F2P % (Z:.c:.a) % (Z:.c:.eps)
     ||| rf1del <<< _F2P % (Z:.c:.eps) % (Z:.c:.eps)
     ||| rf2amino <<< _F0P % (Z:.c:.a)
     ||| rf2del <<< _F0P % (Z:.c:.eps)
     ||| stayamino <<< _F1P % (Z:.c:.a) % (Z:.c:.eps) % (Z:.c:.eps)
     ||| staydel <<< _F1P % (Z:.c:.eps) % (Z:.c:.eps) % (Z:.c:.eps) ... h)
  ,(_F2P
   ,delamino <<< _F2P % (Z:.eps:.a)
     ||| nilnil <<< (Z:.e:.e)
     ||| rf1amino <<< _F0P % (Z:.c:.a) % (Z:.c:.eps)
     ||| rf1del <<< _F0P % (Z:.c:.eps) % (Z:.c:.eps)
     ||| rf2amino <<< _F1P % (Z:.c:.a)
     ||| rf2del <<< _F1P % (Z:.c:.eps)
     ||| stayamino <<< _F2P % (Z:.c:.a) % (Z:.c:.eps) % (Z:.c:.eps)
     ||| staydel <<< _F2P % (Z:.c:.eps) % (Z:.c:.eps) % (Z:.c:.eps) ... h))
