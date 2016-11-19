module Math.HiddenMarkovModel.Test where

import qualified Math.HiddenMarkovModel as HMM
import qualified Math.HiddenMarkovModel.Normalized as Normalized
import qualified Math.HiddenMarkovModel.Private as Priv
import qualified Math.HiddenMarkovModel.Distribution as Distr

import qualified Numeric.Container as NC
import qualified Data.Packed.Matrix as Matrix
import qualified Data.Packed.Vector as Vector
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import qualified System.Random as Rnd

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Traversable as Trav
import qualified Data.Foldable as Fold
import qualified Data.List as List
import qualified Data.Map as Map
import Data.NonEmpty ((!:))


hmm :: HMM.Discrete Double Char
hmm =
   HMM.Cons {
      HMM.initial = Vector.fromList [0.1, 0.2, 0.3, 0.4],
      HMM.transition =
         Matrix.fromLists $
            [0.7, 0.1, 0.0, 0.2] :
            [0.1, 0.6, 0.1, 0.0] :
            [0.1, 0.2, 0.7, 0.0] :
            [0.1, 0.1, 0.2, 0.8] :
            [],
      HMM.distribution =
         Distr.Discrete $ Map.fromList $
            ('a', Vector.fromList [1,0,0,0]) :
            ('b', Vector.fromList [0,1,0,1]) :
            ('c', Vector.fromList [0,0,1,0]) :
            []
   }


sequ :: NonEmpty.T [] Char
sequ = 'a' !: take 20 (HMM.generate hmm (Rnd.mkStdGen 42))

possibleStates :: Char -> [Distr.State]
possibleStates c =
   map Distr.State $ List.findIndices id $
   map
      (\p ->
         case p of
            0 -> False
            1 -> True
            _ -> error "invalid emission probability (must be 0 or 1)") $
   Vector.toList $
   Map.findWithDefault (error "invalid character") c $
   case HMM.distribution hmm of Distr.Discrete m -> m

{- |
Should all be equal.
-}
sequLikelihood :: ((Double, Double), Double, Double, NonEmpty.T [] Double)
sequLikelihood =
   ((Priv.forward hmm sequ, Priv.backward hmm sequ),
    exp $ Normalized.logLikelihood hmm sequ,
    sum $
       map (NonEmpty.product . HMM.probabilitySequence hmm) $
       Trav.mapM (\c -> map (flip (,) c) $ possibleStates c) sequ,
    NonEmptyC.zipWith NC.dot
       (Priv.alpha hmm sequ)
       (Priv.beta hmm $ NonEmpty.tail sequ))

{- |
Should all be one.
-}
sequLikelihoodNormalized :: NonEmpty.T [] Double
sequLikelihoodNormalized =
   let (calphas,betas) = Normalized.alphaBeta hmm sequ
   in  NonEmptyC.zipWith NC.dot (fmap snd calphas) betas


{- |
Lists should be equal, but the first list contains one less element.
-}
zetas ::
   ([Vector Double],
    NonEmpty.T [] (Vector Double),
    NonEmpty.T [] (Vector Double))
zetas =
   let (recipLikelihood, alphas, betas) = Priv.alphaBeta hmm sequ
   in  (Priv.zetaFromXi hmm $
           Priv.xiFromAlphaBeta hmm recipLikelihood sequ alphas betas,
        Priv.zetaFromAlphaBeta recipLikelihood alphas betas,
        uncurry Normalized.zetaFromAlphaBeta $
        Normalized.alphaBeta hmm sequ)

{- |
Quick test of zetas - result should be @(True, very small, very small)@.
-}
zetasDiff :: (Bool, Double, Double)
zetasDiff =
   case zetas of
      (z0,z1,z2) ->
         (length z0 == length (NonEmpty.tail z1) &&
          length z0 == length (NonEmpty.tail z2),
          maximum $ map NC.normInf $ zipWith NC.sub z0 $ NonEmpty.init z1,
          NonEmpty.maximum $ fmap NC.normInf $ NonEmptyC.zipWith NC.sub z1 z2)

{- |
Lists should be equal
-}
xis :: ([Matrix Double], [Matrix Double])
xis =
   let (recipLikelihood, alphas, betas) = Priv.alphaBeta hmm sequ
   in  (Priv.xiFromAlphaBeta hmm recipLikelihood sequ alphas betas,
        uncurry (Normalized.xiFromAlphaBeta hmm sequ) $
        Normalized.alphaBeta hmm sequ)

{- |
Quick test of xis - result should be @(True, very small)@.
-}
xisDiff :: (Bool, Double)
xisDiff =
   case xis of
      (x0,x1) ->
         (length x0 == length x1,
          maximum $ map (NC.normInf . Matrix.flatten) $ zipWith NC.sub x0 x1)


reveal :: Bool
reveal =
   Normalized.reveal hmm sequ == Priv.reveal hmm sequ


trainUnsupervised ::
   (HMM.DiscreteTrained Double Char,
    HMM.DiscreteTrained Double Char)
trainUnsupervised =
   (Priv.trainUnsupervised hmm sequ,
    Normalized.trainUnsupervised hmm sequ)

trainUnsupervisedDiff :: (Double, Double, (Bool, Double))
trainUnsupervisedDiff =
   case trainUnsupervised of
      (hmm0,hmm1) ->
         (NC.normInf $ Matrix.flatten $ NC.sub
             (Priv.trainedTransition hmm0) (Priv.trainedTransition hmm1),
          NC.normInf $ NC.sub
             (Priv.trainedInitial hmm0) (Priv.trainedInitial hmm1),
          case (Priv.trainedDistribution hmm0, Priv.trainedDistribution hmm1) of
             (Distr.DiscreteTrained m0, Distr.DiscreteTrained m1) ->
                (Map.size m0 == Map.size m1,
                 Fold.maximum $ fmap NC.normInf $
                    Map.intersectionWith NC.sub m0 m1))


nonEmptyScanr :: Int -> [Int] -> Bool
nonEmptyScanr x xs =
   Normalized.nonEmptyScanr (-) x xs == NonEmpty.scanr (-) x xs
