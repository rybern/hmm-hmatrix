{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{- |
Counterparts to functions in "Math.HiddenMarkovModel.Private"
that normalize interim results.
We need to do this in order to prevent
to round very small probabilities to zero.
-}
module Math.HiddenMarkovModel.Normalized where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import Math.HiddenMarkovModel.Private
          (T(..), Trained(..), emission, matrixMaxMul, sumTransitions)
import Math.HiddenMarkovModel.Distribution (State(State))
import Math.HiddenMarkovModel.Utility (normalizeFactor, normalizeProb)

import qualified Numeric.Container as NC
import qualified Data.Packed.Development as Dev
import qualified Data.Packed.Vector as Vector
import Numeric.Container ((<>))
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import qualified Control.Functor.HT as Functor

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Foldable as Fold
import qualified Data.List as List
import Data.Traversable (Traversable, mapAccumL)
import Data.Tuple.HT (mapFst, mapSnd, swap)

posterior ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f, NonEmptyC.Zip f, NonEmptyC.Reverse f) =>
   T distr prob ->
   NonEmpty.T f emission ->
   NonEmpty.T f (Vector prob)
posterior hmm xs = fmap normalizeProb $ zetaFromAlphaBeta alphas betas
  where (alphas, betas) = alphaBeta hmm xs


{- |
Logarithm of the likelihood to observe the given sequence.
We return the logarithm because the likelihood can be so small
that it may be rounded to zero in the choosen number type.
-}
logLikelihood ::
   (Distr.EmissionProb distr, Floating prob,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob -> NonEmpty.T f emission -> prob
logLikelihood hmm = Fold.sum . fmap (log . fst) . alpha hmm

alpha ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob ->
   NonEmpty.T f emission -> NonEmpty.T f (prob, Vector prob)
alpha hmm (NonEmpty.Cons x xs) =
   let normMulEmiss y = normalizeFactor . NC.mul (emission hmm y)
   in  NonEmpty.scanl
          (\(_,alphai) xi -> normMulEmiss xi (transition hmm <> alphai))
          (normMulEmiss x (initial hmm))
          xs

beta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f, NonEmptyC.Reverse f) =>
   T distr prob ->
   f (prob, emission) -> NonEmpty.T f (Vector prob)
beta hmm =
   nonEmptyScanr
      (\(ci,xi) betai ->
         NC.scale (recip ci) $ NC.mul (emission hmm xi) betai <> transition hmm)
      (NC.constant 1 (NC.dim $ initial hmm))

alphaBeta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f, NonEmptyC.Zip f, NonEmptyC.Reverse f) =>
   T distr prob ->
   NonEmpty.T f emission ->
   (NonEmpty.T f (prob, Vector prob), NonEmpty.T f (Vector prob))
alphaBeta hmm xs =
   let calphas = alpha hmm xs
   in  (calphas,
        beta hmm $ NonEmpty.tail $ NonEmptyC.zip (fmap fst calphas) xs)


xiFromAlphaBeta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f, NonEmptyC.Zip f) =>
   T distr prob ->
   NonEmpty.T f emission ->
   NonEmpty.T f (prob, Vector prob) ->
   NonEmpty.T f (Vector prob) ->
   f (Matrix prob)
xiFromAlphaBeta hmm xs calphas betas =
   let (cs,alphas) = Functor.unzip calphas
   in  NonEmptyC.zipWith4
          (\x alpha0 c1 beta1 ->
             NC.scale (recip c1) $
             NC.mul
                (NC.outer (NC.mul (emission hmm x) beta1) alpha0)
                (transition hmm))
          (NonEmpty.tail xs)
          (NonEmpty.init alphas)
          (NonEmpty.tail cs)
          (NonEmpty.tail betas)

zetaFromAlphaBeta ::
   (NC.Container Vector prob, NonEmptyC.Zip f) =>
   NonEmpty.T f (prob, Vector prob) ->
   NonEmpty.T f (Vector prob) ->
   NonEmpty.T f (Vector prob)
zetaFromAlphaBeta calphas betas =
   NonEmptyC.zipWith (NC.mul . snd) calphas betas


{- |
Reveal the state sequence
that led most likely to the observed sequence of emissions.
It is found using the Viterbi algorithm.
-}
reveal ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f, NonEmptyC.Reverse f) =>
   T distr prob -> NonEmpty.T f emission -> NonEmpty.T f State
reveal hmm (NonEmpty.Cons x xs) =
   fmap State $
   uncurry (nonEmptyScanr Dev.at') $
   mapFst NC.maxIndex $
   mapAccumL
      (\alphai xi ->
         swap $ mapSnd (NC.mul (emission hmm xi)) $
         matrixMaxMul (transition hmm) $ normalizeProb alphai)
      (NC.mul (emission hmm x) (initial hmm)) xs


{- |
Variant of NonEmpty.scanr with less stack consumption.
-}
nonEmptyScanr ::
   (Traversable f, NonEmptyC.Reverse f) =>
   (a -> b -> b) -> b -> f a -> NonEmpty.T f b
nonEmptyScanr f x =
   NonEmptyC.reverse . NonEmpty.scanl (flip f) x . NonEmptyC.reverse


{- |
Consider a superposition of all possible state sequences
weighted by the likelihood to produce the observed emission sequence.
Now train the model with respect to all of these sequences
with respect to the weights.
This is done by the Baum-Welch algorithm.
-}
trainUnsupervised ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> NonEmpty.T [] emission -> Trained tdistr prob
trainUnsupervised hmm xs =
   let (alphas, betas) = alphaBeta hmm xs
       zetas = zetaFromAlphaBeta alphas betas

   in  Trained {
          trainedInitial = NonEmpty.head zetas,
          trainedTransition =
             sumTransitions hmm $ xiFromAlphaBeta hmm xs alphas betas,
          trainedDistribution =
             Distr.accumulateEmissions $ map (zip (NonEmpty.flatten xs)) $
             List.transpose $ map Vector.toList $ NonEmpty.flatten zetas
       }
