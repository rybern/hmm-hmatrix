{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
module Math.HiddenMarkovModel (
   T(..), Distr.State, state,
   Discrete, DiscreteTrained,
   Gaussian, GaussianTrained,
   uniform,
   generate,
   probabilitySequence,
   Normalized.logLikelihood,
   Normalized.reveal,

   Trained(..),
   trainSupervised,
   Normalized.trainUnsupervised,
   mergeTrained, finishTraining, trainMany,
   deviation,

   toCSV,
   fromCSV,
   ) where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel.Normalized as Normalized
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Private
          (T(..), Trained(..), mergeTrained, toCells, parseCSV)
import Math.HiddenMarkovModel.Distribution (State(State))
import Math.HiddenMarkovModel.Utility
          (randomItemProp, normalizeProb, attachOnes)

import qualified Numeric.LinearAlgebra.Algorithms as Algo
import qualified Numeric.Container as NC
import qualified Data.Packed.Matrix as Matrix
import qualified Data.Packed.Vector as Vector
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import qualified Text.CSV.Lazy.String as CSV

import qualified System.Random as Rnd

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.State as MS
import qualified Control.Monad.HT as Monad

import qualified Data.NonEmpty as NonEmpty
import qualified Data.Array as Array
import Data.Traversable (Traversable, mapAccumL)
import Data.Foldable (Foldable)
import Data.Array (accumArray)



state :: Int -> State
state = State


type DiscreteTrained prob symbol = Trained (Distr.DiscreteTrained prob symbol) prob
type Discrete prob symbol = T (Distr.Discrete prob symbol) prob

type GaussianTrained a = Trained (Distr.GaussianTrained a) a
type Gaussian a = T (Distr.Gaussian a) a


{- |
Create a model with uniform probabilities
for initial vector and transition matrix
given a distribution for the emissions.
You can use this as a starting point for 'Normalized.trainUnsupervised'.
-}
uniform ::
   (Distr.Info distr, Distr.Probability distr ~ prob) =>
   distr -> T distr prob
uniform distr =
   let n = Distr.numberOfStates distr
       c = recip $ fromIntegral n
   in  Cons {
          initial = NC.constant c n,
          transition = NC.konst c (n,n),
          distribution = distr
       }


probabilitySequence ::
   (Traversable f, Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> f (State, emission) -> f prob
probabilitySequence hmm =
   snd
   .
   mapAccumL
      (\index (State s, e) ->
         (NC.atIndex (transition hmm) . flip (,) s,
          index s * Distr.emissionStateProb (distribution hmm) e (State s)))
      (NC.atIndex (initial hmm))

generate ::
   (Rnd.RandomGen g, Ord prob, Rnd.Random prob,
    Distr.Generate distr, Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> g -> [emission]
generate hmm =
   MS.evalState $
   flip MS.evalStateT (initial hmm) $
   Monad.repeat $ MS.StateT $ \v0 -> do
      s <- randomItemProp $ zip [0..] (Vector.toList v0)
      x <- Distr.generate (distribution hmm) (State s)
      return (x, takeColumn s $ transition hmm)

takeColumn :: (Matrix.Element a) => Int -> Matrix a -> Vector a
takeColumn n  =  Matrix.flatten . Matrix.extractColumns [n]



{- |
Contribute a manually labeled emission sequence to a HMM training.
-}
trainSupervised ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   Int -> NonEmpty.T [] (State, emission) -> Trained tdistr prob
trainSupervised n xs =
   let getState (State s, _x) = s
   in  Trained {
          trainedInitial = NC.assoc n 0 [(getState (NonEmpty.head xs), 1)],
          trainedTransition =
             Matrix.trans $ NC.accum (NC.konst 0 (n,n)) (+) $
             attachOnes $ NonEmpty.mapAdjacent (,) $ fmap getState xs,
          trainedDistribution =
             Distr.accumulateEmissions $ map attachOnes $ Array.elems $
             accumArray (flip (:)) [] (State 0, State (n-1)) $ NonEmpty.flatten xs
       }

finishTraining ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob) =>
   Trained tdistr prob -> T distr prob
finishTraining hmm =
   Cons {
      initial = normalizeProb $ trainedInitial hmm,
      transition =
         Matrix.fromColumns $ map normalizeProb $
         Matrix.toColumns $ trainedTransition hmm,
      distribution = Distr.normalize $ trainedDistribution hmm
   }

trainMany ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob,
    Foldable f) =>
   (trainingData -> Trained tdistr prob) ->
   NonEmpty.T f trainingData -> T distr prob
trainMany train =
   finishTraining . NonEmpty.foldl1Map mergeTrained train





{- |
Compute maximum deviation between initial and transition probabilities.
You can use this as abort criterion for unsupervised training.
We omit computation of differences between the emission probabilities.
This simplifies matters a lot and
should suffice for defining an abort criterion.
-}
deviation ::
   (Algo.Field prob, Ord prob) => T distr prob -> T distr prob -> prob
deviation hmm0 hmm1 =
   deviationVec (initial hmm0) (initial hmm1)
   `max`
   deviationVec (transition hmm0) (transition hmm1)

deviationVec ::
   (Ord a, NC.Container c a) =>
   c a -> c a -> a
deviationVec x y =
   let d = NC.sub x y
   in  NC.maxElement d `max` negate (NC.minElement d)


toCSV ::
   (Distr.CSV distr, Algo.Field prob, Show prob) =>
   T distr prob -> String
toCSV hmm =
   CSV.ppCSVTable $ snd $ CSV.toCSVTable $ HMMCSV.padTable "" $
   toCells hmm

fromCSV ::
   (Distr.CSV distr, Algo.Field prob, Read prob) =>
   String -> ME.Exceptional String (T distr prob)
fromCSV =
   MS.evalStateT parseCSV . map HMMCSV.fixShortRow . CSV.parseCSV
