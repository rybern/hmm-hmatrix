{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
module Math.HiddenMarkovModel.Private where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Distribution (State(State))

import qualified Numeric.LinearAlgebra.Algorithms as Algo
import qualified Numeric.LinearAlgebra.Util as LinAlg
import qualified Numeric.Container as NC
import qualified Data.Packed.Development as Dev
import qualified Data.Packed.Matrix as Matrix
import qualified Data.Packed.Vector as Vector
import Numeric.Container ((<>))
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import Control.DeepSeq (NFData, rnf)
import Foreign.Storable (Storable)

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Semigroup as Sg
import qualified Data.List as List
import Data.Traversable (Traversable, mapAccumL)
import Data.Tuple.HT (mapPair, mapFst, mapSnd, swap)


{- |
A Hidden Markov model consists of a number of (hidden) states
and a set of emissions.
There is a vector for the initial probability of each state
and a matrix containing the probability for switching
from one state to another one.
The 'distribution' field points to probability distributions
that associate every state with emissions of different probability.
Famous distribution instances are discrete and Gaussian distributions.
See "Math.HiddenMarkovModel.Distribution" for details.

The transition matrix is transposed
with respect to popular HMM descriptions.
But I think this is the natural orientation, because this way
you can write \"transition matrix times probability column vector\".

The type has two type parameters,
although the one for the distribution would be enough.
However, replacing @prob@ by @Distr.Probability distr@
would prohibit the derived Show and Read instances.
-}
data T distr prob =
   Cons {
      initial :: Vector prob,
      transition :: Matrix prob,
      distribution :: distr
   }
   deriving (Show, Read)

instance
   (NFData distr, NFData prob, Storable prob) =>
      NFData (T distr prob) where
   rnf hmm = rnf (initial hmm, transition hmm, distribution hmm)


emission ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> emission -> Vector prob
emission  =  Distr.emissionProb . distribution


forward ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob -> NonEmpty.T f emission -> prob
forward hmm = NC.sumElements . NonEmpty.last . alpha hmm

alpha ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob ->
   NonEmpty.T f emission -> NonEmpty.T f (Vector prob)
alpha hmm (NonEmpty.Cons x xs) =
   NonEmpty.scanl
      (\alphai xi -> NC.mul (emission hmm xi) (transition hmm <> alphai))
      (NC.mul (emission hmm x) (initial hmm))
      xs


backward ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob -> NonEmpty.T f emission -> prob
backward hmm (NonEmpty.Cons x xs) =
   NC.sumElements $
   NC.mul (initial hmm) $
   NC.mul (emission hmm x) $
   NonEmpty.head $ beta hmm xs

beta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob ->
   f emission -> NonEmpty.T f (Vector prob)
beta hmm =
   NonEmpty.scanr
      (\xi betai -> NC.mul (emission hmm xi) betai <> transition hmm)
      (NC.constant 1 (NC.dim $ initial hmm))


alphaBeta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob ->
   NonEmpty.T f emission ->
   (prob, NonEmpty.T f (Vector prob), NonEmpty.T f (Vector prob))
alphaBeta hmm xs =
   let alphas = alpha hmm xs
       betas = beta hmm $ NonEmpty.tail xs
       recipLikelihood = recip $ NC.sumElements $ NonEmpty.last alphas
   in  (recipLikelihood, alphas, betas)



xiFromAlphaBeta ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> prob ->
   NonEmpty.T [] emission ->
   NonEmpty.T [] (Vector prob) ->
   NonEmpty.T [] (Vector prob) ->
   [Matrix prob]
xiFromAlphaBeta hmm recipLikelihood xs alphas betas =
   zipWith3
      (\x alpha0 beta1 ->
         NC.scale recipLikelihood $
         NC.mul
            (NC.outer (NC.mul (emission hmm x) beta1) alpha0)
            (transition hmm))
      (NonEmpty.tail xs)
      (NonEmpty.init alphas)
      (NonEmpty.tail betas)

zetaFromXi ::
   (NC.Product prob) => T distr prob -> [Matrix prob] -> [Vector prob]
zetaFromXi hmm xis =
   map (NC.constant 1 (Matrix.rows $ transition hmm) <>) xis

zetaFromAlphaBeta ::
   (NC.Container Vector prob) =>
   prob ->
   NonEmpty.T [] (Vector prob) ->
   NonEmpty.T [] (Vector prob) ->
   NonEmpty.T [] (Vector prob)
zetaFromAlphaBeta recipLikelihood alphas betas =
   fmap (NC.scale recipLikelihood) $
   NonEmptyC.zipWith NC.mul alphas betas


{- |
In constrast to Math.HiddenMarkovModel.reveal
this does not normalize the vector.
This is slightly simpler but for long sequences
the product of probabilities might be smaller
than the smallest representable number.
-}
reveal ::
   (Distr.EmissionProb distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission,
    Traversable f) =>
   T distr prob -> NonEmpty.T f emission -> NonEmpty.T f State
reveal hmm (NonEmpty.Cons x xs) =
   fmap State $
   uncurry (NonEmpty.scanr Dev.at') $
   mapFst NC.maxIndex $
   mapAccumL
      (\alphai xi ->
         swap $ mapSnd (NC.mul (emission hmm xi)) $
         matrixMaxMul (transition hmm) alphai)
      (NC.mul (emission hmm x) (initial hmm)) xs

matrixMaxMul ::
   (NC.Container Vector a) =>
   Matrix a -> Vector a -> (Vector Int, Vector a)
matrixMaxMul m v =
   mapPair (Vector.fromList, Vector.fromList) $ unzip $
   map ((\x -> (NC.maxIndex x, NC.maxElement x)) . NC.mul v) $
   Matrix.toRows m



{- |
A trained model is a temporary form of a Hidden Markov model
that we need during the training on multiple training sequences.
It allows to collect knowledge over many sequences with 'mergeTrained',
even with mixed supervised and unsupervised training.
You finish the training by converting the trained model
back to a plain modul using 'finishTraining'.

You can create a trained model in three ways:

* supervised training using an emission sequence with associated states,

* unsupervised training using an emission sequence and an existing Hidden Markov Model,

* derive it from state sequence patterns, cf. "Math.HiddenMarkovModel.Pattern".
-}
data Trained distr prob =
   Trained {
      trainedInitial :: Vector prob,
      trainedTransition :: Matrix prob,
      trainedDistribution :: distr
   }
   deriving (Show, Read)

instance
   (NFData distr, NFData prob, Storable prob) =>
      NFData (Trained distr prob) where
   rnf hmm =
      rnf (trainedInitial hmm, trainedTransition hmm, trainedDistribution hmm)


sumTransitions ::
   (NC.Container Vector e) =>
   T distr e -> [Matrix e] -> Matrix e
sumTransitions hmm =
   List.foldl' NC.add (NC.konst 0 $ LinAlg.size $ transition hmm)

{- |
Baum-Welch algorithm
-}
trainUnsupervised ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob, Distr.Emission distr ~ emission) =>
   T distr prob -> NonEmpty.T [] emission -> Trained tdistr prob
trainUnsupervised hmm xs =
   let (recipLikelihood, alphas, betas) = alphaBeta hmm xs
       zetas = zetaFromAlphaBeta recipLikelihood alphas betas

   in  Trained {
          trainedInitial = NonEmpty.head zetas,
          trainedTransition =
             sumTransitions hmm $
             xiFromAlphaBeta hmm recipLikelihood xs alphas betas,
          trainedDistribution =
             Distr.accumulateEmissions $ map (zip (NonEmpty.flatten xs)) $
             List.transpose $ map Vector.toList $ NonEmpty.flatten zetas
       }


mergeTrained ::
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob) =>
   Trained tdistr prob -> Trained tdistr prob -> Trained tdistr prob
mergeTrained hmm0 hmm1 =
   Trained {
      trainedInitial = NC.add (trainedInitial hmm0) (trainedInitial hmm1),
      trainedTransition =
         NC.add (trainedTransition hmm0) (trainedTransition hmm1),
      trainedDistribution =
         Distr.combine
            (trainedDistribution hmm0) (trainedDistribution hmm1)
   }

instance
   (Distr.Estimate tdistr, Distr.Distribution tdistr ~ distr,
    Distr.Probability distr ~ prob) =>
      Sg.Semigroup (Trained tdistr prob) where
   (<>) = mergeTrained


toCells ::
   (Distr.CSV distr, Algo.Field prob, Show prob) =>
   T distr prob -> [[String]]
toCells hmm =
   (HMMCSV.cellsFromVector $ initial hmm) :
   (HMMCSV.cellsFromMatrix $ transition hmm) ++
   [] :
   (Distr.toCells $ distribution hmm)

parseCSV ::
   (Distr.CSV distr, Algo.Field prob, Read prob) =>
   HMMCSV.CSVParser (T distr prob)
parseCSV = do
   v <- HMMCSV.parseNonEmptyVectorCells
   m <- HMMCSV.parseSquareMatrixCells $ Vector.dim v
   HMMCSV.skipEmptyRow
   distr <- Distr.parseCells $ Vector.dim v
   return $ Cons {
      initial = v,
      transition = m,
      distribution = distr
   }
