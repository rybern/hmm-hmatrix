{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{- |
This module provides a simple way to train
the transition matrix and initial probability vector
using simple patterns of state sequences.

You may create a trained model using semigroup combinators like this:

> let a = atom $ HMM.state 0
>     b = atom $ HMM.state 1
>     distr =
>        Distr.DiscreteTrained $ Map.fromList $
>        ('a', Vector.fromList [1,2]) :
>        ('b', Vector.fromList [4,3]) :
>        ('c', Vector.fromList [0,1]) :
>        []
> in  finish 2 distr $ replicate 5 $ replicate 10 a <> replicate 20 b
-}
module Math.HiddenMarkovModel.Pattern (
   T,
   atom,
   append,
   replicate,
   finish,
   ) where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel as HMM
import Math.HiddenMarkovModel.Private (Trained(..))
import Math.HiddenMarkovModel.Distribution (State(State))

import qualified Numeric.LinearAlgebra.Algorithms as Algo
import qualified Numeric.Container as NC
import qualified Data.Packed.Vector as Vector
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import qualified Data.Map as Map
import Data.Semigroup (Semigroup, (<>), stimes)

import Prelude hiding (replicate)


newtype T prob = Cons (Int -> (State, Matrix prob, State))

atom ::
   (NC.Container Vector prob) =>
   State -> T prob
atom s = Cons $ \n -> (s, NC.konst 0 (n,n), s)


instance (Algo.Field prob) => Semigroup (T prob) where
   (<>) = append
   stimes k = replicate $ fromIntegral k


infixl 5 `append`

append ::
   (NC.Container Vector prob) =>
   T prob -> T prob -> T prob
append (Cons f) (Cons g) =
   Cons $ \n ->
      case (f n, g n) of
         ((sai, ma, sao), (sbi, mb, sbo)) ->
            (sai, increment (sbi,sao) 1 $ NC.add ma mb, sbo)

replicate ::
   (NC.Container Vector prob) =>
   Int -> T prob -> T prob
replicate ki (Cons f) =
   Cons $ \n ->
      case f n of
         (si, m, so) ->
            let k = fromIntegral ki
            in  (si, increment (si,so) (k-1) $ NC.scale k m, so)

increment ::
   (NC.Container Vector a) =>
   (State, State) -> a -> Matrix a -> Matrix a
increment (State i, State j) x m  =  NC.accum m (+) [((i,j), x)]


finish ::
   (NC.Container Vector prob) =>
   Int -> tdistr -> T prob -> Trained tdistr prob
finish n tdistr (Cons f) =
   case f n of
      (State si, m, _so) ->
         Trained {
            trainedInitial = NC.assoc n 0 [(si,1)],
            trainedTransition = m,
            trainedDistribution = tdistr
         }


_example :: HMM.DiscreteTrained Double Char
_example =
   let a = atom $ HMM.state 0
       b = atom $ HMM.state 1
       distr =
          Distr.DiscreteTrained $ Map.fromList $
          ('a', Vector.fromList [1,2]) :
          ('b', Vector.fromList [4,3]) :
          ('c', Vector.fromList [0,1]) :
          []
   in  finish 2 distr $ replicate 5 $ replicate 10 a <> replicate 20 b
