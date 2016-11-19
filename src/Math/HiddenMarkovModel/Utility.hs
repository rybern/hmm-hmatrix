{-# LANGUAGE FlexibleContexts #-}
module Math.HiddenMarkovModel.Utility where

import qualified Numeric.Container as NC
import Data.Packed.Vector (Vector)

import qualified System.Random as Rnd

import qualified Control.Monad.Trans.State as MS


normalizeProb ::
   (NC.Container Vector a) => Vector a -> Vector a
normalizeProb = snd . normalizeFactor

normalizeFactor ::
   (NC.Container Vector a) => Vector a -> (a, Vector a)
normalizeFactor xs =
   let c = NC.sumElements xs
   in  (c, NC.scale (recip c) xs)

-- see htam:Stochastic
randomItemProp ::
   (Rnd.RandomGen g, Rnd.Random b, Num b, Ord b) =>
   [(a,b)] -> MS.State g a
randomItemProp props =
   let (keys,ps) = unzip props
   in  do p <- MS.state (Rnd.randomR (0, sum ps))
          return $
             fst $ head $ dropWhile ((0<=) . snd) $
             zip keys $ tail $ scanl (-) p ps

attachOnes :: (Num b) => [a] -> [(a,b)]
attachOnes = map (flip (,) 1)
