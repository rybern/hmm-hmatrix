{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
module Math.HiddenMarkovModel.Distribution (
   State(..),
   Emission, Probability, Trained,
   Info(..), Generate(..), EmissionProb(..), Estimate(..),

   Discrete(..), DiscreteTrained(..),
   Gaussian(..), GaussianTrained(..), gaussian,

   CSV(..), HMMCSV.CSVParser, CSVSymbol(..),
   ) where

import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Utility (randomItemProp, normalizeProb)

import qualified Numeric.LinearAlgebra.HMatrix as HMatrix
import qualified Numeric.LinearAlgebra.Algorithms as Algo
import qualified Numeric.Container as NC
import qualified Data.Packed.Matrix as Matrix
import qualified Data.Packed.Vector as Vector
import Numeric.Container ((<>))
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)
import Foreign.Storable (Storable)

import qualified System.Random as Rnd

import qualified Text.CSV.Lazy.String as CSV
import Text.Read.HT (maybeRead)
import Text.Printf (printf)

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.State as MS
import Control.DeepSeq (NFData, rnf)
import Control.Monad (liftM2)

import qualified Data.NonEmpty as NonEmpty
import qualified Data.Foldable as Fold
import qualified Data.Map as Map
import qualified Data.Set as Set
import qualified Data.Array as Array
import qualified Data.List as List
import Data.Foldable (foldMap)
import Data.Tuple.HT (mapFst)
import Data.Array (Array, Ix, listArray, (!))
import Data.Map (Map)
import Data.Maybe (listToMaybe)


newtype State = State Int
   deriving (Eq, Ord, Show, Read, Ix)

instance Enum State where
   toEnum = State
   fromEnum (State n) = n

instance NFData State where
   rnf (State n) = rnf n


type family Probability distr
type family Emission distr
type family Trained distr


class
   (NC.Container Vector (Probability distr), NC.Product (Probability distr)) =>
      Info distr where
   numberOfStates :: distr -> Int

class
   (NC.Container Vector (Probability distr), NC.Product (Probability distr)) =>
      Generate distr where
   generate ::
      (Rnd.RandomGen g, Probability distr ~ prob, Emission distr ~ emission) =>
      distr -> State -> MS.State g emission

class
   (NC.Container Vector (Probability distr), NC.Product (Probability distr)) =>
      EmissionProb distr where
   {-
   This function could be implemented generically in terms of emissionStateProb
   but that would require an Info constraint.
   -}
   emissionProb :: distr -> Emission distr -> Vector (Probability distr)
   emissionStateProb :: distr -> Emission distr -> State -> Probability distr
   emissionStateProb distr e (State s) = NC.atIndex (emissionProb distr e) s

class
   (EmissionProb (Distribution tdistr),
    Trained (Distribution tdistr) ~ tdistr) =>
      Estimate tdistr where
   type Distribution tdistr
   accumulateEmissions ::
      (Distribution tdistr ~ distr, Probability distr ~ prob) =>
      [[(Emission distr, prob)]] -> tdistr
   -- could as well be in Semigroup class
   combine :: tdistr -> tdistr -> tdistr
   normalize :: (Distribution tdistr ~ distr) => tdistr -> distr



newtype Discrete prob symbol = Discrete (Map symbol (Vector prob))
   deriving (Show)

newtype DiscreteTrained prob symbol = DiscreteTrained (Map symbol (Vector prob))
   deriving (Show)

type instance Probability (Discrete prob symbol) = prob
type instance Emission (Discrete prob symbol) = symbol

type instance Trained (Discrete prob symbol) = DiscreteTrained prob symbol


instance (NFData prob, NFData symbol) => NFData (Discrete prob symbol) where
   rnf (Discrete m) = rnf m

instance
   (NFData prob, NFData symbol) =>
      NFData (DiscreteTrained prob symbol) where
   rnf (DiscreteTrained m) = rnf m

instance
   (NC.Container Vector prob, NC.Product prob, Ord symbol) =>
      Info (Discrete prob symbol) where
   numberOfStates (Discrete m) = Vector.dim $ snd $ Map.findMin m

instance
   (NC.Container Vector prob, NC.Product prob, Ord symbol,
    Ord prob, Rnd.Random prob) =>
      Generate (Discrete prob symbol) where
   generate (Discrete m) (State state) =
      randomItemProp $ Map.toAscList $ fmap (flip NC.atIndex state) m

instance
   (NC.Container Vector prob, NC.Product prob, Ord symbol) =>
      EmissionProb (Discrete prob symbol) where
   emissionProb (Discrete m) =
      mapLookup "emitDiscrete: unknown emission symbol" m

instance
   (NC.Container Vector prob, NC.Product prob, Ord symbol) =>
      Estimate (DiscreteTrained prob symbol) where
   type Distribution (DiscreteTrained prob symbol) = Discrete prob symbol
   accumulateEmissions grouped =
      let set = Set.toAscList $ foldMap (Set.fromList . map fst) grouped
          emi = Map.fromAscList $ zip set [0..]
      in  DiscreteTrained $ Map.fromAscList $ zip set $
          transposeVectorList $
          map
             (NC.accum (NC.konst 0 (length set)) (+) .
              map (mapFst
                 (mapLookup "estimateDiscrete: unknown emission symbol" emi)))
             grouped
   combine (DiscreteTrained distr0) (DiscreteTrained distr1) =
      DiscreteTrained $ Map.unionWith NC.add distr0 distr1
   normalize (DiscreteTrained distr) =
      Discrete $ Map.fromAscList $ zip (Map.keys distr) $
      transposeVectorList $ map normalizeProb $
      transposeVectorList $ Map.elems distr

transposeVectorList :: (Matrix.Element a) => [Vector a] -> [Vector a]
transposeVectorList = Matrix.toRows . Matrix.fromColumns

mapLookup :: (Ord k) => String -> Map.Map k a -> k -> a
mapLookup msg dict x =
   Map.findWithDefault (error msg) x dict


newtype Gaussian a = Gaussian (Array State (Vector a, Matrix a, a))
   deriving (Show)

newtype GaussianTrained a = GaussianTrained (Map State (Vector a, Matrix a, a))
   deriving (Show)

type instance Probability (Gaussian a) = a
type instance Emission (Gaussian a) = Vector a

type instance Trained (Gaussian a) = GaussianTrained a


instance (NFData a, Storable a) => NFData (Gaussian a) where
   rnf (Gaussian params) = rnf params

instance (NFData a, Storable a) => NFData (GaussianTrained a) where
   rnf (GaussianTrained params) = rnf params

instance (Algo.Field a) => Info (Gaussian a) where
   numberOfStates (Gaussian params) = Array.rangeSize $ Array.bounds params

instance (Algo.Field a) => Generate (Gaussian a) where
   generate (Gaussian allParams) state = do
      let (center, covarianceChol, _c) = allParams ! state
      seed <- MS.state Rnd.random
      return $
         NC.add center $
         NC.cmap realToFrac
               (NC.randomVector seed NC.Gaussian (Vector.dim center))
            <> covarianceChol

instance (HMatrix.Numeric a, Algo.Field a) => EmissionProb (Gaussian a) where
   emissionProb (Gaussian allParams) x =
      Vector.fromList $ map (emissionProbGen x) $ Array.elems allParams
   emissionStateProb (Gaussian allParams) x s =
      emissionProbGen x $ allParams ! s

emissionProbGen ::
   (HMatrix.Numeric a, Algo.Field a) =>
   Vector a -> (Vector a, Matrix a, a) -> a
emissionProbGen =
   let cholSolve m x = Matrix.flatten $ Algo.cholSolve m $ Matrix.asColumn x
   in  \x (center, covarianceChol, c) ->
         let x0 = NC.sub x center
         in  c * exp ((-1/2) * NC.dot x0 (cholSolve covarianceChol x0))


instance (HMatrix.Numeric a, Algo.Field a) => Estimate (GaussianTrained a) where
   type Distribution (GaussianTrained a) = Gaussian a
   accumulateEmissions =
      let params xs =
             let center =
                    NonEmpty.foldl1Map NC.add (uncurry $ flip NC.scale) xs
                 covariance =
                    NonEmpty.foldl1Map NC.add (\(x,c) -> NC.scale c $ NC.outer x x) xs
             in  (center, covariance, Fold.sum $ fmap snd xs)
      in  GaussianTrained . fmap params . Map.mapMaybe NonEmpty.fetch .
          Map.fromList . zip [State 0 ..]
   combine (GaussianTrained distr0) (GaussianTrained distr1) =
      let comb (center0, covariance0, weight0)
               (center1, covariance1, weight1) =
             (NC.add center0 center1,
              NC.add covariance0 covariance1,
              weight0 + weight1)
      in  GaussianTrained $ Map.unionWith comb distr0 distr1
   {-
     Sum_i (xi-mi) * (xi-mi)^T
   = Sum_i xi*xi^T + Sum_i mi*mi^T - Sum_i xi*mi^T - Sum_i mi*xi^T
   = Sum_i xi*xi^T - Sum_i mi*mi^T
   = Sum_i xi*xi^T - n * mi*mi^T
   -}
   normalize (GaussianTrained distr) =
      let params (centerSum, covarianceSum, weight) =
             let c = recip weight
                 center = NC.scale c centerSum
             in  (center,
                  NC.sub (NC.scale c covarianceSum) (NC.outer center center))
      in  Gaussian $
          Array.array (fst $ Map.findMin distr, fst $ Map.findMax distr) $
          Map.toList $ fmap (gaussianParameters . params) distr

gaussian ::
   (Algo.Field prob) =>
   [(Vector prob, Matrix prob)] -> Gaussian prob
gaussian =
   consGaussian . map gaussianParameters

gaussianParameters ::
   (Algo.Field prob) =>
   (Vector prob, Matrix prob) -> (Vector prob, Matrix prob, prob)
gaussianParameters (center, covariance) =
   gaussianFromCholesky center $ Algo.chol covariance

consGaussian :: [(Vector a, Matrix a, a)] -> Gaussian a
consGaussian xs =
   Gaussian $ listArray (State 0, State $ length xs - 1) xs

gaussianFromCholesky ::
   (Algo.Field prob) =>
   Vector prob -> Matrix prob -> (Vector prob, Matrix prob, prob)
gaussianFromCholesky center covarianceChol =
   let covarianceSqrtDet = NC.prodElements $ Matrix.takeDiag covarianceChol
   in  (center, covarianceChol,
        recip (sqrt (2*pi) ^ Vector.dim center * covarianceSqrtDet))



class CSV distr where
   toCells :: distr -> [[String]]
   parseCells :: Int -> HMMCSV.CSVParser distr

class (Ord symbol) => CSVSymbol symbol where
   cellFromSymbol :: symbol -> String
   symbolFromCell :: String -> Maybe symbol

instance CSVSymbol Char where
   cellFromSymbol = (:[])
   symbolFromCell = listToMaybe

instance CSVSymbol Int where
   cellFromSymbol = show
   symbolFromCell = maybeRead


instance
   (Algo.Field prob, Show prob, Read prob, CSVSymbol symbol) =>
      CSV (Discrete prob symbol) where
   toCells (Discrete m) =
      map
         (\(symbol, probs) ->
            cellFromSymbol symbol : HMMCSV.cellsFromVector probs) $
      Map.toAscList m
   parseCells n =
      fmap (Discrete . Map.fromList) $
      HMMCSV.manyRowsUntilEnd $ parseSymbolProb n

parseSymbolProb ::
   (Algo.Field prob, Read prob, CSVSymbol symbol) =>
   Int -> CSV.CSVRow -> HMMCSV.CSVParser (symbol, Vector prob)
parseSymbolProb n row =
   case row of
      [] -> MT.lift $ ME.throw "missing symbol"
      c:cs ->
         liftM2 (,)
            (let str = CSV.csvFieldContent c
             in  MT.lift $ ME.fromMaybe (printf "unknown symbol %s" str) $
                 symbolFromCell str)
            (do v <- HMMCSV.parseVectorFields cs
                HMMCSV.assert (n == Vector.dim v)
                   (printf "number of states (%d) and size of probability vector (%d) mismatch"
                      n (Vector.dim v))
                return v)


instance (Algo.Field a, Eq a, Show a, Read a) => CSV (Gaussian a) where
   toCells (Gaussian params) =
      List.intercalate [[]] $
      map
         (\(center, covarianceChol, _) ->
            HMMCSV.cellsFromVector center :
            HMMCSV.cellsFromMatrix covarianceChol) $
      Array.elems params
   parseCells n = do
      gs <- HMMCSV.manySepUntilEnd parseSingleGaussian
      HMMCSV.assert (length gs == n) $
         printf "number of states (%d) and number of Gaussians (%d) mismatch"
            n (length gs)
      return $ consGaussian gs

parseSingleGaussian ::
   (Algo.Field prob, Eq prob, Read prob) =>
   HMMCSV.CSVParser (Vector prob, Matrix prob, prob)
parseSingleGaussian = do
   center <- HMMCSV.parseNonEmptyVectorCells
   covarianceChol <- HMMCSV.parseSquareMatrixCells $ Vector.dim center
   HMMCSV.assert (isUpperTriang covarianceChol) $
      "matrices must be upper triangular"
   return $ gaussianFromCholesky center covarianceChol


{-
Maybe this test is too strict.
It would also be ok, and certainly more intuitive
to use an orthogonal but not normalized matrix.
We could get such a matrix from the eigensystem.
-}
isUpperTriang :: (Algo.Field a, Eq a) => Matrix a -> Bool
isUpperTriang m =
   let cleared = Matrix.mapMatrixWithIndex (\(i,j) x -> if i>j then x else 0) m
   in  NC.minElement cleared == 0 &&
       NC.maxElement cleared == 0
