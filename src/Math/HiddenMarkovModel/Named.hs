module Math.HiddenMarkovModel.Named (
   T(..),
   Discrete,
   Gaussian,
   fromModelAndNames,
   toCSV,
   fromCSV,
   ) where

import qualified Math.HiddenMarkovModel.Distribution as Distr
import qualified Math.HiddenMarkovModel.Private as HMM
import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Distribution (State(..))
import Math.HiddenMarkovModel.Utility (attachOnes)

import qualified Numeric.LinearAlgebra.Algorithms as Algo
import qualified Data.Packed.Vector as Vector

import qualified Text.CSV.Lazy.String as CSV
import Text.Printf (printf)

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.State as MS
import Control.DeepSeq (NFData, rnf)
import Foreign.Storable (Storable)

import qualified Data.Map as Map
import qualified Data.List as List
import Data.Tuple.HT (swap)
import Data.Map (Map)


{- |
A Hidden Markov Model with names for each state.

Although 'nameFromStateMap' and 'stateFromNameMap' are exported
you must be careful to keep them consistent when you alter them.
-}
data T distr prob =
   Cons {
      model :: HMM.T distr prob,
      nameFromStateMap :: Map State String,
      stateFromNameMap :: Map String State
   }
   deriving (Show, Read)

type Discrete prob symbol = T (Distr.Discrete prob symbol) prob
type Gaussian a = T (Distr.Gaussian a) a


instance
   (NFData distr, NFData prob, Storable prob) =>
      NFData (T distr prob) where
   rnf hmm = rnf (model hmm, nameFromStateMap hmm, stateFromNameMap hmm)


fromModelAndNames :: HMM.T distr prob -> [String] -> T distr prob
fromModelAndNames md names =
   let m = Map.fromList $ zip [State 0 ..] names
   in  Cons {
          model = md,
          nameFromStateMap = m,
          stateFromNameMap = inverseMap m
       }

inverseMap :: Map State String -> Map String State
inverseMap =
   Map.fromListWith (error "duplicate label") .
   map swap . Map.toList


toCSV ::
   (Distr.CSV distr, Algo.Field prob, Show prob) =>
   T distr prob -> String
toCSV hmm =
   CSV.ppCSVTable $ snd $ CSV.toCSVTable $ HMMCSV.padTable "" $
      Map.elems (nameFromStateMap hmm) : HMM.toCells (model hmm)

fromCSV ::
   (Distr.CSV distr, Algo.Field prob, Read prob) =>
   String -> ME.Exceptional String (T distr prob)
fromCSV =
   MS.evalStateT parseCSV . map HMMCSV.fixShortRow . CSV.parseCSV

parseCSV ::
   (Distr.CSV distr, Algo.Field prob, Read prob) =>
   HMMCSV.CSVParser (T distr prob)
parseCSV = do
   names <- HMMCSV.parseStringList =<< HMMCSV.getRow
   let duplicateNames =
         Map.keys $ Map.filter (> (1::Int)) $
         Map.fromListWith (+) $ attachOnes names
    in HMMCSV.assert (null duplicateNames) $
          "duplicate names: " ++ List.intercalate ", " duplicateNames
   md <- HMM.parseCSV
   let n = length names
       m = Vector.dim (HMM.initial md)
    in HMMCSV.assert (n == m) $
          printf "got %d state names for %d state" n m
   return $ fromModelAndNames md names
