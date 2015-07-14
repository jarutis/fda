(ns fda.utils
  (:import (org.apache.commons.math3.stat.descriptive.rank Percentile
                                                           Percentile$EstimationType)
           (org.apache.commons.math3.special Erf))
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s])
  (:gen-class))

(defn diff
  "Calculates lagged difference of specified order"
  [series & {:keys [order] :or {order 1}}]
  (m/sub (drop order series) (drop-last order series)))

(defn linspace
  "Creates a sequence of equaly spaced numbers from *xmin* to *xmax* with count
  specified in *steps*."
  [xmin xmax steps]
  (m/add xmin (m/mul (- xmax xmin) (m/div (range steps) (dec steps)))))

(defn logspace
  "Creates a sequence of logarithmicaly spaced numbers from *xmin* to *xmax* with
  count specified in *steps*."
  [xmin xmax steps]
  (m/pow 10 (linspace xmin xmax steps)))

(defn least-squares-estimator
  "Least squares estimator."
  ([x]
   (m/mmul (m/inverse (m/mmul (m/transpose x) x)) (m/transpose x)))
  ([x weigths]
   (m/mmul (m/inverse (m/mmul (m/transpose x) weigths x))
           (m/mmul (m/transpose x) weigths))))

(defn least-squares
  "Least squares with optional weights."
  [x y]
  (m/mmul (least-squares-estimator x) y))

(defn sse
  "Sum of Squared Errors"
  [y y-estimate]
  (s/sum (m/pow (m/sub y y-estimate) 2)))

(defn gcv
  "Calculates value of generalized cross-validation (GCV) rule. Small GCV value
  indicates a balance between better goodness of fit and low complexity. Here,
  better goodness of fit often means less bias while smaller model complexity
  usually means less rough or smaller variance."
  [y y-estimate degrees-of-freedom]
  (m/div (sse y y-estimate)
         (m/pow (- 1 (/ degrees-of-freedom (m/ecount y))) 2)))

(defn ser2
  "Calculates squared standard error of regression for leas-squares fit"
  [y y-estimate degrees-of-freedom]
  (m/div (sse y y-estimate)
         (- (count y) degrees-of-freedom)))

(defn factorial
  ([n acc] (if (< n 2) acc (recur (dec n) (* acc n))))
  ([n] (factorial n 1N)))

(defn quantile
  "Wrapper around apache commons math Percetile class"
  [data & {:keys [probs type]
           :or {probs [0.0 0.25 0.5 0.75 1.0]
                type :R_5}}]
  (let [estimation-type (Percentile$EstimationType/valueOf (name type))
        percentile (.withEstimationType (Percentile. ) estimation-type)]
    (map #(.evaluate percentile (double-array data) (* 100 %)) probs)))

(defn normal-quantile
  "Quantile of normal distribution"
  [x]
  (* (m/sqrt 2) (Erf/erfInv (- (* 2 x) 1))))

(defn z-score
  "Z-score based on confidence level"
  [level]
  (normal-quantile (- 1 (/ level 2))))
