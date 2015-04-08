(ns fda.utils-test
  (:use midje.sweet
        fda.utils)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [fda.smoothing.rs-test :refer [non-progesterone5 example24-fit]]))

(facts "about lagged differences"
  (diff [1 3 4 5]) => [2 1 1]
  (diff [1 4 5 7 8] :order 2) => [4 3 3])

(facts "about linear spacing"
  (fact "linspace function creates equaly spaced points within interval"
    (s/sum (diff (diff (linspace 0.0 65 999)))) => (roughly 0 1E-4))
  (m/to-nested-vectors (linspace 0.0 1.0 3)) => [0.0 0.5 1.0])

(facts "about logarithmic spacing"
  (m/to-nested-vectors (logspace (m/log10 1.0) (m/log10 100.0) 3)) => [1.0 10.0 100.0])

(facts "about factorial"
  (fact "it is a multiple of all numbers up to the specified number"
    (factorial 4) => (* 1 2 3 4)
    (factorial 1) => 1
    (factorial 10) => (* 1 2 3 4 5 6 7 8 9 10))
  (fact "factorial of 0 is 1"
    (factorial 0) => 1))

(facts "about least squares"
  (let [solution  (least-squares (m/array [[2 0] [-1 1] [0 2]]) (m/array [1 0 -1]))]
    (m/to-nested-vectors solution) => (just [(roughly 1/3) (roughly -1/3)])))

(facts "about quantile"
  (fact "it computes matlab style quantiles"
    (quantile [0 1 2 3 4 5 6 7 8 9 10] :probs [0.2]) => (just [(roughly 1.7)])
    (quantile [0 1 2 3 4 5 6 7 8 9 10] :probs [0.1]) => (just [(roughly 0.6)])))

(facts "about normal distribution quantile"
  (normal-quantile 0.975) => (roughly 1.96)
  (normal-quantile 0.95) => (roughly 1.645)
  (normal-quantile 0.995) => (roughly 2.575))

(facts "about sum of squared errors"
  (sse [1 2 3] [2 2 2]) => 2.0)

(facts "about squared standard error of regression"
  (ser2 (m/get-column non-progesterone5 1) example24-fit 8) => (roughly 0.015 1E-4))
 
(facts "about gcv rule"
  (gcv (m/get-column non-progesterone5 1) example24-fit 8) => (roughly 0.52976))
