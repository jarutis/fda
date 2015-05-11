(ns fda.smoothing.lpk-test
  (:use midje.sweet
        fda.smoothing.lpk)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [fda.utils :as utils]))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Reproduced examples from Jin-Ting Zhang - Analysis of Variance for Functional
;; Data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example 2.1

(def non-progesterone1
  [[-8  -0.10651824]
   [-7   0.04606832]
   [-6  -0.02592671]
   [-5  -0.09207704]
   [-4  -0.05878623]
   [-3  -0.14683805]
   [-2   0.31639465]
   [-1   0.38676637]
   [ 0   0.21613959]
   [ 1   0.49649667]
   [ 2   0.74412508]
   [ 3   1.14207080]
   [ 4   1.46769922]
   [ 5   2.10771918]
   [ 6   2.43389387]
   [ 7   3.00657616]
   [ 8   2.97612922]
   [ 9   3.03166915]
   [ 10  2.75816422]
   [ 11  2.91728551]
   [ 12  2.77770041]
   [ 13  3.20467519]
   [ 14  2.86933860]
   [ 15  1.84179351]])

(def example21-fit032
  [-0.106518 0.046068 -0.025926 -0.092077 -0.058786 -0.146838 0.316394 0.386766
   0.216139 0.496496 0.744125 1.142070 1.467699 2.107719 2.433893 3.006576
   2.976129 3.031669 2.758164 2.917285 2.777700 3.204675 2.869338 1.841793])

(def example21-fit522
  [-0.006203 -0.062398 -0.087043 -0.076300 -0.027154 0.062281 0.192455 0.362126
   0.568208 0.805745 1.068008 1.346694 1.632203 1.913983 2.180903 2.421677
   2.625289 2.781412 2.880729 2.915106 2.877529 2.761763 2.561736 2.270676])

(def example21-fit130
  [-0.098250 -0.005220 -0.027929 -0.080182 -0.096702 0.002335 0.197738 0.306095
   0.334956 0.474064 0.753487 1.111864 1.544724 2.039343 2.519843 2.878568
   3.016537 2.967017 2.871099 2.856362 2.954654 3.035915 2.773713 1.854062])

(facts "about scaled difference"
  (m/to-nested-vectors (scaled-diff [3 4 5][1 2 3] 2)) =>
  (just (map roughly [1 1 1])))

(facts "about local taylor expansion"
  (m/to-nested-vectors (taylor-series [1 2] 2)) =>
  (just (map just (m/emap roughly [[1 1 1][1 2 4]]))))

(facts "about bandwidth selection for LPK smoother"
  (select-bandwidth non-progesterone1 :degree 2) => (roughly 1.307))

(facts "about reconstructing smooth curves"
  (fact "smooth data with 0.32 bandwidth is the same as in figure 2.1"
    (:y-estimate (meta (fit non-progesterone1 :bandwidth 0.32))) =>
    (just (map roughly example21-fit032)))
  (fact "smooth data with 1.30 bandwidth is the same as in figure 2.1"
    (:y-estimate (meta (fit non-progesterone1 :bandwidth 1.307))) =>
    (just (map roughly example21-fit130)))
  (fact "smooth data with 5.20 bandwidth is the same as in figure 2.1"
    (:y-estimate (meta (fit non-progesterone1 :bandwidth 5.22))) =>
    (just (map roughly example21-fit522))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example 2.2

(def example22-h
  [1.0541666 1.1324857 1.2166235 1.3070123 1.4041165 1.5084350 1.6205039 1.7408989
   1.8702386 2.0091876 2.1584598 2.3188221 2.4910985 2.6761742 2.8750000])

(def example22-gcv
  [0.034610 0.033479 0.032739 0.032487 0.032776 0.033624 0.035018 0.036920 0.039275
   0.042007 0.045031 0.048253 0.051581 0.054930 0.058238])

(fact "candidates for bandwidth selection match"
  (m/to-nested-vectors (generate-bandwidths non-progesterone1)) =>
  (just (map roughly example22-h)))

(fact "GCV matches example 2.2"
  (let [gcv (fn [bandwidth]
              (let [[x y] (m/slices non-progesterone1 1)
                    smoother-fn (smoother x 2 (:gaussian kernels))
                    smoother-matrix (m/matrix (map #(smoother-fn % bandwidth) x))
                    y-estimate (m/mmul smoother-matrix y)
                    dof (m/trace smoother-matrix)]
                (/ (utils/gcv y y-estimate dof) (count non-progesterone1))))]
    (map gcv example22-h)) => (just (map roughly example22-gcv)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example 2.3

(def non-progesterone20
  [[-7  -2.237961]
   [-6  -2.854023]
   [-5  -2.815112]
   [-4  -2.385496]
   [-3  -2.622849]
   [-2  -2.262391]
   [-1  -2.995194]
   [ 0  -3.270483]
   [ 1  -1.719572]
   [ 2  -1.683703]
   [ 3  -1.092581]
   [ 4  -0.118210]
   [ 5   0.707598]
   [ 6  -0.212485]
   [ 7   0.103174]
   [ 8   0.394925]
   [ 9   0.174976]
   [ 10 -3.108084]
   [ 11  0.273360]
   [ 12  0.239464]
   [ 13 -0.093736]
   [ 14 -0.243517]])

(def example23-iterations
  [[-2.417741 -2.389000 -2.370280 -2.350739]
   [-2.571041 -2.581454 -2.590545 -2.602505]
   [-2.635410 -2.623511 -2.623549 -2.628145]
   [-2.682079 -2.637098 -2.616630 -2.598819]
   [-2.729513 -2.680480 -2.652839 -2.624393]
   [-2.746513 -2.733700 -2.723006 -2.710705]
   [-2.671900 -2.713753 -2.732151 -2.750339]
   [-2.446485 -2.524239 -2.559904 -2.594040]
   [-2.046988 -2.119705 -2.149674 -2.175598]
   [-1.507892 -1.543102 -1.552777 -1.559866]
   [-0.918200 -0.910120 -0.897635 -0.884349]
   [-0.391938 -0.356980 -0.331802 -0.303168]
   [-0.028080  0.018617  0.039343  0.064516]
   [ 0.120504  0.196158  0.206114  0.211448]
   [ 0.060581  0.223527  0.243945  0.238733]
   [-0.143807  0.171781  0.237046  0.245049]
   [-0.387266  0.100176  0.225083  0.257256]
   [-0.553563  0.047256  0.210981  0.259311]
   [-0.559382  0.025722  0.181645  0.229088]
   [-0.398829  0.012476  0.112648  0.142887]
   [-0.165924 -0.050074 -0.027342 -0.021115]
   [-0.036614 -0.229442 -0.259595 -0.266158]])

(fact "lpk fit matches example 2.3"
  (:y-estimate (meta (fit non-progesterone20))) =>
  (just (map roughly (m/get-column example23-iterations 0))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example 2.4

(facts "about outlier adjustmens"
  (fact "first iteration matches example 2.4"
    (:y-estimate (meta (fit non-progesterone20 :iterations 1))) =>
    (just (map roughly (m/get-column example23-iterations 1)))))
