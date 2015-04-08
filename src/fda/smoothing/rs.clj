(ns fda.smoothing.rs
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [clojure.core.matrix.linear :as ml]
            [fda.utils :refer :all]))


(defn knot-positions
  "Calculates knot positions for regression spline smoother"
  [x number-of-knots & {:keys [type] :or {type :uniform}}]
  (condp = type
    :uniform (m/array (rest (drop-last (linspace (m/emin x)
                                                 (m/emax x)
                                                 (+ number-of-knots 2)))))
    :quantile (m/array (quantile x :probs (m/div (range 1 (inc number-of-knots))
                                                 (inc number-of-knots))))))

(defn truncated-power-basis
  [degree & {:keys [knots] :or {knots []}}]
  (fn [x]
    (let [zero-or-more #(if (< % 0) 0  %)
          power-part (m/emap #(m/pow x %) (m/array (range (inc degree))))
          truncated-part (m/emap #(m/pow (zero-or-more (- x %)) degree) knots)]
      (m/join power-part truncated-part))))

(defn knot-numbers
  [data]
  (let [max-number (min 50 (/ (m/row-count data) 2))]
    (range (inc (Math/floor max-number)))))

(declare fit)
(defn select-knot-number
  "Selects optimal number of knots"
  [data degree knot-spacing]
  (let [[x y] (m/slices data 1)
        knot-numbers (knot-numbers data)
        y-estimate (fn [n] (fit data :degree degree :number-of-knots n
                                     :knot-spacing knot-spacing))
        gcv (pmap #(gcv y (y-estimate %) (+ % degree 1)) knot-numbers)]
    (val (apply min-key key (zipmap gcv knot-numbers)))))

(defn significance
  [data & {:keys [degree knot-spacing number-of-knots level]
           :or {degree 2
                knot-spacing :uniform
                level 0.05}}]
  (let [[x y] (m/slices data 1)
        z-score (normal-quantile (- 1 (/ level 2)))
        number-of-knots (or number-of-knots (select-knot-number data degree knot-spacing))
        knots (knot-positions x number-of-knots :type knot-spacing)
        series-fn (truncated-power-basis degree :knots knots)
        series (m/array (map series-fn x))
        least-squares-fit (least-squares-estimator series)
        y-estimate (m/to-nested-vectors (m/mmul series least-squares-fit y))
        ser2 (ser2 y y-estimate (+ 1 number-of-knots degree))
        se (-> series
               (m/mmul least-squares-fit)
               (m/diagonal)
               (m/mul ser2)
               (m/sqrt)
               (m/to-nested-vectors))]
    {:fit y-estimate
     :se se
     :sig-lower (m/sub y-estimate (m/mul z-score se))
     :sig-upper (m/add y-estimate (m/mul z-score se))}))

(defn fit
  "Regression spline smoother"
  [data & {:keys [t degree knot-spacing number-of-knots]
           :or {degree 2
                knot-spacing :uniform}}]
  (let [[x y] (m/slices data 1)
        number-of-knots (or number-of-knots (select-knot-number data degree knot-spacing))
        knots (knot-positions x number-of-knots :type knot-spacing)
        series-fn (truncated-power-basis degree :knots knots)
        series (m/emap series-fn x)
        t-series (if t (m/emap series-fn t) series)
        least-squares-fit (ml/least-squares (m/matrix series) y)]
    (m/mmul t-series least-squares-fit)))
