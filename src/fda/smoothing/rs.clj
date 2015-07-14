(ns fda.smoothing.rs
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [clojure.core.matrix.linear :as ml]
            [fda.utils :as utils]))


(defn knot-positions
  "Calculates knot positions for regression spline smoother"
  [x knot-count & {:keys [type] :or {type :uniform}}]
  (condp = type
    :uniform (m/array (rest (drop-last (utils/linspace (m/emin x)
                                                       (m/emax x)
                                                       (+ knot-count 2)))))
    :quantile (m/array (utils/quantile x :probs (m/div (range 1 (inc knot-count))
                                                       (inc knot-count))))))

(defn truncated-power-basis
  [degree & {:keys [knots] :or {knots []}}]
  (fn [x]
    (let [zero-or-more #(if (< % 0) 0  %)
          power-part (m/emap #(m/pow x %) (m/array (range (inc degree))))
          truncated-part (m/emap #(m/pow (zero-or-more (- x %)) degree) knots)]
      (m/join power-part truncated-part))))

(defn- max-knots
  [data]
  (Math/floor (min 50 (/ (m/row-count data) 2))))

(defn- rs-fit
  [data degree knot-count knot-spacing]
  (let [[x y] (m/slices data 1)
        knots (knot-positions x knot-count :type knot-spacing)
        series-fn (truncated-power-basis degree :knots knots)
        series (m/emap series-fn x)
        ls-fit (ml/least-squares (m/matrix series) y)]
    (with-meta (m/mmul series ls-fit)
      {:ls-fit ls-fit
       :fn series-fn
       :series series})))

(defn gcv
  [data degree knot-count knot-spacing]
  (let [y-estimate (rs-fit data degree knot-count knot-spacing)
        y (m/get-column data 1)]
    (utils/gcv y y-estimate (+ knot-count degree 1))))

(defn knot-counts
  "Selects optimal number of knots"
  [data degree knot-spacing]
  (let [[x y] (m/slices data 1)
        data (if (= (m/dimensionality data) 2) [data] data)
        knot-counts (range (inc (m/emax (map max-knots data))))
        gcv #(gcv %1 %2 degree knot-spacing)
        mean-gcv (fn [knot-count] (s/mean (map #(gcv % knot-count) data)))]
    {:knot-count knot-counts
     :gcv (map mean-gcv knot-counts)}))

(defn select-knot-count
  [knot-counts]
  (:knot-count (reduce (fn [b1 b2] (if (< (:gcv b1) (:gcv b2)) b1 b2))
                       knot-counts)))

(defn significance
  [y y-estimate series least-squares-fit dof level]
  (let [score (utils/z-score level)
        ser2 (utils/ser2 y y-estimate dof)
        se (-> series
               (m/mmul least-squares-fit)
               (m/diagonal)
               (m/mul ser2)
               (m/sqrt)
               (m/to-nested-vectors))]
    {:level level
     :se se
     :sig-lower (m/sub y-estimate (m/mul score se))
     :sig-upper (m/add y-estimate (m/mul score se))}))

(defn fit
  "Regression spline smoother"
  [data & {:keys [degree knot-spacing knot-count iterations level]
           :or {degree 2
                knot-spacing :uniform
                level 0.05}}]
  (let [[x y] (m/slices data 1)
        knot-counts (if knot-count [{:knot-count knot-count
                                     :gcv (gcv data knot-count degree knot-spacing)}]
                        (knot-counts data degree knot-spacing))
        best-knot-count (select-knot-count knot-counts)
        y-estimate (rs-fit data degree best-knot-count knot-spacing)
        {ls-fit :ls-fit series-fn :fn series :series} (meta y-estimate)
        dof (+ 1 best-knot-count degree)]
    (with-meta (fn [t] (m/mget (m/mmul (series-fn t) ls-fit)))
      {:x x
       :y y
       :degree degree
       :knot-spacing knot-spacing
       :knot-counts knot-count
       :selected-knot-count best-knot-count
       :y-estimate (m/to-nested-vectors y-estimate)
       :y-adjustments 0
       :significance (significance y y-estimate series ls-fit dof level)})))
