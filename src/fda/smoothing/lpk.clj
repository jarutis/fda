(ns fda.smoothing.lpk
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [clojure.core.reducers :as r]
            [fda.utils :as utils])
  (:gen-class :name fda.smoothing.Lpk))

(def kernels {:uniform (fn [x] (m/add (m/mul x 0) 1/2))
              :gaussian (fn [x] (m/div (m/exp (m/mul -0.5 (m/pow x 2)))
                                      (m/pow (* 2 Math/PI) 1/2)))})

(defn scaled-diff
  [x t bandwidth]
  (m/div (m/sub (m/array x) (m/array t)) bandwidth))

(defn taylor-series
  [diff degree]
  (m/transpose (m/array (map #(m/pow diff %) (range (inc degree))))))

(defn smoother
  "Calculates component for Local Polynomial Kernel smoother matrix of a
  function observed at obs-points."
  [obs-points degree kernel]
  (fn [t bandwidth]
    (let [diff (scaled-diff obs-points t bandwidth)
          w (m/diagonal-matrix (kernel diff))
          x (taylor-series diff degree)
          col (utils/least-squares-estimator x w)]
      (m/get-row col 0))))

(defn generate-bandwidths
  "Generates possible bandwidths for evaluation with GCV rule."
  [data]
  (let [steps 15
        x (m/get-column data 0)
        xrange (- (m/emax x) (m/emin x))
        min-bandwidth (* 1.1 (/ xrange (count (set x))))
        max-bandwidth (/ xrange 8)]
    (utils/logspace (m/log10 min-bandwidth) (m/log10 max-bandwidth) steps)))

(defn significance
  "Calculates significance level at design points."
  [y y-estimate smoother-matrix level]
  (let [score (utils/z-score level)
        hsig2 (m/div (utils/sse y y-estimate)
                     (m/sub (count y) (m/trace smoother-matrix)))
        se (-> smoother-matrix
               (m/mmul (m/transpose smoother-matrix))
               (m/diagonal)
               (m/mul hsig2)
               (m/sqrt)
               (m/to-nested-vectors))]
    {:level level
     :se se
     :lower (m/add y-estimate (m/mul score se))
     :upper (m/sub y-estimate (m/mul score se))}))

(defn gcv
  [data bandwidth degree kernel]
  (let [[x y] (m/slices data 1)
        smoother-fn (smoother x degree kernel)
        smoother-matrix (m/matrix (map #(smoother-fn % bandwidth) x))
        y-estimate (m/mmul smoother-matrix y)
        dof (m/trace smoother-matrix)]
    (m/div (utils/gcv y y-estimate dof) (count y))))

(defn bandwidths
  "Selects a set possible bandwidths and calculates corresponding GCV values"
  [data & {:keys [degree kernel] :or {degree 2 kernel (:gaussian kernels)}}]
  (let [data (if (= (m/dimensionality data) 2) [data] data)
        bandwidths (generate-bandwidths (first (concat data)))
        gcv #(gcv %1 %2 degree kernel)
        mean-gcv (fn [bandwidth] (s/mean (map #(gcv % bandwidth) data)))]
    {:bandwidths bandwidths
     :gcv (pmap mean-gcv bandwidths)}))

(defn select-bandwidth
  [bandwidths]
  (:bandwidth (reduce (fn [b1 b2] (if (< (:gcv b1) (:gcv b2)) b1 b2)))))

(defn adjust-outliers
  "Reconstruct smooth function and adjust data for outliers by substituting
  outlier elements with estimated values"
  [data degree kernel iterations]
  (let [[x y] (m/slices data 1)
        bandwidth (:best (select-bandwidth data :degree degree :kernel kernel))
        smoother-fn #((smoother x degree kernel) % bandwidth)
        smoother-matrix (m/array (map smoother-fn x))]
    (loop [y-new [y] iteration iterations]
      (let [y-estimate (m/mmul smoother-matrix (last y-new))
            residuals (m/sub y y-estimate)
            sd (s/sd residuals)
            y-robust (m/emap #(if (< (m/abs %1) (m/mul 2 sd)) %2 %3)
                             residuals y y-estimate)]
        (if (<= iteration 0) y-new
            (recur (conj y-new y-robust) (dec iteration)))))))

(defn fit
  "Reconstruct smooth function using Local Polynomial Kernel method"
  [data & {:keys [degree kernel bandwidth iterations level]
             :or {degree 2
                  kernel (:gaussian kernels)
                  level 0.05}}]
  (let [y-adjustments (adjust-outliers data degree kernel iterations)
        adjusted-data (m/transpose [(m/get-column data 0) (last y-adjustments)])
        [x y] (m/slices adjusted-data 1)
        bandwidths (if bandwidth {:bandwidths bandwidth
                                  :gcv (gcv data bandwidth degree kernel)}
                       (bandwidths adjusted-data :degree degree :kernel kernel))
        best-bandwidth (select-bandwidth bandwidths)
        smoother-fn #((smoother x degree kernel) % best-bandwidth)
        smoother-matrix (m/array (map smoother-fn x))
        y-estimate (m/mmul smoother-matrix y)]
    (with-meta (fn [t] (m/mget (m/mmul (m/array (smoother-fn t)) y)))
      {:x x
       :y (m/get-column data 1)
       :y-estimate (m/to-nested-vectors y-estimate)
       :degree degree
       :kernel kernel
       :selected-bandwidth best-bandwidth
       :bandwidths bandwidths
       :y-adjustments y-adjustments
       :significance (significance y y-estimate smoother-matrix level)})))

