(ns fda.smoothing.lpk
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [clojure.core.reducers :as r]
            [fda.utils :as utils]))

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
  (let [quantile (utils/normal-quantile (- 1 (/ level 2)))
        hsig2 (m/div (utils/sse y y-estimate)
                     (m/sub (count y) (m/trace smoother-matrix)))
        uysig (->> (m/mmul smoother-matrix (m/transpose smoother-matrix))
                   (m/diagonal)
                   (m/mul hsig2)
                   (m/sqrt))]
    {:lower (m/to-nested-vectors (m/add y-estimate (m/mul quantile uysig)))
     :upper (m/to-nested-vectors (m/sub y-estimate (m/mul quantile uysig)))}))

(defn select-bandwidth
  "Selects the best bandwitdh for a set of curves based on GCV rule."
  [data & {:keys [degree kernel] :or {degree 2 kernel (:gaussian kernels)}}]
  (let [data (if (= (m/dimensionality data) 2) [data] data)
        bandwidths (generate-bandwidths (first (concat data)))
        gcv (fn [data bandwidth]
              (let [[x y] (m/slices data 1)
                    smoother-fn (smoother x degree kernel)
                    smoother-matrix (m/matrix (map #(smoother-fn % bandwidth) x))
                    y-estimate (m/mmul smoother-matrix y)
                    dof (m/trace smoother-matrix)]
                (utils/gcv y y-estimate dof)))
        mean-gcv (fn [bandwidth] (s/mean (map #(gcv % bandwidth) data)))
        gcv-values (pmap mean-gcv bandwidths)]
    (val (apply min-key key (zipmap gcv-values bandwidths)))))

(defn adjust-outliers
  "Reconstruct smooth function and adjust data for outliers by substituting
  outlier elements with estimated values"
  [data degree kernel iterations]
  (let [[x y] (m/slices data 1)
        bandwidth (select-bandwidth data :degree degree :kernel kernel)
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
  (let [[x y] (m/slices data 1)
        quantile (utils/normal-quantile (- 1 (/ level 2)))
        y-adjustments (if iterations (adjust-outliers data degree kernel iterations)
                          [y])
        y-adjusted (last y-adjustments)
        adjusted-data (m/transpose [x y-adjusted])
        bandwidth (or bandwidth (select-bandwidth adjusted-data :degree degree :kernel kernel))
        smoother-fn #((smoother x degree kernel) % bandwidth)
        smoother-matrix (m/array (map smoother-fn x))
        y-estimate (m/mmul smoother-matrix y-adjusted)]
    (with-meta (fn [t] (m/mget (m/mmul (m/array (smoother-fn t)) y-adjusted)))
      {:x x
       :y y
       :y-estimate (m/to-nested-vectors y-estimate)
       :y-adjustments y-adjustments
       :significance (significance y y-estimate smoother-matrix level)})))

