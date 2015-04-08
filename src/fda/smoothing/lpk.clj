(ns fda.smoothing.lpk
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.stats :as s]
            [taoensso.timbre :as timbre]
            [clojure.core.reducers :as r]
            [fda.utils :refer :all]))

(timbre/refer-timbre)

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
  [obs-points bandwidth degree kernel]
  (fn [t]
    (let [diff (scaled-diff obs-points t bandwidth)
          w (m/diagonal-matrix (kernel diff))
          x (taylor-series diff degree)
          col (least-squares-estimator x w)]
      (m/get-row col 0))))

(defn smoother-matrix
  [obs-points bandwidth degree kernel]
  (map (smoother obs-points bandwidth degree kernel) obs-points))

(defn generate-bandwidths
  "Generates possible bandwidths for evaluation with GCV rule."
  [data observation-count]
  (let [steps 15
        x (m/get-column data 0)
        xrange (- (m/emax x) (m/emin x))
        min-bandwidth (* 1.1 (/ xrange observation-count))
        max-bandwidth (/ xrange 8)]
        (logspace (m/log10 min-bandwidth) (m/log10 max-bandwidth) steps)))

(declare fit)
(defn select-bandwidth
  "Selects the best bandwitdh for one curve based on GCV rule."
  [data & {:keys [degree kernel]
           :or {degree 2
                kernel (:gaussian kernels)}}]
  (let [[x y] (m/slices data 1)
        bandwidths (m/to-nested-vectors (generate-bandwidths data (count data)))
        gcv (fn [bandwidth]
              (let [smoother (smoother x bandwidth degree kernel)
                    smoother-matrix (m/matrix (map smoother x))
                    y-estimate (m/mmul smoother-matrix y)
                    dof (m/trace smoother-matrix)]
                (gcv y y-estimate dof)))
        gcv-values (pmap gcv bandwidths)]
    (val (apply min-key key (zipmap gcv-values bandwidths)))))

(defn select-group-bandwidth
  "Selects the best bandwitdh for a set of curves based on GCV rule."
  [data & {:keys [degree kernel]
             :or {degree 2
                  kernel (:gaussian kernels)}}]
  (let [[x y] (m/slices data 1)
        observation-count (s/mean (map #(count %) data))
        bandwidths (generate-bandwidths (concat data) observation-count)
        gcv (fn [data bandwidth]
              (let [smoother (smoother x bandwidth degree kernel)
                    smoother-matrix (m/matrix (map smoother x))
                    y-estimate (m/mmul smoother-matrix y)
                    dof (m/trace smoother-matrix)]
                (gcv y y-estimate dof)))
        mean-gcv (fn [bandwidth] (s/mean (map #(gcv % bandwidth) data)))
        gcv-values (pmap mean-gcv bandwidths)]
    (val (apply min-key key (zipmap gcv-values bandwidths)))))

(defn significance
  "Calculates significance level at design points."
  [data & {:keys [degree kernel bandwidth level]
             :or {degree 2
                  kernel (:gaussian kernels)
                  level 0.05}}]
  (let [[x y] (m/slices data 1)
        bandwidth (or bandwidth (select-bandwidth data :degree degree :kernel kernel))
        quantile (normal-quantile (- 1 (/ level 2)))
        smoother (smoother x bandwidth degree kernel)
        smoother-matrix (m/matrix (map smoother x))
        y-estimate (m/mmul smoother-matrix y)
        hsig2 (m/div (s/sum (m/pow (m/sub y y-estimate) 2))
                     (m/sub (count data) (m/trace smoother-matrix)))
        uysig (->> (m/mmul smoother-matrix (m/transpose smoother-matrix))
                   (m/diagonal)
                   (m/mul hsig2)
                   (m/sqrt))]
    {:position x
     :lower (m/add y-estimate (m/mul quantile uysig))
     :upper (m/sub y-estimate (m/mul quantile uysig))}))

(defn fit
  "Reconstruct smooth function using Local Polynomial Kernel method"
  [data & {:keys [t degree kernel bandwidth]
             :or {degree 2
                  kernel (:gaussian kernels)}}]
  (let [[x y] (m/slices data 1)
        bandwidth (or bandwidth (select-bandwidth data :degree degree :kernel kernel))
        t (or t x)]
    (m/to-nested-vectors
     (m/mmul (m/array (map (smoother x bandwidth degree kernel) t)) y))))

(defn robust-fit
  "Reconstruct smooth function and adjust data for outliers by substituting
  outlier elements with estimated values"
  [data & {:keys [t degree kernel bandwidth iterations]
             :or {degree 2
                  kernel (:gaussian kernels)
                  iterations 3}
           :as args}]
  (let [[x y] (m/slices data 1)
        fit #(fit %1 :t %2 :degree degree :kernel kernel :bandwidth bandwidth)]
    (loop [y-new y iteration iterations]
      (let [data-new (m/transpose [x y-new])
            y-estimate (fit data-new x)
            residuals (m/sub y y-estimate)
            sd (s/sd residuals)
            y-robust (m/emap #(if (< (m/abs %1) (* 2 sd)) %2 %3)
                             residuals y y-estimate)]
        (if (<= iteration 0)
          {:data (m/transpose [x y-robust])
           :fit-at-x y-estimate
           :fit (if t (fit data-new t) y-estimate)}
          (recur y-robust (dec iteration)))))))
