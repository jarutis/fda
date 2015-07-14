(defproject fda "0.0.1-SNAPSHOT"
  :description "FDA library"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :plugins [[lein-gorilla "0.3.4"]
            [cider/cider-nrepl "0.9.0-SNAPSHOT"]]
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.apache.commons/commons-math3 "3.5"]
                 [net.mikera/core.matrix "0.33.2"]
                 [net.mikera/core.matrix.stats "0.5.0"]
                 [net.mikera/vectorz-clj "0.29.0"]]

  :profiles {:dev {:dependencies [[midje "1.6.3"]
                                  [com.taoensso/timbre "3.4.0"]
                                  [incanter "1.9.0"]
                                  [net.sourceforge.jmatio/jmatio "1.0"]]
                   :resource-paths ["dev-resources/zhang_fanova/matlabdat"]}})
