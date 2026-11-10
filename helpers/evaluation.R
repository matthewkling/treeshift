# Model-evaluation settings: minimum sample size, metric labels, and the fields
# its results depend on.

# Below this many held-out plots, metrics are not reported at all. A small unit
# can easily contain single-digit plot counts, and an AUC over eight plots is
# worse than no answer because it looks like one.
eval_min_n <- 10

# Display labels for the metrics cv_performance() returns, in the order they
# should be read. Both metric sets are listed; only the applicable ones appear
# for a given response (binary -> n/auc/tss/tss_threshold/brier,
# continuous -> n/rmse/mae/bias/r2).
eval_metric_labels <- c(n             = "n",
                        auc           = "AUC",
                        tss           = "TSS",
                        tss_threshold = "TSS thresh.",
                        brier         = "Brier",
                        rmse          = "RMSE",
                        mae           = "MAE",
                        bias          = "bias",
                        r2            = "R\u00b2")

# Fields the model-evaluation results depend on, and how to name them to the
# user. Note clim_agg is absent: eval_var_set() extracts climate from the raw
# rasters at plot coordinates, never from the aggregated modeling grid, so the
# grain slider does not affect evaluation at all.
eval_stamp_labels <- c(track      = "focal entity type",
                       focal      = "focal entity",
                       measure    = "community measure",
                       fia_var    = "FIA variable",
                       bbox       = "model bounding box",
                       clim_vars  = "climate variables",
                       theta_clim = "climate bandwidth",
                       theta_geog = "geographic bandwidth",
                       polygon    = "focal polygon",
                       site       = "focal site")
