# Spatial scopes the model-evaluation module can score over. These are the
# EVALUATION set, never the analog pool: predictions are always built from
# every available analog, and the scope only decides which held-out plots the
# metrics are computed from. "Best bandwidth here" and "best bandwidth overall"
# are different answers, and this is what lets the module give the first one.
#
# Deliberately not the same three scopes as model_delta_region: a focal site is
# a point, and a point is not an evaluation set. The meaningful site-scoped
# quantity is the plots within the geographic cutoff of the click -- the same
# neighborhood the model draws its analogs from.
# The areas a module can summarize over, named once so a code means the same
# thing everywhere. This previously diverged: "site" meant the analog
# neighborhood in the evaluation and climate-relationship modules and the single
# clicked cell in the composition and climate-model modules. Nothing shared them
# yet, so nothing broke -- but the first helper used by both would have silently
# done the wrong one.
#
# "site" is the neighborhood; "cell" is the one clicked cell.
area_names <- c("Model bounding box"      = "domain",
                "Focal polygon"           = "polygon",
                "Buffered polygon"        = "buffered",
                "Focal site neighborhood" = "site",
                "Focal site"              = "cell")

# Modules take the subset that is meaningful for them. Evaluation and climate
# relationships exclude the single cell: one plot is not an evaluation set, and
# one cell is not a relationship.
eval_scope_names <- area_names[area_names != "cell"]
