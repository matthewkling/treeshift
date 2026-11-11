# Page-level CSS and <head> tags.

ui_head <- function() {
      tagList(
            tags$style(HTML("
  /* Only tiles belonging to the provider tile layer */
  #map .leaflet-layer.basemap-gray img.leaflet-tile {
    filter: grayscale(100%) brightness(100%) opacity(35%);
    -webkit-filter: grayscale(100%) brightness(100%) opacity(35%);
  }
")),

            tags$style(HTML("
  /* Notices come from two mechanisms that cannot share markup: inline divs in
     the sidebar, and validate() messages, which Shiny renders as plain text in
     place of an output. Styling them together is what makes them read as one
     kind of thing.

     Pitched between the two extremes it replaced: colored enough to register as
     something to read rather than blending into body text, but without the
     solid fill of an alert box. Most of these report an unmet precondition -- no polygon
     selected, no focal site set -- which is a normal state on the way to using
     the tool, and a full-width warning flag overstates it. */
  .ts-notice,
  .shiny-output-error-validation {
    display: block;
    color: #9a3412;
    background-color: transparent;
    border: none;
    border-left: 3px solid #c2410c;
    border-radius: 0;
    padding: 0.35rem 0 0.35rem 0.7rem;
    margin: 0.5rem 0;
    font-size: 0.85rem;
    line-height: 1.4;
    font-style: normal;
    text-align: left;
  }
")),

            tags$style(HTML("
  /* Raster legend: title beside the color bar rather than stacked above it.
     A multi-line title over a narrow gradient leaves a lot of dead width, and
     the title is the part that wraps. Targets .info.legend, which addLegend()
     produces; circleLegend() builds its own .leaflet-control container and is
     unaffected. */
  /* Gradient legends only: a discrete legend's swatch rows are themselves the
     container's children, so laying them out in a row strings every category
     across the screen. Those opt out via .legend-discrete. */
  #map .info.legend:not(.legend-discrete) {
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 10px;
  }
  #map .info.legend > div:first-child {
    max-width: 11em;
    min-width: 6em;
    line-height: 1.2;
    text-align: left;
  }
")),

            tags$style(HTML("
  /* Remove padding on the main content wrapper(s) */
  .bslib-page-sidebar main,
  .bslib-page-sidebar .main {
    padding: 0 !important;
  }
")),

            tags$head(
                  tags$style(HTML("
              /* Make the tooltip wider */
              .tooltip-inner {
                  max-width: 350px !important;
                  text-align: left;
                  text-transform: none !important;
                  line-height: 1.2 !important;
                  padding: 8px 12px;
                  font-size: 0.85rem;
              }

              /* Subtle styling for definition list descriptions */
              .tooltip-inner dd {
                  color: #aaa;
                  margin-left: 0;
                  margin-bottom: 0.5em;
              }

              /* Optional: make terms stand out more */
              .tooltip-inner dt {
                  color: #fff;
                  font-weight: 600;
              }
      "))
            ),

            tags$style(HTML("
  /* Update map button (server/model_settings.R). Solid dark when settings
     are waiting to be applied; when the map is up to date it is disabled and
     drawn as a quiet outline rather than as a faded solid button. */
  .ts-update { margin-bottom: 0.5rem; }
  #update_model {
    /* the header's dark gray; Flatly's own btn-dark is a gray-green */
    background-color: #212529;
    border-color: #212529;
    color: #fff;
  }
  #update_model:not(:disabled):hover { background-color: #000; border-color: #000; }
  #update_model:disabled {
    background-color: transparent;
    border-color: #ced4da;
    color: #6c757d;
    opacity: 1;
  }
  .ts-pending-note {
    font-size: 0.85rem;
    color: #495057;
    margin-top: 0.4rem;
  }

  /* Wash over the map while settings are pending. Above the legends
     (Leaflet's controls sit at z-index 1000) so they fade with the layers
     they describe; clicks pass through to the map. */
  .ts-map-stale {
    position: absolute;
    inset: 0;
    z-index: 1001;
    background: rgba(255, 255, 255, 0.65);
    pointer-events: none;
  }
"))
      )
}
