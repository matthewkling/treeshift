# Page-level CSS and <head> tags.

ui_head <- function() {
      tagList(
            ## Report Back and Forward to the server (server/url_state.R).
            ## Shiny's own url_search can't be used: it ignores a value equal
            ## to the one it last sent, which is what Back to the opening
            ## link looks like. As an event, every navigation arrives.
            tags$script(HTML("
  window.addEventListener('popstate', function() {
    if (window.Shiny && Shiny.setInputValue)
      Shiny.setInputValue('history_nav', window.location.search, {priority: 'event'});
  });
")),

            ## The browser and the page that linked here, for the usage log
            ## (server/usage_log.R), which can't read them from the request
            ## behind Shiny Server
            tags$script(HTML("
  $(document).on('shiny:connected', function() {
    Shiny.setInputValue('ts_client', {user_agent: navigator.userAgent,
                                      referrer: document.referrer});
  });
")),

            ## A module tab (ui/right_sidebar.R) opens the modules sidebar
            ## and expands its module.
            tags$script(HTML("
  document.addEventListener('click', function(e) {
    var tab = e.target.closest('.ts-module-tab');
    if (!tab) return;
    var layout = tab.closest('.bslib-sidebar-layout');
    window.bslib.Sidebar.getInstance(layout).toggle('open');
    var item = document.querySelector('#modules_acc > .accordion-item[data-value=\"' +
                                      CSS.escape(tab.dataset.module) + '\"]');
    var button = item && item.querySelector('.accordion-button');
    if (button && button.classList.contains('collapsed')) button.click();
  });
")),

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
  /* Run model button (server/model_settings.R). Solid dark when settings
     are waiting to be applied; when the model is up to date it is disabled and
     drawn as a quiet outline rather than as a faded solid button. */
  /* Under the model settings heading; stays at the top of the sidebar while
     a long panel scrolls beneath it. The sidebar's own background is a
     translucent tint; laid over the page color it is opaque, so the panels
     don't show through. */
  .ts-update { position: sticky; top: 0; z-index: 3;
               background: linear-gradient(var(--_sidebar-bg), var(--_sidebar-bg)),
                           var(--bs-body-bg);
               margin-inline: calc(-1 * var(--_padding, 1rem));
               padding: 0.5rem var(--_padding, 1rem); }
  /* Group headings over the sidebar panels (ui/left_sidebar.R) */
  .ts-group-title { font-size: 0.8rem; font-weight: 600; text-transform: uppercase;
                    letter-spacing: 0.05em; color: var(--_sidebar-fg); margin: 0; }
  /* Same spacing, relative to the text, as the module headers' icons */
  .ts-group-title > i { margin-right: 0.75em; }
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

  /* Module tabs down the map's right edge (ui/right_sidebar.R), in place of
     the sidebar's own toggle while the sidebar is collapsed. The toggle is
     hidden with visibility, not display: bslib ends the open/close animation
     when the toggle's icon finishes turning, which an undisplayed icon never
     does, leaving the sidebar's contents hidden. */
  .ts-map-area { min-width: 0; }
  .ts-module-tabs {
    display: flex; flex-direction: column; gap: 2px;
    flex: 0 0 auto; width: 40px; padding: 4px 0;
    overflow-y: auto;
    background: var(--bs-body-bg);
    border-left: 1px solid var(--bs-border-color);
  }
  .bslib-sidebar-layout:not(.sidebar-collapsed) > .ts-modules-main .ts-module-tabs { display: none; }
  .bslib-sidebar-layout.sidebar-collapsed:has(> .ts-modules-main) > .collapse-toggle { visibility: hidden; }
  .ts-module-tab {
    display: flex; flex-direction: column; align-items: center; gap: 6px;
    padding: 10px 0; border: 0; border-radius: 0;
    background: transparent; color: #495057;
    font-size: 0.8rem; line-height: 1;
  }
  .ts-module-tab span { writing-mode: vertical-rl; white-space: nowrap; }
  .ts-module-tab:hover, .ts-module-tab:focus-visible { background: rgba(0, 0, 0, 0.06); color: #212529; }
"))
      )
}
