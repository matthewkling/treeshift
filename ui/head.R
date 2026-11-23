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

            ## The wash over a stale map (output$map_pending,
            ## server/model_settings.R) goes as soon as Run model applies the
            ## settings, before the map is redrawn. A copy holds it from the
            ## click until the app next goes idle, by which time the new map
            ## has been sent, or its error shown.
            tags$script(HTML("
  (function() {
    var hold = null, busy = false;
    $(document).on('click', '#update_model', function() {
      var area = document.querySelector('.ts-map-area');
      if (hold || !area) return;
      hold = document.createElement('div');
      hold.className = 'ts-map-stale';
      area.appendChild(hold);
      busy = false;
    });
    $(document).on('shiny:busy', function() { if (hold) busy = true; });
    $(document).on('shiny:idle', function() {
      if (hold && busy) { hold.remove(); hold = null; }
    });
  })();
")),

            ## Progress messages sit in the middle of the map rather than in
            ## the corner with other notifications. They stay in Shiny's
            ## notification panel, where Shiny finds them to update and
            ## remove them, and are pinned to the map's center (CSS below),
            ## which this keeps current as the sidebars open and close. While
            ## the app is busy for more than a moment with something that has
            ## no message of its own (a model run's comes from withProgress(),
            ## server/usage_log.R), a generic one in the same style.
            tags$script(HTML("
  (function() {
    function center() {
      var m = document.getElementById('map');
      var r = m && m.getBoundingClientRect();
      var root = document.documentElement.style;
      // (The map has no size while another tab is shown: the window's center)
      root.setProperty('--ts-map-cx', r && r.width ? (r.left + r.width / 2) + 'px' : '50vw');
      root.setProperty('--ts-map-cy', r && r.height ? (r.top + r.height / 2) + 'px' : '50vh');
    }
    $(document).on('shiny:connected', function() {
      center();
      new ResizeObserver(center).observe(document.getElementById('map'));
      window.addEventListener('resize', center);
      document.addEventListener('shown.bs.tab', center);
    });

    var timer = null;
    function panel() {
      var p = document.getElementById('shiny-notification-panel');
      if (!p) {
        p = document.createElement('div');
        p.id = 'shiny-notification-panel';
        document.body.appendChild(p);
      }
      return p;
    }
    $(document).on('shiny:busy', function() {
      clearTimeout(timer);
      timer = setTimeout(function() {
        if (document.getElementById('ts-busy')) return;
        var n = document.createElement('div');
        n.id = 'ts-busy';
        n.className = 'shiny-notification shiny-progress-notification';
        n.innerHTML = '<div class=\"progress active\"><div class=\"progress-bar\"></div></div>' +
                      '<div class=\"progress-text\"><span class=\"progress-message\">Updating</span></div>';
        panel().appendChild(n);
      }, 700);
    });
    $(document).on('shiny:idle', function() {
      clearTimeout(timer);
      var n = document.getElementById('ts-busy');
      if (n) n.remove();
    });
  })();
")),

            ## A link to an About section (about_link(), ui/about.R) opens the
            ## About tab and scrolls to the section once the tab is shown.
            tags$script(HTML("
  document.addEventListener('click', function(e) {
    var link = e.target.closest('.ts-about-link');
    if (!link) return;
    e.preventDefault();
    var target = document.getElementById(link.getAttribute('href').slice(1));
    var tab = document.querySelector('#main_nav a[data-value=\"about\"]');
    if (!target || !tab) return;
    var go = function() { target.scrollIntoView(); };
    if (tab.classList.contains('active')) return go();
    tab.addEventListener('shown.bs.tab', go, {once: true});
    bootstrap.Tab.getOrCreateInstance(tab).show();
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

            ## Plot hover readouts (plot_tooltip(), app.R)
            tags$style(HTML("
  .ts-plot-tooltip {
    position: absolute; z-index: 10; pointer-events: none; width: max-content; max-width: 320px;
    background: rgba(255, 255, 255, 0.95); border: 1px solid #bbb; border-radius: 4px;
    padding: 4px 8px; font-size: 12px; line-height: 1.35;
  }
")),

            tags$style(HTML("
  /* Only tiles belonging to the provider tile layer */
  #map .leaflet-layer.basemap-gray img.leaflet-tile {
    filter: grayscale(100%) brightness(100%) opacity(35%);
    -webkit-filter: grayscale(100%) brightness(100%) opacity(35%);
  }
")),

            tags$style(HTML("
  /* Warnings: stale results, settings too large to run, layers withdrawn.
     Colored enough to register as something to read, without the solid fill
     of an alert box. */
  .ts-notice {
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

  /* Unmet preconditions -- no focal site, no polygon, no plots in range --
     are a normal state on the way to using the tool, so they're muted rather
     than flagged: a module's own notice (server/module_status.R) and the
     validate() messages Shiny shows in place of an output. */
  .ts-module-notice,
  .shiny-output-error-validation {
    display: flex;
    gap: 0.5rem;
    align-items: baseline;
    color: #495057;
    background-color: #f1f3f5;
    border: none;
    border-radius: 4px;
    padding: 0.5rem 0.75rem;
    margin: 0.5rem 0;
    font-size: 0.85rem;
    line-height: 1.4;
    font-style: normal;
    text-align: left;
  }
  /* A validate() message takes over its output's own element, so it would
     otherwise keep a plot's full height */
  .shiny-output-error-validation { height: auto !important; }
  .ts-module-notice > i { color: #868e96; }
")),

            tags$style(HTML("
  /* The map's legend panel (helpers/cartography.R) and value readout
     (server/map_readout.R), in one style */
  .ts-legend, .ts-map-readout {
    padding: 8px 10px; border-radius: 5px;
    background: rgba(255, 255, 255, 0.9); box-shadow: 0 0 15px rgba(0, 0, 0, 0.2);
    font: 12px/1.35 'Helvetica Neue', Arial, Helvetica, sans-serif; color: #333;
  }
  .ts-legend { width: 220px; }
  .ts-legend-section + .ts-legend-section {
    border-top: 1px solid #dee2e6; margin-top: 8px; padding-top: 8px;
  }
  .ts-legend-title { font-weight: 700; margin-bottom: 6px; }
  .ts-legend-muted { font-weight: 400; color: #6c757d; }
  /* Inset so the end labels, centered on their ticks, stay in the box */
  .ts-legend-bar { height: 10px; margin: 0 12px; border-radius: 2px; }
  .ts-legend-ticks { position: relative; height: 18px; margin: 0 12px; }
  .ts-legend-ticks span {
    position: absolute; top: 5px; transform: translateX(-50%);
    white-space: nowrap; font-size: 11px;
  }
  .ts-legend-ticks span::before {
    content: ''; position: absolute; left: 50%; top: -5px; height: 4px;
    border-left: 1px solid #666;
  }
  .ts-legend-row { display: flex; align-items: center; gap: 6px; margin: 2px 0; }
  /* Outlined, so a pale swatch (the gray 'other') shows against the panel */
  .ts-legend-row i {
    flex: 0 0 auto; width: 14px; height: 12px;
    box-shadow: inset 0 0 0 1px rgba(0, 0, 0, 0.15);
  }
  .ts-legend-circles { display: flex; justify-content: space-between; margin-bottom: 4px; }
  .ts-legend-circles > div {
    display: flex; flex-direction: column; align-items: center; gap: 1px; font-size: 11px;
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

  /* Progress messages, pinned to the map's center (see the script above).
     The generic busy message gives way to a model run's own. */
  #shiny-notification-panel .shiny-notification:has(.shiny-progress-notification),
  #ts-busy {
    position: fixed;
    left: var(--ts-map-cx, 50vw);
    top: var(--ts-map-cy, 50vh);
    transform: translate(-50%, -50%);
    width: min(360px, 90vw);
    margin: 0;
    opacity: 1;
  }
  #shiny-notification-panel:has(.shiny-notification:not(#ts-busy) .shiny-progress-notification) #ts-busy { display: none; }
  /* A bar with no value set (withProgress(value = NULL), and the busy
     message) shows that work is going on, not how far along it is: full
     width, with moving stripes. A bar given a value fills as usual. Shiny
     hides the bar until it has a value, with an inline style. */
  .shiny-progress-notification .progress:has(.progress-bar:not([style*='width'])) {
    display: flex !important;
  }
  .shiny-progress-notification .progress-bar:not([style*='width']) {
    width: 100%;
    background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.25) 25%, transparent 25%,
                                      transparent 50%, rgba(255, 255, 255, 0.25) 50%,
                                      rgba(255, 255, 255, 0.25) 75%, transparent 75%, transparent);
    background-size: 1rem 1rem;
    animation: progress-bar-stripes 1s linear infinite;
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
  /* The copy that outlasts it (see the script above) isn't stacked on it */
  .ts-map-area:has(#map_pending .ts-map-stale) > .ts-map-stale { display: none; }

  /* The value under the cursor (server/map_readout.R), above the scale bar
     and styled like the legend (above). Below the stale wash, so it fades with
     the map it reads. */
  .ts-map-readout {
    position: absolute; left: 10px; bottom: 36px; z-index: 1000;
    max-width: 260px; pointer-events: none;
  }
  .ts-map-readout-where { font-weight: 700; }

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
