# About tab: background and documentation. Each section is a Markdown file in
# docs/about/, shown in filename order and starting with a "## " heading;
# "### " headings within a file are its subsections. The text can be edited
# without touching R.
#
# The tab fills the area below the header (app.R) and scrolls within it, so
# the header stays in view. A contents list in a left-hand column stays in
# view too. It links to every section and subsection, and marks the one being
# read. Only the current section's subsections are shown, so the list stays
# short as the page grows. Below the md breakpoint the column stacks above the
# text and doesn't stick.

about_slug <- function(x) gsub("^-+|-+$", "", gsub("[^a-z0-9]+", "-", tolower(x)))

# One Markdown file as HTML, with ids on its subsection headings, plus the
# entries for the contents list
about_section <- function(file) {
      id <- paste0("about-", sub("^[0-9]+_", "", tools::file_path_sans_ext(basename(file))))
      html <- commonmark::markdown_html(paste(readLines(file, encoding = "UTF-8", warn = FALSE),
                                              collapse = "\n"),
                                        extensions = "autolink")
      title <- sub("^.*?<h2>(.*?)</h2>.*$", "\\1", html)
      h3 <- regmatches(html, gregexpr("<h3>.*?</h3>", html))[[1]]
      sub_titles <- sub("^<h3>(.*)</h3>$", "\\1", h3)
      # (paste0() of a zero-length vector still returns one string)
      sub_ids <- if (length(h3)) paste0(id, "-", about_slug(gsub("<[^>]+>", "", sub_titles)))
      for (k in seq_along(h3)) {
            html <- sub(h3[k], sprintf('<h3 id="%s">%s</h3>', sub_ids[k], sub_titles[k]),
                        html, fixed = TRUE)
      }
      list(id = id, title = title, sub_ids = sub_ids, sub_titles = sub_titles, html = html)
}

ui_about <- function(dir = "docs/about") {
      sections <- lapply(sort(list.files(dir, pattern = "\\.md$", full.names = TRUE)),
                         about_section)

      toc <- tags$nav(
            id = "about-toc", class = "nav flex-column",
            lapply(sections, function(s) tagList(
                  tags$a(class = "nav-link", href = paste0("#", s$id), HTML(s$title)),
                  if (length(s$sub_ids))
                        tags$nav(class = "nav flex-column",
                                 lapply(seq_along(s$sub_ids), function(k)
                                       tags$a(class = "nav-link", href = paste0("#", s$sub_ids[k]),
                                              HTML(s$sub_titles[k]))))
            )))

      ## The scrolling pane spans the full width, so its scrollbar is at the
      ## window edge; the text is centered within it
      div(
            id = "about-scroll",
            class = "html-fill-item",
            style = "overflow-y: auto;",
            div(
                  id = "about-page",
                  class = "container-xxl py-4",
                  ## Monochrome, matching the header: the current entry in the
                  ## contents is dark with a rule beside it, like the active tab
                  ## HTML(): tags$style() escapes plain text, turning the ">"
                  ## combinator into "&gt;" and voiding that rule
                  tags$style(HTML("
                    /* bslib pads fillable tabs; the pane runs to the window edges */
                    .tab-pane:has(> #about-scroll) { padding: 0; }
                    #about-page section h2 { font-size: 1.4rem; margin-top: 1.5rem; }
                    #about-page section, #about-page section h3 { scroll-margin-top: 1rem; }
                    #about-page section h3 { font-size: 1.05rem; font-weight: 600; margin-top: 1.25rem; }
                    #about-page section a { color: #212529; text-decoration-color: #adb5bd; }
                    #about-page section a:hover { text-decoration-color: #212529; }
                    #about-toc .nav-link { color: #6c757d; padding: 0.2rem 0.75rem;
                                           border-left: 2px solid #dee2e6; font-size: 0.95rem; }
                    #about-toc .nav-link:hover { color: #212529; }
                    #about-toc .nav-link.active { color: #212529; font-weight: 600;
                                                  border-left-color: #212529; }
                    #about-toc .nav .nav-link { padding-left: 1.5rem; font-size: 0.85rem; }
                    /* !important: bslib sets .nav:not(.nav-hidden) to flex !important */
                    #about-toc .nav { display: none !important; }
                    #about-toc > .nav-link.active + .nav { display: flex !important; }
                    @media (min-width: 768px) {
                      #about-toc { position: sticky; top: 1rem; }
                    }")),
                  ## Marks the entry being read: the last heading scrolled up to the
                  ## top of the pane, or the last one once the pane can't scroll
                  ## further. (Bootstrap's ScrollSpy marks whichever heading most
                  ## recently entered its band, which with short subsections is
                  ## often the one after the one being read.) While the About tab
                  ## is hidden this does nothing.
                  tags$script(HTML("
                    (function() {
                      var queued = false;
                      function mark() {
                        queued = false;
                        var pane = document.getElementById('about-scroll');
                        var toc = document.getElementById('about-toc');
                        if (!pane || !toc || toc.offsetParent === null) return;
                        var links = Array.from(toc.querySelectorAll('a.nav-link'));
                        var atEnd = pane.scrollTop + pane.clientHeight >= pane.scrollHeight - 2;
                        var below = pane.getBoundingClientRect().top + 60;
                        var current = links[0];
                        links.forEach(function(a) {
                          var target = document.getElementById(a.getAttribute('href').slice(1));
                          if (target && (atEnd || target.getBoundingClientRect().top <= below)) current = a;
                        });
                        links.forEach(function(a) { a.classList.remove('active'); });
                        current.classList.add('active');
                        var list = current.parentElement;
                        if (list !== toc) list.previousElementSibling.classList.add('active');
                      }
                      function queue() { if (!queued) { queued = true; requestAnimationFrame(mark); } }
                      // (scroll doesn't bubble; capturing catches the pane's)
                      document.addEventListener('scroll', queue, { passive: true, capture: true });
                      window.addEventListener('resize', queue);
                      window.addEventListener('load', queue);
                      document.addEventListener('shown.bs.tab', queue);
                    })();")),
                  div(class = "row",
                      div(class = "col-md-3 mb-4", toc),
                      div(class = "col-md-9", style = "max-width: 820px;",
                          tags$h1("About TreeShift", class = "h2 mb-3"),
                          lapply(sections, function(s)
                                tags$section(id = s$id, class = "mb-4", HTML(s$html))))))
      )
}
