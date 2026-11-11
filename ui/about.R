# About tab: background and documentation. Each section is a Markdown file in
# docs/about/, shown in filename order and starting with a "## " heading, so
# the text can be edited without touching R. The list of sections at the top
# links to each one.

ui_about <- function(dir = "docs/about") {
      files <- sort(list.files(dir, pattern = "\\.md$", full.names = TRUE))
      ids <- sub("^[0-9]+_", "", tools::file_path_sans_ext(basename(files)))
      text <- lapply(files, readLines, encoding = "UTF-8", warn = FALSE)
      headings <- vapply(text, function(x) sub("^##\\s+", "", x[grepl("^##\\s", x)][1]), "")

      div(
            id = "about-page",
            class = "container py-4",
            style = "max-width: 820px;",
            ## Monochrome links, matching the header, in place of the
            ## theme's green
            tags$style("#about-page section h2 { font-size: 1.4rem; margin-top: 1.5rem; }
                        #about-page a { color: #212529; text-decoration-color: #adb5bd; }
                        #about-page a:hover { text-decoration-color: #212529; }"),
            tags$h1("About TreeShift", class = "h2 mb-3"),
            tags$nav(class = "mb-4",
                     tags$ul(class = "list-inline",
                             lapply(seq_along(ids), function(i)
                                   tags$li(class = "list-inline-item me-3",
                                           tags$a(href = paste0("#about-", ids[i]), headings[i]))))),
            lapply(seq_along(ids), function(i)
                  tags$section(id = paste0("about-", ids[i]), class = "mb-4",
                               HTML(commonmark::markdown_html(paste(text[[i]], collapse = "\n"),
                                                                    extensions = "autolink"))))
      )
}
