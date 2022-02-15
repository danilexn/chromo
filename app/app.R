# ChroMo
# Web App for unsupervised analysis of nuclear oscillations
#
# MIT License
#
# Copyright (c) 2021 Daniel Leon-Perinan (danilexn)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

packages <-
  c(
    "shiny",
    "shinycssloaders",
    "promises",
    "future",
    "ipc",
    "forecast",
    "dplyr",
    "ggplot2",
    "dplR",
    "lmtest",
    "purrr",
    "tidyr",
    "reshape2",
    "tseries",
    "segclust2d",
    "segmenTier",
    "mclust",
    "lattice",
    "TTR",
    "egg",
    "readxl",
    "RCurl",
    "shinythemes",
    "ggpubr",
    "htmlwidgets",
    "plotly",
    "rapportools",
    "tie",
    "pcalg",
    "VLTimeCausality",
    "network",
    "sna",
    "networkD3",
    "igraph",
    "GGally",
    "PerformanceAnalytics",
    "rjson"
  )

# Load all packages
lapply(packages, require, character.only = TRUE)

# How the futures will be processed
plan(sequential)

# Load essential components
source("ui/ui.R")
source("server/server.R")
source("utils/utils.R")

# ChroMo-specific components
sourceFolder("plots")
sourceFolder("analysis")
sourceFolder("calculators")

# Limit upload size
options(shiny.maxRequestSize = 5 * 1024 ^ 2) # 5 Mb max.

# Attach UI and Server components
shinyApp(ui = ui, server = server, enableBookmarking = "server")