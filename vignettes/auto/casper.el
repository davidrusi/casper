(TeX-add-style-hook "casper"
 (lambda ()
    (LaTeX-add-bibliographies
     "references")
    (LaTeX-add-labels
     "sec:intro"
     "tab:sampledata"
     "sec:import"
     "sec:preprocess"
     "fig:plotprocess1"
     "sec:knownvar"
     "eq:rpkm"
     "sec:plots"
     "fig:plot2"
     "fig:plot3")
    (TeX-run-style-hooks
     "natbib"
     "hyperref"
     "graphicx"
     "bm"
     "amssymb"
     "amsmath"
     "latex2e"
     "art12"
     "article"
     "a4paper"
     "12pt")))

