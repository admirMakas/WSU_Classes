(TeX-add-style-hook "wfem"
 (lambda ()
    (LaTeX-add-index-entries
     "rod3"
     "Elements!Rod"
     "rod3t"
     "Elements!Rod!Thermal"
     "beam3"
     "Elements!Beam")
    (LaTeX-add-environments
     "theorem"
     "corollary"
     "definition")
    (LaTeX-add-labels
     "listing:#1"
     "example.txt"
     "command:nodes"
     "sec:bcnconstr"
     "staticloads"
     "sec:actions"
     "sec:elements"
     "el:rod3"
     "el:rod3t"
     "el:beam3"
     "el:inertia")
    (TeX-add-symbols
     '("lstref" 1)
     '("includelisting" 2)
     '("varg" 1)
     '("filename" 1)
     '("variable" 1)
     '("command" 1)
     '("sarg" 1))
    (TeX-run-style-hooks
     "graphicx"
     "pdftex"
     "pdfsync"
     "hyperref"
     "colorlinks"
     "hhline"
     "amsmath"
     "listings"
     "latex2e"
     "art12"
     "article"
     "12pt")))

