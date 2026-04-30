When answering mathematical questions, write the response to notes/math_output.html using KaTeX, and print open notes/math_output.html to the terminal. The responses should be concatenated, rather than replacing anything that's already in notes/math_output.html.

IMPORTANT: Subagents (Task tool) must NEVER write to or overwrite notes/math_output.html. If a subagent needs to produce HTML output, it should write to a separate file (e.g. notes/literature_review.html) and the main agent will handle any merging.

This repository is built for an academic research project examining the impact of punctuated individual infectiousness profiles on epidemic dynamics. 

The main analysis is run by sourcing `code/run_analysis.R` --- this calls all the essential files in order. The writeup is in `writeup/main.pdf` and is built by `writeup/main.tex`, which calls other `.tex` files in the `writeup/` directory. Some findings, that may be partially out of date, are in `notes/findings*.md`. 