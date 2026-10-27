#!/usr/bin/env python3
"""Build the manuscript PDF with Pandoc and pdfLaTeX."""
from pathlib import Path
import subprocess

here = Path(__file__).resolve().parent
stem = "gaussian_pnp_stability"
subprocess.run([
    "pandoc", str(here / f"{stem}.md"),
    "--from=markdown+tex_math_dollars+tex_math_single_backslash+autolink_bare_uris+pipe_tables",
    "--standalone", "--pdf-engine=pdflatex", "--resource-path=" + str(here),
    "-V", "geometry:margin=1in", "-V", "fontsize:11pt",
    "-V", "colorlinks:true", "-V", "linkcolor:purple", "-V", "urlcolor:purple",
    "-H", str(here / "pdf_header.tex"), "-o", str(here / f"{stem}.pdf"),
], check=True, cwd=here)
print("Built", stem + ".pdf")
