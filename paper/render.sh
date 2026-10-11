#!/usr/bin/env bash
set -euo pipefail

# manuscript.md is the authored source. The generated standalone .tex also
# opens in the Codex editor with its built-in compiler and live PDF preview.
# This local command exports a shareable PDF, retaining draft labelling until
# the completion contract in results-requirements.md is satisfied.
paper_root=$(git rev-parse --show-toplevel)
cd "$paper_root"
build_dir="$paper_root/results/paper-build"
output_dir="$paper_root/results/paper"
mkdir -p "$build_dir" "$output_dir"
pandoc paper/manuscript.md --standalone --to latex --output paper/manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$build_dir" \
  paper/manuscript.tex > "$build_dir/export-pass1.log"
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$build_dir" \
  paper/manuscript.tex > "$build_dir/export-pass2.log"
cp "$build_dir/manuscript.pdf" "$output_dir/manuscript-in-progress.pdf"
printf '%s\n' "$output_dir/manuscript-in-progress.pdf"
