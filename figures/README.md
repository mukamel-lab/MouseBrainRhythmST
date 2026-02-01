# figures/

Figure generation.

- `source/` figure scripts (R) and templates
- `tables/` table scripts and exported CSV/TSV (small)
- `output/` rendered figures (PNG/PDF/SVG). Commit only final, lightweight assets.

Convention: each figure has one script in `source/` that reads from `analysis/results/` and writes to `output/`.
