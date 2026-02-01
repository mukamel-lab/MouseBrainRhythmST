# shiny_apps/

Shiny apps for interactive exploration.

- `app/` main `app.R` (or `ui.R`/`server.R`)
- `modules/` reusable Shiny modules
- `www/` static assets (CSS, images)
- `config/` app configuration (paths, feature flags)

Recommendation: apps should load small demo data by default and optionally mount full data via a configurable path.
