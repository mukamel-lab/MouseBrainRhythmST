# preprocessing/

Transform raw inputs into analysis-ready objects.

- `scripts/` QC, filtering, aggregation, normalization
- `config/` sample sheets + parameters
- `output/` written artifacts (not tracked)
- `logs/` run logs (optional)

Entry point should be a single script in `scripts/` that produces versioned outputs.
