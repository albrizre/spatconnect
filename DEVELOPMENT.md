# Development workflow

A typical local workflow is:

```r
install.packages(c("devtools", "roxygen2", "testthat"))

devtools::document()
devtools::test()
devtools::check()
```

The public API is intentionally small: `global_connectivity()` and
`local_connectivity()`. The package uses the terminology `superlevel` and
`sublevel` throughout.

Before the public GitHub release, replace `YOUR_GITHUB_USERNAME` in `README.md`.

For GitHub:

```bash
git init
git add .
git commit -m "Initial spatconnect package"
git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/spatconnect.git
git push -u origin main
```

The included GitHub Actions workflow runs `R CMD check` on pushes and pull
requests.
