# Uploading to GitHub

Suggested repository name: **HapReturn**.

Suggested description:

> Recover hifiasm hap1/hap2 provenance from curated merged assemblies.

Create an empty GitHub repository under your own account. Upload the contents
of the HapReturn source directory so that `README.md`, `pyproject.toml` and
`hapreturn/` are at the repository root. Include the hidden `.github/`,
`.gitignore` and `.gitattributes` files; Git is the easiest way to retain them.

```bash
git init
git add .
git commit -m "Initial HapReturn 0.1.0 release"
git branch -M main
git remote add origin https://github.com/YOUR_ACCOUNT/HapReturn.git
git push -u origin main
```

Replace YOUR_ACCOUNT with your GitHub account. The workflow in
`.github/workflows/tests.yml` will run synthetic checks on Linux and Windows.

For a GitHub Release, create tag `v0.1.0` and attach the source archive,
source distribution (`.tar.gz`) and wheel (`.whl`). The wheel installs with
`python -m pip install hapreturn-0.1.0-py3-none-any.whl`. Publishing to GitHub
does not automatically publish the project to PyPI.

Before a public release, review the MIT license, contributor attribution and
the initial-release validation limits described in README.md. No account URL,
DOI or publication metadata has been invented in this package.
