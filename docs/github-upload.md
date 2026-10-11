# Build and release

From the repository root:

```bash
python -m pip install build
python -m build
python -m unittest discover -s tests -v
```

The `dist/` directory contains the wheel and source distribution. Attach these
files to a versioned GitHub release. Install a downloaded wheel with:

```bash
python -m pip install hapreturn-0.5.0-py3-none-any.whl
```

The GitHub Actions workflow runs the test suite on Linux and Windows.
