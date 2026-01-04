.. _release:

***************
Release Process
***************


The release process for Roentgen is as follows:

#. Update the `CHANGELOG.rst` file with the new version number and a summary of changes.
#. Commit the changes to the `CHANGELOG.rst` file and make sure that there are no other uncommitted changes.
#. Determine the new version number according to `semantic versioning <https://semver.org>`.
#. Create a new local git tag for the release using the format `vX.Y.Z`, where `X.Y.Z` is the new version number.
    `git tag vX.Y.Z`
#. Clean out the `dist/` directory to remove any previous build artifacts.

#. Build the distribution packages using `uv build`.
#. Test the distribution packages by installing them in a clean virtual environment and running the test suite.
    `uv run --isolated --no-project --with dist/*.whl pytest`
#. Test that the distribution package can be installed and imported and that the version number is correct.
    `uv run --isolated --no-project --with dist/*.whl python -c "import roentgen; print(roentgen.__version__)"`
#. Upload the distribution packages to TestPyPI using `uv publish --index testpypi --token <TESTPYPI_TOKEN> `.
#. Verify the release on TestPyPI (https://test.pypi.org/project/roentgen/)
#. If you need to make changes, delete the local git tag using `git tag -d vX.Y.Z`, make the necessary changes, and repeat steps above.
#. Push the git tag to the remote repository. This will trigger the GitHub Actions workflow to build and upload the package to PyPI.
