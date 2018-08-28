Python Code Style Guide for `mirage`
===================================

This document serves as a style guide for all `mirage` software development.  Any requested contribution to the `mirage` code repository should be checked against this guide, and any violation of the guide should be fixed before the code is committed to
the `master` branch.  Please refer to the accompanying [`example.py`](https://github.com/spacetelescope/mirage/blob/master/style_guide/example.py) script for a example code that abides by this style guide.

Prerequisite Reading
--------------------

It is assumed that the reader of this style guide has read and is familiar with the following:

- The [PEP8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/)
- The [PEP257 Docstring Conventions Style Guide](https://www.python.org/dev/peps/pep-0257/)
- The [`numpydoc` docstring convention](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)


Workflow
--------

All software development for the `mirage` project should follow a continuous integration workflow.  Before committing any code changes, use `flake8` to check the code against `PEP8` standards.  Also check that your code is conforming to this style guide.


Version Numbers and Tags
------------------------

Any changes pushed to the `master` branch should be tagged with a version number.  The version number convention is `x.y.z`, where

    x = The main version number.  Increase when making incompatible API changes.
    y = The feature number.  Increase when change contains a new feature with or without bug fixes.
    z = The hotfix number. Increase when change only contains bug fixes.

Currently the version number is set in setup.py. Updating the version number should be one of the last things you change prior to merging a
feature branch into `master`. The branch's author and reviewer should agree on the version number prior to merging.

Security
--------

The following items should never be committed in the `mirage` source code or GitHub issues/pull requests:

- Account credentials of any kind (e.g. database usernames and passwords)
- Internal directory structures or filepaths
- Machine names
- Proprietary data


`mirage`-specific Code Standards
------------------------------

`mirage` code shall adhere to the `PEP8` conventions save for the following exceptions:

 - Lines of code need not to be restricted to 79 characters.  However, it is encouraged to break up obnoxiously long lines into several lines if it benefits the overall readability of the code

 Additionally, the code shall adhere to the following special guidelines:

 - Function and class definitions should be placed in alphabetical order in the module
 - It is encouraged to annotate variables and functions using the [`typing`](https://docs.python.org/3/library/typing.html) module (see [PEP 483](https://www.python.org/dev/peps/pep-0483/), [PEP 484](https://www.python.org/dev/peps/pep-0484/), and [PEP 526](https://www.python.org/dev/peps/pep-0526/)).


`mirage`-Specific Documentation Standards
---------------------------------------

`mirage` code shall adhere to the `PEP257` and `numpydoc` conventions.  The following are further recommendations:

- Each module should have at minimum a description, `Authors` and `Use` section.
- Each function/method should have at minimum a description, `Parameters` (if necessary), and `Returns` (if necessary) sections

Acknowledgements
----------------

This style guide as well as [`example.py`](https://github.com/spacetelescope/mirage/blob/master/style_guide/example.py) were adapted from those used by the [`jwql` project](https://github.com/spacetelescope/jwql).