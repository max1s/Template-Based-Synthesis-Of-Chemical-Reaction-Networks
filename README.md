# CRNSynthesis

This is a python package for the Syntax-Guide Synthesis of Chemical Reaction Networks to meet formal specifications, using the method described in [Syntax-Guided Optimal Synthesis for Chemical Reaction Networks](https://link.springer.com/chapter/10.1007/978-3-319-63390-9_20)..

A graphical web-interface to this library is provided by [crn-designer](https://github.com/jamesscottbrown/crn-designer).

## Installation

Installing this package requires a Python installation: it should be compatible with both Python 2.7 and Python 3. 

After downloading the code (either as a ZIP file or at the command-line using ``https://github.com/max1s/CRNSynthesis.git``), change directory into the top-level ``CRNSynthesis/`` directory and install using ``python setup.py install``.

This will automatically install required dependencies ([listed here](https://github.com/max1s/CRNSynthesis/blob/master/setup.py#L12)).

You can either install the solver that you wish to use (iSAT or dReach) in a location that is on your ``$PATH`` (for example, by executing the command ``export PATH=$PATH:/some/path/dReal-3.16.09.01-linux/bin`` in the terminal), or provide the path to the correpsonding binary as an argument the the solver caller constructor (``isat_path=`` to ``SolverCallerISAT``, or ``dreal_path=`` to ``SolverCallerDReal``).

## Usage

Examples of usage are contained in the [``examples/``](https://github.com/max1s/CRNSynthesis/tree/master/examples) directory. 

API docs (in HTML) are contained in the ``docs/`` directory. These were generated using the ``rst`` branch of [peterjc's fork of pdoc](https://github.com/peterjc/pdoc.git), using ``pdoc --html --overwrite --docformat=restructuredtext CRNSynthesis``.
