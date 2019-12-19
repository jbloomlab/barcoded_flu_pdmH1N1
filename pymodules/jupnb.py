"""Utilities for use with Jupyter notebooks."""


import os
import subprocess

import papermill


def run_nb_to_html(input_nb, output_nb, parameters):
    """Runs parameterized Jupyter notebook and generates HTML output.

    `input_nb` is parameterized with `parameters` using
    `papermill <https://papermill.readthedocs.io/>_, and then
    run in the current directory to create `output_nb`. That output
    notebook is then converted to HTML to create a file in the same
    location as `output_nb` but with extension ``.html``.

    Parameters
    ----------
    input_nb : str
        Name of Jupyter notebook to parameterize and run.
    output_nb : str
        Name of output Jupyter notebook.
    parameters : dict
        Parameters for Jupyter notebook.

    """
    if os.path.splitext(input_nb)[1] != '.ipynb':
        raise ValueError(f"`input_nb` lacks extension `.ipynb`: {input_nb}")
    if os.path.splitext(output_nb)[1] != '.ipynb':
        raise ValueError(f"`output_nb` lacks extension `.ipynb`: {output_nb}")

    papermill.execute_notebook(
            input_path=input_nb,
            output_path=output_nb,
            cwd=os.getcwd(),
            parameters=parameters,
            )

    # https://github.com/ipython-contrib/jupyter_contrib_nbextensions/issues/901
    subprocess.check_call(['jupyter', 'nbconvert', output_nb,
                           '--to', 'html_embed', '--template', 'toc2'])

    output_html = os.path.splitext(output_nb)[0] + '.html'
    assert os.path.isfile(output_html), f"failed to generate {output_html}"
