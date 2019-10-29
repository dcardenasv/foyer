from __future__ import division

import os
import sys
import warnings
from distutils.spawn import find_executable
from subprocess import PIPE, Popen

from foyer.utils.tempdir import temporary_directory
from foyer.utils.tempdir import temporary_cd

ANTECHAMBER = find_executable('antechamber')

def antechamber_atomtyping():

    _check_antechamber(ANTECHAMBER)

    # Check valid ff name

    # Confirm single connected molecule
    # Also possible that connectivity info doesn't exist

    # Work within a temporary directory
    # to clean up after antechamber
    with temporary_directory() as tmpdir:
        with temporary_cd(tmpdir):
            # Save the existing molecule to file
            molecule.save('ante_in.mol2')
            # Call antechamber
            command = ( 'antechamber -i ante_in.mol2 -fi mol2 '
                         '-o ante_out.mol2 fo mol2 -at '
                         + forcefield_type + '-c ' + total_charge )

            proc = Popen(command, stdout=PIPE, stderr=PIPE,
                         universal_newlines=True, shell=True)

            out, err = proc.communicate()

            # Error handling here
            if 'ERROR' in out or proc.returncode != 0:
                _antechamber_error(out,err)


def _antechamber_error(out, err):
    """Log antechamber output to file. """
    with open('../ante_errorlog.txt', 'w') as log_file:
        log_file.write(out)
    raise RuntimeError("Antechamber failed. See 'ante_errorlog.txt'")

def _check_antechamber(ANTECHAMBER):
    if not ANTECHAMBER:
        msg = "Antechamber not found."
        raise IOError(msg)
