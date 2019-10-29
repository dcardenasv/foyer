from __future__ import division

import os
import sys
import warnings

import parmed as pmd
import networkx as nx

from distutils.spawn import find_executable
from subprocess import PIPE, Popen
from foyer.utils.tempdir import temporary_directory
from foyer.utils.tempdir import temporary_cd

from foyer.exceptions import FoyerError
from foyer.utils.io import import_, has_mbuild

ANTECHAMBER = find_executable('antechamber')

def ante_atomtyping(molecule, atype_style, net_charge=0.0,
                           multiplicity=1):
    """Perform atomtyping by calling antechamber

    Parameters
    ----------
    molecule : parmed.Structure or mbuild.Compound
        Molecular structure to perform atomtyping on
    atype_style : str
        Style of atomtyping. Options include 'gaff', 'gaff2',
        'amber', 'bcc', 'sybyl'.
    net_charge : float, optional, default=0.0
        Net charge of the molecule
    multiplicity : int, optional, default=1
        Multiplicity, 2S + 1

    Returns
    -------
    typed_molecule : parmed.Structure
        The molecule with antechamber atomtyping applied
    """
    _check_antechamber(ANTECHAMBER)

    # Check valid atomtype name
    supported_atomtypes = ['gaff', 'gaff2', 'amber', 'bcc', 'sybyl' ]
    if atype_style not in supported_atomtypes:
        raise FoyerError( 'Unsupported atomtyping style requested. '
            'Please select from {}'.format(supported_atomtypes))

    # Check parmed structure. Convert if from mbuild
    if not isinstance(molecule, pmd.Structure):
        if has_mbuild:
            mb = import_('mbuild')
            if isinstance(molecule, mb.Compound):
                molecule = molecule.to_parmed()
        else:
            raise FoyerError('Unknown molecule format: {}\n'
                             'Supported formats are: '
                             '"parmed.Structure" and '
                             '"mbuild.Compound"'.format(molecule))

    # Confirm single connected molecule
    # Also possible that connectivity info is missing. Either case is a problem.
    graph_edges = []
    for bond in molecule.bonds:
        graph_edges.append([bond.atom1.idx,bond.atom2.idx])
    bond_graph = nx.Graph()
    bond_graph.add_edges_from(graph_edges)
    if not nx.is_connected(bond_graph):
        raise ValueError("Antechamber requires connectivity information and "
                         "only supports single molecules (i.e., all atoms "
                         "in the molecule are connected by bonds.")

    # Work within a temporary directory
    # to clean up after antechamber
    with temporary_directory() as tmpdir:
        with temporary_cd(tmpdir):
            # Save the existing molecule to file
            molecule.save('ante_in.mol2')
            # Call antechamber
            command = ( 'antechamber -i ante_in.mol2 -fi mol2 '
                         '-o ante_out.mol2 -fo mol2 ' +
                         '-at ' + atype_style + ' ' +
                         '-nc ' + str(net_charge) + ' ' +
                         '-m ' + str(multiplicity) +  ' ' +
                         '-s 2' )

            proc = Popen(command, stdout=PIPE, stderr=PIPE,
                         universal_newlines=True, shell=True)

            out, err = proc.communicate()

            # Error handling here
            if 'ERROR' in out or proc.returncode != 0:
                _antechamber_error(out,err)

            # Now read in the mol2 file with atomtyping
            typed_molecule = pmd.load_file('ante_out.mol2',structure=True)

    # And return it
    return typed_molecule

def _antechamber_error(out, err):
    """Log antechamber output to file. """
    with open('../ante_errorlog.txt', 'w') as log_file:
        log_file.write(out)
    raise RuntimeError("Antechamber failed. See 'ante_errorlog.txt'")

def _check_antechamber(ANTECHAMBER):
    if not ANTECHAMBER:
        msg = "Antechamber not found."
        raise IOError(msg)
