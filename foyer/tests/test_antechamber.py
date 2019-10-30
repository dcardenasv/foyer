
import pytest
import parmed as pmd

from foyer.antechamber import ante_atomtyping
from foyer.tests.utils import get_fn
from foyer.utils.io import has_mbuild
from foyer.exceptions import FoyerError

from distutils.spawn import find_executable
ANTECHAMBER = find_executable('antechamber')

@pytest.mark.skipif(ANTECHAMBER is not None, reason="antechamber is installed")
def test_check_antechamber():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    with pytest.raises(IOError):
        typed = ante_atomtyping(ethane,'gaff')

# TODO: Break this into mult. functions and check atypes
@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
def test_valid_atype_style():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    typed = ante_atomtyping(ethane,'gaff')
    typed = ante_atomtyping(ethane,'gaff2')
    typed = ante_atomtyping(ethane,'amber')
    typed = ante_atomtyping(ethane,'bcc')
    typed = ante_atomtyping(ethane,'sybyl')

@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
def test_invalid_atype_style():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    with pytest.raises(FoyerError,match=r"Unsupported atomtyping style.*"):
        typed = ante_atomtyping(ethane,'gaffff')

@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
def test_ethane_gaff_atypes():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    typed = ante_atomtyping(ethane,'gaff')
    assert sum((1 for at in typed.atoms if at.type == 'c3')) == 2
    assert sum((1 for at in typed.atoms if at.type == 'hc')) == 6

@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
def test_invalid_structure():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    with pytest.raises(FoyerError,match=r"Unknown molecule format.*"):
        typed = ante_atomtyping(ethane.topology,'gaff')

@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_convert_mbuild():
    import mbuild as mb
    ethane = mb.load(get_fn('ethane.mol2'))
    print(ethane.__class__)
    typed = ante_atomtyping(ethane,'gaff')

@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
def test_multiple_molecules():
    ethane = pmd.load_file(get_fn('ethane.mol2'),structure=True)
    # Remove a bond to break the molecule apart
    ethane.bonds.remove(ethane.bonds[0])
    print(ethane.bonds)
    with pytest.raises(FoyerError,match=r"Antechamber requires connectivity.*"):
        typed = ante_atomtyping(ethane,'gaff')

# TODO: Write this test. What to do about error logfile handling? 
#@pytest.mark.skipif(ANTECHAMBER is None, reason="antechamber is not installed")
#def test_ante_error():

