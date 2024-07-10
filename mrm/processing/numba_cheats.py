import atexit
import inspect
from pathlib import Path
import random
import shutil
import string

from dustgoggles.dynamic import define, get_codechild
from hostess.utilities import curry
import numba as nb


def _maybermtree(path):
    try:
        shutil.rmtree(path)
    except (OSError, FileNotFoundError):
        pass


def init_fake_nb_path(remove_on_exit=True):
    fakepre = f"{''.join(random.choices(string.ascii_lowercase, k=15))}"
    fakemod = Path(f".nbtemp/{fakepre}/")
    fakefile = fakemod / f"{fakepre}.py"
    (fakemod / "__pycache__").mkdir(parents=True)
    fakefile.touch()
    if remove_on_exit is True:
        atexit.register(curry(_maybermtree, fakemod))
    return str(fakefile.absolute())


def njit_dyncached(func, fakefile=None, globals_=None, sig=None):
    fakefile = fakefile if fakefile is not None else init_fake_nb_path()
    code = get_codechild(
        compile(inspect.getsource(func), fakefile, "exec")
    )
    defined = define(code, globals_)
    dispatcher = nb.njit(defined, cache=True)
    if sig is None:
        return dispatcher, fakefile
    return dispatcher.compile(sig), fakefile


# DUMB TEST
"""
from hostess.subutils import piped


def thing(c):
    # noinspection PyUnresolvedReferences
    return A + B * c


def create_fakery():
    return njit_dyncached(
        thing, globals_={'A': 2, 'B': 3}, sig=nb.float32(nb.float32)
    )[1]


def roundtrip(fakefile):
    jitted, _ = njit_dyncached(thing, fakefile, sig=nb.float32(nb.float32))
    return int(nb.float32(jitted(4)))


def test_dyncache():
    fakefile = create_fakery()
    assert piped(roundtrip)(fakefile) == 14


test_dyncache()
"""