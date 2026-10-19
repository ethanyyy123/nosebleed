"""Minimal compatibility shim for the removed stdlib module `imp` (only what PyMOL 2.5 needs)."""
import importlib.util, os
def find_module(name, path=None):
    spec = importlib.util.find_spec(name)
    if spec is None:
        raise ImportError(name)
    loc = spec.origin if spec.origin and spec.origin != "namespace" else None
    if loc is None and spec.submodule_search_locations:
        loc = list(spec.submodule_search_locations)[0]
    return (None, os.path.dirname(loc) if loc and loc.endswith("__init__.py") else loc, ("", "", 5))
