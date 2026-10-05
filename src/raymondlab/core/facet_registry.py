"""The metadata validator: checks a metadata dict against the lab's registry of keys.

The registry is the data file standard/facets.json. It lists every key a
sidecar may use, where the key may appear, its type and its allowed values,
and it lists the keys each kind of sidecar must carry. This code holds no key
list and no required list of its own.
"""

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

# A type that mentions one of these words is a structured block. Its inner
# keys are described in prose, so they are not registered one by one.
_BLOCK_WORDS = ("block", "typed", "map", "note")


@dataclass(frozen=True)
class Facet:
    """One registered metadata key: one entry of facets.json."""
    key: str                 # e.g. "sex", "surgery[].eye" ("[]" marks a list of objects)
    scope: tuple[str, ...]   # where the key may appear: "session", "subject", "calib run", ...
    origin: str              # who supplies the value: "captured" (typed by a person), "derived", ...
    type: str                # e.g. "string", "int", "enum", "list[string]", "typed block"
    allowed: tuple[str, ...]  # for type "enum": the permitted values. Empty means any value passes.
    status: str              # "confirmed", "provisional" or "open"


@lru_cache(maxsize=None)
def _read_file(path: Path | None) -> dict:
    """Parse facets.json once per path. Callers must not change the result."""
    if path is None:
        path = files("raymondlab.standard") / "facets.json"
    return json.loads(path.read_text(encoding="utf-8"))


def load_registry(path: Path | None = None) -> dict[str, list[Facet]]:
    """Load facets.json and return {key: [Facet, ...]}.

    path defaults to the file shipped inside this package.
    The "note" field of each entry is for people and is not loaded.
    The same key may appear in several entries with different scopes,
    so each value is a list.
    """
    registry: dict[str, list[Facet]] = {}
    for entry in _read_file(path)["facets"]:
        facet = Facet(
            key=entry["key"],
            scope=tuple(entry["scope"]),
            origin=entry["origin"],
            type=entry["type"],
            allowed=tuple(entry.get("allowed", ())),
            status=entry["status"],
        )
        registry.setdefault(facet.key, []).append(facet)
    return registry


def required_keys(kind: str, path: Path | None = None) -> list[str]:
    """The keys a sidecar of this kind must carry, from the "required" block of facets.json.

    kind is "common" (every session sidecar) or a schema_type such as
    "eye-calibration". Raises KeyError, naming the kind, when the file has no
    list for it.
    """
    required = _read_file(path).get("required", {})
    if kind not in required:
        raise KeyError(f"no required-key list for kind: {kind}")
    return list(required[kind]["keys"])


def flatten(obj: dict, prefix: str = "") -> list[tuple[str, object]]:
    """Turn a nested dict into (dotted_path, value) pairs, in the path style the registry uses.

      {"a": {"b": 1}}     -> [("a.b", 1)]
      {"s": [{"d": 1}]}   -> [("s[].d", 1)]     a list of dicts becomes "[]"
      {"x": [1, 2]}       -> [("x", [1, 2])]    a list of plain values is one leaf
      {"s": []}           -> [("s[]", [])]      an empty list is one leaf, named with "[]"
    """
    pairs: list[tuple[str, object]] = []
    for name, value in obj.items():
        path = f"{prefix}.{name}" if prefix else name
        if isinstance(value, dict):
            pairs.extend(flatten(value, path))
        elif isinstance(value, list) and not value:
            pairs.append((f"{path}[]", value))
        elif isinstance(value, list) and all(isinstance(item, dict) for item in value):
            for item in value:
                pairs.extend(flatten(item, f"{path}[]"))
        else:
            pairs.append((path, value))
    return pairs


def _ancestors(path: str) -> list[str]:
    """The shorter paths above this one: "a.b.c" gives ["a.b", "a"]."""
    parts = path.split(".")
    return [".".join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def _in_scope(entries: list[Facet], scopes: set[str]) -> list[Facet]:
    return [f for f in entries if scopes.intersection(f.scope)]


def validate(obj: dict, scopes: set[str], registry: dict[str, list[Facet]] | None = None) -> list[str]:
    """Check every key in obj against the registry. Return one error string per problem; [] means valid.

    scopes are the registry's scope names, exactly as written there:
    "session", "subject", "block", "signal", "device", "calib run", ...
    A sidecar mixes several of them.

    A flattened path is accepted when one of these holds:
      - it is registered exactly, and one of its entries has a scope in scopes;
      - a registered ancestor path has a type that mentions "block", "typed",
        "map" or "note" (a structured block whose inner keys are described in
        prose, not registered one by one);
      - it is an empty list "X[]" and some registered key starts with "X[]."
        and has an entry with a scope in scopes. Nothing inside an empty list
        can be wrong, and a new sidecar starts with empty lists.
    If the matching entry's type is "enum" and its allowed list is not empty,
    str(value) must be in that list. (str() because a channel number is
    stored as the integer 1 but registered as "1".)
    A value of None skips that list check: the standard allows "null" for keys such as
    state.qc[].reason (null when the run is usable). The key must still be registered
    and in scope.
    """
    if registry is None:
        registry = load_registry()
    errors = []
    for path, value in flatten(obj):
        matches = _in_scope(registry.get(path, []), scopes)
        if matches:
            if value is not None and not _enum_accepts(matches, value):   # None means "no value"
                errors.append(f"value not allowed for {path}: {value!r}")
            continue
        if any(f.type and any(w in f.type for w in _BLOCK_WORDS)
               for ancestor in _ancestors(path) for f in registry.get(ancestor, [])):
            continue
        if value == [] and path.endswith("[]") and _has_inner_key(registry, path, scopes):
            continue
        errors.append(f"unregistered key: {path}")
    return errors


def _enum_accepts(matches: list[Facet], value: object) -> bool:
    """True when some matching entry is not a limited enum, or lists str(value)."""
    return any(f.type != "enum" or not f.allowed or str(value) in f.allowed for f in matches)


def _has_inner_key(registry: dict[str, list[Facet]], list_path: str, scopes: set[str]) -> bool:
    """True when a registered key inside this list has an entry with a scope in scopes."""
    inner = list_path + "."
    return any(key.startswith(inner) and _in_scope(entries, scopes) for key, entries in registry.items())


def _is_present(obj: object, parts: list[str]) -> bool:
    """Walk the path parts into obj. See check_required for the rules."""
    if not isinstance(obj, dict):
        return False
    name, rest = parts[0], parts[1:]
    if name.endswith("[]"):
        items = obj.get(name[:-2])
        if not isinstance(items, list):
            return False
        return all(_is_present(item, rest) for item in items) if rest else True
    value = obj.get(name)
    if value is None:
        return False
    return _is_present(value, rest) if rest else True


def check_required(obj: dict, keys: list[str]) -> list[str]:
    """Check that every key in `keys` is present in obj. Return one error per missing key; [] means all present.

    A path is dotted: "state.lifecycle[]" means obj["state"]["lifecycle"].
    A plain path "weight" must be a key of obj with a value that is not None.
    A path "runs[]" must be a list (an empty list passes).
    A path "runs[].run_index" must be a list whose every entry has a "run_index" key (an empty list passes).
    A missing parent counts as the path missing.
    """
    return [f"missing required key: {key}" for key in keys if not _is_present(obj, key.split("."))]
