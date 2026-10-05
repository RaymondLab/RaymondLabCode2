"""The 31 characters that rig slots and minted ids are built from.

Ids use a reduced set of digits and letters, so a person who reads or types an
id cannot confuse one character for another.
"""

# Digits 2-9 and lowercase letters, minus 0/1/i/l/o, which are easy to misread.
ALPHABET = "23456789abcdefghjkmnpqrstuvwxyz"
BASE = len(ALPHABET)  # 31

_INDEX = {c: i for i, c in enumerate(ALPHABET)}


def is_slot(s: str) -> bool:
    """True when s is a rig slot: exactly 2 characters, both of them in ALPHABET."""
    return len(s) == 2 and all(c in _INDEX for c in s)


def encode3(n: int) -> str:
    """Write the number n as 3 ALPHABET characters: 0 -> "222", 29790 -> "zzz".

    Raise ValueError when n is outside [0, 31**3).
    """
    if not 0 <= n < BASE**3:
        raise ValueError(f"sequence {n} is outside 0..{BASE**3 - 1}")
    chars = []
    for _ in range(3):
        n, r = divmod(n, BASE)
        chars.append(ALPHABET[r])
    return "".join(reversed(chars))


def decode3(s: str) -> int:
    """The inverse of encode3: read 3 ALPHABET characters back into a number.

    Raise ValueError on a wrong length or a character that is not in ALPHABET.
    """
    if len(s) != 3:
        raise ValueError(f"sequence {s!r} must be exactly 3 characters")
    n = 0
    for c in s:
        if c not in _INDEX:
            raise ValueError(f"character {c!r} in {s!r} is not in the id alphabet")
        n = n * BASE + _INDEX[c]
    return n
