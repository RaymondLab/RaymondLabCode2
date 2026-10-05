import pytest

from raymondlab.core import alphabet


def test_alphabet_has_31_safe_characters():
    assert alphabet.ALPHABET == "23456789abcdefghjkmnpqrstuvwxyz"
    assert alphabet.BASE == 31
    for bad in "01ilo":
        assert bad not in alphabet.ALPHABET


def test_encode3_known_values():
    assert alphabet.encode3(0) == "222"
    assert alphabet.encode3(1) == "223"
    assert alphabet.encode3(31) == "232"
    assert alphabet.encode3(29790) == "zzz"


def test_encode3_rejects_out_of_range():
    with pytest.raises(ValueError):
        alphabet.encode3(-1)
    with pytest.raises(ValueError):
        alphabet.encode3(31**3)


def test_round_trip_every_value():
    seen = set()
    for n in range(31**3):
        s = alphabet.encode3(n)
        assert len(s) == 3
        assert alphabet.decode3(s) == n
        seen.add(s)
    assert len(seen) == 31**3


def test_decode3_rejects_bad_characters_and_length():
    with pytest.raises(ValueError):
        alphabet.decode3("2o2")
    with pytest.raises(ValueError):
        alphabet.decode3("22")
    with pytest.raises(ValueError):
        alphabet.decode3("2222")


@pytest.mark.parametrize("slot,ok", [("7k", True), ("22", True), ("zz", True),
                                     ("o1", False), ("7K", False), ("7", False),
                                     ("7kk", False), ("", False)])
def test_is_slot(slot, ok):
    assert alphabet.is_slot(slot) is ok
