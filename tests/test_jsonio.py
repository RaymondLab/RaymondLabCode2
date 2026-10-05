import json
from pathlib import Path

from raymondlab.core import jsonio


def test_write_read_write_gives_identical_bytes(tmp_path: Path):
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"
    obj = {"z": 1, "a": {"y": [1, 2], "b": "é"}}
    jsonio.write_json(a, obj)
    jsonio.write_json(b, jsonio.read_json(a))
    assert a.read_bytes() == b.read_bytes()


def test_layout_sorted_keys_two_space_indent_trailing_newline(tmp_path: Path):
    p = tmp_path / "x.json"
    jsonio.write_json(p, {"b": 1, "a": {"d": 2, "c": 3}})
    text = p.read_text(encoding="utf-8")
    assert text == '{\n  "a": {\n    "c": 3,\n    "d": 2\n  },\n  "b": 1\n}\n'
    assert b"\r" not in p.read_bytes()


def test_non_ascii_is_kept_as_is(tmp_path: Path):
    p = tmp_path / "x.json"
    jsonio.write_json(p, {"name": "é"})
    assert "é" in p.read_text(encoding="utf-8")
    assert "\\u00e9" not in p.read_text(encoding="utf-8")


def test_no_tmp_file_left_behind(tmp_path: Path):
    p = tmp_path / "x.json"
    jsonio.write_json(p, {"a": 1})
    assert sorted(q.name for q in tmp_path.iterdir()) == ["x.json"]


def test_a_stale_tmp_file_does_not_touch_the_target(tmp_path: Path):
    # Simulate a crash after the tmp file was written but before the replace.
    p = tmp_path / "x.json"
    jsonio.write_json(p, {"a": 1})
    before = p.read_bytes()
    p.with_suffix(".json.tmp").write_text("half-written", encoding="utf-8")
    assert p.read_bytes() == before
    assert jsonio.read_json(p) == {"a": 1}


def test_read_json_raises_on_bad_file(tmp_path: Path):
    p = tmp_path / "bad.json"
    p.write_text("{not json", encoding="utf-8")
    try:
        jsonio.read_json(p)
    except json.JSONDecodeError:
        return
    raise AssertionError("read_json swallowed the error")
