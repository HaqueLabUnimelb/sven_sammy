"""Write a copy of an archived script with literal path replacements, and print the diff.

Each OLD string must occur exactly once in the source file.

Usage: python make_command_copy.py SRC DST OLD1 NEW1 [OLD2 NEW2 ...]
"""
import difflib
import sys

src, dst = sys.argv[1], sys.argv[2]
pairs = sys.argv[3:]
if not pairs or len(pairs) % 2:
    sys.exit("ERROR: give OLD NEW pairs")

with open(src, encoding="utf-8") as fh:
    text = fh.read()
orig = text
for old, new in zip(pairs[0::2], pairs[1::2]):
    n = text.count(old)
    if n != 1:
        sys.exit(f"ERROR: expected exactly one occurrence of {old!r}, found {n}")
    text = text.replace(old, new)

with open(dst, "w", encoding="utf-8") as fh:
    fh.write(text)

print("\n".join(difflib.unified_diff(orig.splitlines(), text.splitlines(),
                                     fromfile=src, tofile=dst, lineterm="", n=0)))
