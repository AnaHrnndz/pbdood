#!/usr/bin/env python
"""
check_names.py -- Validate and sanitize FASTA sequence names for downstream
phylogenetics.

The Newick tree format reserves the characters  : ; , ( ) [ ]  ' "  and
whitespace, so any of them inside a sequence ID corrupts the trees. For example,
FastTree reads the part after a ':' as a branch length, so '386585.gene:10367786'
is split into name '386585.gene' + length '10367786'.

This script rewrites each ID keeping only safe characters -- letters, digits,
'.', '_', '-' and the configured taxid delimiter -- and replacing anything else
with '_' (or '-' when the delimiter itself is '_'). The taxid delimiter is always
preserved, so taxid extraction (OG_Delineation --sp_delimitator) keeps working
whatever delimiter you use.

Usage:
    check_names.py <input.faa> <output.faa> <map.tsv> [delimiter]

    delimiter   taxid delimiter to preserve (default '.').

It exits with an error if:
  * the delimiter is itself a Newick-reserved character (it can't be both a taxid
    separator and tree-safe), or
  * two different original IDs collapse to the same sanitized ID (which would make
    those sequences indistinguishable downstream).
"""

import re
import sys

NEWICK_RESERVED = set(":;,()[]'\" \t\n\r")


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: check_names.py <input.faa> <output.faa> <map.tsv> [delimiter]")

    infile, outfile, mapfile = sys.argv[1], sys.argv[2], sys.argv[3]
    delimiter = sys.argv[4] if len(sys.argv) > 4 else "."

    # The delimiter must survive into the trees, so it cannot be a Newick metachar.
    if delimiter in NEWICK_RESERVED:
        sys.exit(
            "[check_names] ERROR: the taxid delimiter '{}' is a Newick-reserved "
            "character; it can't be both a taxid separator and tree-safe. Choose a "
            "different --ogd_sp_delimitator or reformat the IDs.".format(delimiter)
        )

    # Safe = alphanumerics + . _ - + the delimiter. Everything else -> replacement.
    unsafe = re.compile("[^A-Za-z0-9._\\-" + re.escape(delimiter) + "]")
    # Don't turn replaced chars into the delimiter itself (would confuse taxid split).
    replacement = "-" if delimiter == "_" else "_"

    n_seqs = 0
    n_changed = 0
    n_no_delim = 0
    clean2orig = {}      # sanitized -> first original seen (to detect collisions)
    collisions = []
    mapping = []         # list of (original, sanitized)

    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_seqs += 1
                # Sequence ID = first whitespace-delimited token (standard FASTA).
                orig = line[1:].strip().split()[0]
                clean = unsafe.sub(replacement, orig)
                if clean != orig:
                    n_changed += 1
                if delimiter not in clean:
                    n_no_delim += 1
                mapping.append((orig, clean))
                if clean in clean2orig and clean2orig[clean] != orig:
                    collisions.append((clean2orig[clean], orig, clean))
                clean2orig.setdefault(clean, orig)
                fout.write(">" + clean + "\n")
            else:
                fout.write(line)

    with open(mapfile, "w") as fm:
        fm.write("original\tsanitized\n")
        for orig, clean in mapping:
            fm.write("{}\t{}\n".format(orig, clean))

    sys.stderr.write(
        "[check_names] {} sequences read, {} IDs sanitized "
        "(delimiter '{}').\n".format(n_seqs, n_changed, delimiter)
    )
    if n_no_delim:
        sys.stderr.write(
            "[check_names] WARNING: {} IDs do not contain the delimiter '{}'; "
            "taxid extraction may fail for those.\n".format(n_no_delim, delimiter)
        )

    if collisions:
        sys.stderr.write(
            "[check_names] ERROR: {} name collision(s) after sanitizing "
            "(different IDs map to the same name):\n".format(len(collisions))
        )
        for a, b, c in collisions[:20]:
            sys.stderr.write("    '{}' and '{}'  ->  '{}'\n".format(a, b, c))
        if len(collisions) > 20:
            sys.stderr.write("    ... and {} more\n".format(len(collisions) - 20))
        sys.stderr.write(
            "[check_names] Make the IDs stay unique after replacing special "
            "characters, then re-run.\n"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()