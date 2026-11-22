"""The message collection: signals designed to be decoded by a stranger who shares no conventions with the sender.

Content tiers. The rasters stop at tier 3 (LAYOUT): a decoder can recover the dimensions from the prime
factorisation of the length and confirm them by the compression gain and the vertical coherence of the reshaped
picture, but everything past that is a reading of shapes. Even the binary numbers drawn in the picture are legible
only under conventions the picture itself chooses (column orientation, bit order, what the marker row means), and a
resemblance of pixels to a figure, an arc or a planetary system is a perceptual judgement whose trials cannot be
counted; anchoring to shared knowledge with a trials correction is therefore not available, and anything more is a
tier-5 hypothesis. The counting sequences, the primes and the arithmetic primer reach tier 4 (ANCHORED): their
content is the integers and arithmetic, which any decoder shares, and a match (the successor function, the primes,
true sums and products) can be tested against a small, countable set of alternatives.
"""
from __future__ import annotations

import numpy as np

from stargazer.core import BitStream, ContentTier

from .zoo import ZooItem, expect, int_to_bits, register

# ---------------------------------------------------------------------------------------------------------------
# Drawing helpers (row, column) on a boolean canvas; all deterministic.
# ---------------------------------------------------------------------------------------------------------------


def _fill(canvas: np.ndarray, r0: int, r1: int, c0: int, c1: int) -> None:
    """Fill rows r0..r1-1 and columns c0..c1-1."""
    canvas[r0:r1, c0:c1] = 1


def _binary_column(canvas: np.ndarray, value: int, col: int, top: int, n_bits: int) -> None:
    """Draw value as a vertical binary number, most significant bit at row `top`."""
    canvas[top:top + n_bits, col] = int_to_bits(value, n_bits)


def _line(canvas: np.ndarray, r0: int, c0: int, r1: int, c1: int) -> None:
    """A straight line of pixels between two points, both ends included."""
    n = max(abs(r1 - r0), abs(c1 - c0)) + 1
    rows = np.rint(np.linspace(r0, r1, n)).astype(int)
    cols = np.rint(np.linspace(c0, c1, n)).astype(int)
    canvas[rows, cols] = 1


def _disc(canvas: np.ndarray, rc: int, cc: int, radius: float, filled: bool = True) -> None:
    rows, cols = np.indices(canvas.shape)
    d = np.hypot(rows - rc, cols - cc)
    canvas[(d <= radius) if filled else (np.abs(d - radius) <= 0.5)] = 1


def _triangle(canvas: np.ndarray, apex_r: int, base_r: int, c_left: int, c_right: int) -> None:
    """A filled isosceles triangle with its apex up at apex_r and its base on row base_r."""
    centre = (c_left + c_right) / 2
    half = (c_right - c_left) / 2
    for r in range(apex_r, base_r + 1):
        w = half * (r - apex_r) / max(base_r - apex_r, 1)
        canvas[r, int(np.ceil(centre - w)):int(np.floor(centre + w)) + 1] = 1


# ---------------------------------------------------------------------------------------------------------------
# The 23 x 73 picture in the style of the 1974 Arecibo message (not a reproduction of it).
# ---------------------------------------------------------------------------------------------------------------

RASTER_W, RASTER_H = 23, 73
ATOMIC_NUMBERS = (1, 6, 7, 8, 15)


def arecibo_style_raster() -> tuple[np.ndarray, dict]:
    """The 73-row by 23-column picture and the location of every element in it (rows and columns as half-open
    (start, stop) ranges)."""
    img = np.zeros((RASTER_H, RASTER_W), dtype=np.uint8)
    elements: dict = {}
    number_cols = [1 + 2 * (n - 1) for n in range(1, 11)]
    for n, c in zip(range(1, 11), number_cols):
        _binary_column(img, n, c, top=1, n_bits=4)
        img[5, c] = 1
    elements["numbers_1_to_10"] = {"rows": (1, 5), "cols": (1, 20), "values": list(range(1, 11)),
                                   "columns": number_cols, "n_bits": 4, "msb": "top", "marker_row": 5}
    atom_cols = [3, 7, 11, 15, 19]
    for z, c in zip(ATOMIC_NUMBERS, atom_cols):
        _binary_column(img, z, c, top=8, n_bits=4)
        img[12, c] = 1
    elements["atomic_numbers"] = {"rows": (8, 12), "cols": (3, 20), "values": list(ATOMIC_NUMBERS),
                                  "columns": atom_cols, "n_bits": 4, "msb": "top", "marker_row": 12}
    # Stick figure centred on column 11.
    _fill(img, 16, 19, 10, 13)                       # head
    _fill(img, 19, 28, 11, 12)                       # body
    _fill(img, 21, 22, 7, 16)                        # arms
    _line(img, 28, 10, 32, 6)                        # legs
    _line(img, 28, 12, 32, 16)
    elements["stick_figure"] = {"rows": (16, 33), "cols": (6, 17)}
    # Row of eleven dots with the third one raised.
    dot_cols = list(range(1, 22, 2))
    raised = 2
    for i, c in enumerate(dot_cols):
        img[38 if i == raised else 40, c] = 1
    elements["dots"] = {"rows": (38, 41), "cols": (1, 22), "columns": dot_cols, "row": 40, "raised_index": raised,
                        "raised_row": 38}
    # A dish-like arc opening upwards: vertex at row 60, column 11; ends at row 46.
    prev = None
    for c in range(2, 21):
        r = 60 - int(round((c - 11) ** 2 * 14 / 81))
        img[r, c] = 1
        if prev is not None and abs(prev - r) > 1:      # keep the arc 8-connected on its steep flanks
            lo, hi = sorted((prev, r))
            img[lo + 1:hi, c] = 1
        prev = r
    elements["arc"] = {"rows": (46, 61), "cols": (2, 21), "vertex": (60, 11)}
    # A filled bar at the foot of the picture.
    _fill(img, 66, 69, 4, 19)
    elements["base_bar"] = {"rows": (66, 69), "cols": (4, 19)}
    return img, elements


def _raster_truth(img: np.ndarray, elements: dict) -> dict:
    h, w = img.shape
    return {"width": w, "height": h, "image": img.copy(), "order": "row_major", "first_pixel": "top_left",
            "elements": elements, "length": img.size}


@register("message/raster_23x73")
def raster_23x73(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """An Arecibo-style 23 x 73 raster (1,679 bits, 23 x 73 both prime), sent row by row from the top left.
    truth: width, height, the 2-D image, scan order, and where each element sits."""
    img, elements = arecibo_style_raster()
    return ZooItem("message/raster_23x73", "message", BitStream(img.ravel()), _raster_truth(img, elements),
                   expect(True, ContentTier.LAYOUT, "designed"),
                   "A 23-column by 73-row picture with binary numbers, a figure, dots and an arc.")


# ---------------------------------------------------------------------------------------------------------------
# The 127 x 127 picture with a 5 x 7 dot-matrix font.
# ---------------------------------------------------------------------------------------------------------------

FONT_5X7: dict[str, tuple[str, ...]] = {
    "0": (" ### ", "#   #", "#  ##", "# # #", "##  #", "#   #", " ### "),
    "1": ("  #  ", " ##  ", "  #  ", "  #  ", "  #  ", "  #  ", " ### "),
    "2": (" ### ", "#   #", "    #", "   # ", "  #  ", " #   ", "#####"),
    "3": ("#####", "   # ", "  #  ", "   # ", "    #", "#   #", " ### "),
    "4": ("   # ", "  ## ", " # # ", "#  # ", "#####", "   # ", "   # "),
    "5": ("#####", "#    ", "#### ", "    #", "    #", "#   #", " ### "),
    "6": ("  ## ", " #   ", "#    ", "#### ", "#   #", "#   #", " ### "),
    "7": ("#####", "    #", "   # ", "  #  ", " #   ", " #   ", " #   "),
    "8": (" ### ", "#   #", "#   #", " ### ", "#   #", "#   #", " ### "),
    "9": (" ### ", "#   #", "#   #", " ####", "    #", "   # ", " ##  "),
    "+": ("     ", "  #  ", "  #  ", "#####", "  #  ", "  #  ", "     "),
    "-": ("     ", "     ", "     ", "#####", "     ", "     ", "     "),
    "=": ("     ", "     ", "#####", "     ", "#####", "     ", "     "),
    "x": ("     ", "#   #", " # # ", "  #  ", " # # ", "#   #", "     "),
    ".": ("     ", "     ", "     ", "     ", "     ", " ##  ", " ##  "),
    " ": ("     ",) * 7,
}
GLYPH_W, GLYPH_H, ADVANCE = 5, 7, 6


def glyph(ch: str) -> np.ndarray:
    """The 7 x 5 bitmap of one character of FONT_5X7."""
    return np.array([[1 if p == "#" else 0 for p in row] for row in FONT_5X7[ch]], dtype=np.uint8)


def draw_text(canvas: np.ndarray, text: str, row: int, col: int) -> None:
    """Draw text in FONT_5X7 with its top-left corner at (row, col), one blank column between characters."""
    for i, ch in enumerate(text):
        c = col + i * ADVANCE
        canvas[row:row + GLYPH_H, c:c + GLYPH_W] |= glyph(ch)


TEXT_LINES_127 = (("0123456789", 4, 4), ("1+1=2", 14, 4), ("2+3=5", 24, 4), ("3x4=12", 34, 4), ("7-2=5", 44, 4))


def raster_127_picture() -> tuple[np.ndarray, dict]:
    """The 127 x 127 picture: a border, five lines of dot-matrix text, a separator line and four shapes."""
    n = 127
    img = np.zeros((n, n), dtype=np.uint8)
    img[0, :] = img[-1, :] = img[:, 0] = img[:, -1] = 1
    elements: dict = {"border": {"rows": (0, n), "cols": (0, n), "kind": "outline"}}
    lines = []
    for text, r, c in TEXT_LINES_127:
        draw_text(img, text, r, c)
        lines.append({"text": text, "row": r, "col": c, "rows": (r, r + GLYPH_H),
                      "cols": (c, c + ADVANCE * len(text) - 1)})
    elements["text"] = {"lines": lines, "glyph_width": GLYPH_W, "glyph_height": GLYPH_H, "advance": ADVANCE,
                        "font": {k: glyph(k) for k in FONT_5X7}}
    _disc(img, 30, 100, 14)
    elements["disc"] = {"centre": (30, 100), "radius": 14, "rows": (16, 45), "cols": (86, 115)}
    _fill(img, 60, 61, 4, 123)
    elements["separator"] = {"rows": (60, 61), "cols": (4, 123)}
    _fill(img, 70, 96, 10, 36)
    elements["square"] = {"rows": (70, 96), "cols": (10, 36)}
    _triangle(img, 70, 110, 46, 86)
    elements["triangle"] = {"rows": (70, 111), "cols": (46, 87), "apex": (70, 66)}
    _disc(img, 92, 107, 14, filled=False)
    elements["circle"] = {"centre": (92, 107), "radius": 14, "rows": (78, 107), "cols": (93, 122)}
    return img, elements


@register("message/raster_127x127")
def raster_127x127(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A 127 x 127 raster (both prime) with dot-matrix digits and arithmetic, and geometric shapes, sent row by
    row from the top left. truth: width, height, image, the font, and where each element sits."""
    img, elements = raster_127_picture()
    return ZooItem("message/raster_127x127", "message", BitStream(img.ravel()), _raster_truth(img, elements),
                   expect(True, ContentTier.LAYOUT, "designed"),
                   "A 127 x 127 picture with dot-matrix digits, simple sums and geometric shapes.")


@register("message/raster_23x73_framed")
def raster_23x73_framed(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """The 23 x 73 raster preceded by a 40-bit random header and followed by a 24-bit footer of alternating bits
    (1010...). truth: header and footer lengths and bits, the raster's offset, and the raster truth."""
    img, elements = arecibo_style_raster()
    header = (rng.random(40) < 0.5).astype(np.uint8)
    footer = np.tile(np.array([1, 0], dtype=np.uint8), 12)
    bits = np.concatenate([header, img.ravel(), footer])
    truth = _raster_truth(img, elements)
    truth.update({"header_len": 40, "footer_len": 24, "header": header, "footer": footer, "offset": 40,
                  "length": bits.size})
    return ZooItem("message/raster_23x73_framed", "message", BitStream(bits), truth,
                   expect(True, ContentTier.LAYOUT, "designed"),
                   "The 23 x 73 picture between a 40-bit header and a 24-bit alternating footer.")


@register("message/raster3d_11x13x17")
def raster3d(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A three-dimensional raster of 11 x 13 x 17 voxels (three primes): the twelve edges of the bounding box and a
    solid ellipsoid inside it. Sent in C order: index = (i * 13 + j) * 17 + k, the last axis fastest. truth: dims,
    order, the 3-D array, the ellipsoid centre and semi-axes."""
    dims = (11, 13, 17)
    i, j, k = np.indices(dims)
    at_edge = ((i == 0) | (i == dims[0] - 1)).astype(int) + ((j == 0) | (j == dims[1] - 1)) + \
        ((k == 0) | (k == dims[2] - 1))
    box = at_edge >= 2
    centre, axes = (5.0, 6.0, 8.0), (3.5, 4.5, 6.5)
    ellipsoid = ((i - centre[0]) / axes[0]) ** 2 + ((j - centre[1]) / axes[1]) ** 2 + \
        ((k - centre[2]) / axes[2]) ** 2 <= 1.0
    vol = (box | ellipsoid).astype(np.uint8)
    truth = {"dims": dims, "order": "C", "fastest_axis": 2, "volume": vol, "box_edges": box, "ellipsoid": ellipsoid,
             "ellipsoid_centre": centre, "ellipsoid_semi_axes": axes, "length": vol.size}
    return ZooItem("message/raster3d_11x13x17", "message", BitStream(vol.ravel()), truth,
                   expect(True, ContentTier.LAYOUT, "designed"),
                   "An 11 x 13 x 17 voxel raster: a box outline around a solid ellipsoid.")


# ---------------------------------------------------------------------------------------------------------------
# Arithmetic anchors.
# ---------------------------------------------------------------------------------------------------------------


def unary(values: list[int]) -> np.ndarray:
    """Each value as a run of that many ones, runs separated by single zeros (no trailing zero)."""
    parts = []
    for i, v in enumerate(values):
        if i:
            parts.append(np.zeros(1, np.uint8))
        parts.append(np.ones(v, np.uint8))
    return np.concatenate(parts)


def primes_below(n: int) -> list[int]:
    sieve = np.ones(n, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    return [int(p) for p in np.flatnonzero(sieve)]


_UNARY_TRUTH = {"representation": "unary", "symbol": 1, "separator": 0, "trailing_separator": False}


@register("message/unary_count_1_40")
def unary_count(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Counting in unary: runs of 1, 2, ..., 40 ones separated by single zeros. truth: the sequence and its
    representation."""
    seq = list(range(1, 41))
    bits = unary(seq)
    return ZooItem("message/unary_count_1_40", "message", BitStream(bits),
                   dict(_UNARY_TRUTH, sequence=seq, rule="n", length=bits.size),
                   expect(True, ContentTier.ANCHORED, "designed"), "Counting from 1 to 40 in unary.")


@register("message/unary_primes_300")
def unary_primes(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """The primes below 300 in unary: runs of 2, 3, 5, 7, ... 293 ones separated by single zeros. truth: the
    sequence and its representation."""
    seq = primes_below(300)
    bits = unary(seq)
    return ZooItem("message/unary_primes_300", "message", BitStream(bits),
                   dict(_UNARY_TRUTH, sequence=seq, rule="primes below 300", length=bits.size),
                   expect(True, ContentTier.ANCHORED, "designed"), "The 62 primes below 300 in unary.")


BINARY_DELIMITER = 0  # 00000000: the only 8-bit pattern that is not a number in 1..255


@register("message/binary_count_255")
def binary_count(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """Counting in binary: the integers 1 to 255 as 8-bit fields, most significant bit first, with the delimiter
    00000000 between consecutive numbers. truth: field width, bit order, delimiter, the sequence."""
    seq = list(range(1, 256))
    parts = []
    for i, v in enumerate(seq):
        if i:
            parts.append(int_to_bits(BINARY_DELIMITER, 8))
        parts.append(int_to_bits(v, 8))
    bits = np.concatenate(parts)
    truth = {"field_width": 8, "bit_order": "msb_first", "delimiter": BINARY_DELIMITER, "delimiter_width": 8,
             "delimiter_placement": "between", "sequence": seq, "length": bits.size}
    return ZooItem("message/binary_count_255", "message", BitStream(bits), truth,
                   expect(True, ContentTier.ANCHORED, "designed"),
                   "The integers 1 to 255 in 8-bit binary separated by a zero byte.")


PRIMER_CODES = {"+": 252, "x": 253, "=": 254, ";": 255}  # numbers 0..251 are their own codes


@register("message/arithmetic_primer")
def arithmetic_primer(rng: np.random.Generator, size: str = "small") -> ZooItem:
    """A Lincos-style arithmetic primer: true statements "a + b = c" and "a x b = c" with a, b in 1..12, each
    written as six 8-bit tokens [a, op, b, =, c, ;], most significant bit first. truth: token width, code book,
    the statements as (a, op, b, c)."""
    n = 400 if size == "small" else 2_000
    a = rng.integers(1, 13, n)
    b = rng.integers(1, 13, n)
    ops = np.where(rng.random(n) < 0.5, "+", "x")
    statements = [(int(x), str(op), int(y), int(x + y if op == "+" else x * y)) for x, op, y in zip(a, ops, b)]
    tokens = []
    for x, op, y, c in statements:
        tokens += [x, PRIMER_CODES[op], y, PRIMER_CODES["="], c, PRIMER_CODES[";"]]
    bits = np.concatenate([int_to_bits(t, 8) for t in tokens])
    truth = {"token_width": 8, "bit_order": "msb_first", "code_book": dict(PRIMER_CODES),
             "numbers": "binary value of the token, 0..251", "statement_layout": ("a", "op", "b", "=", "c", ";"),
             "operand_range": (1, 12), "statements": statements, "tokens": np.array(tokens), "length": bits.size}
    return ZooItem("message/arithmetic_primer", "message", BitStream(bits), truth,
                   expect(True, ContentTier.HYPOTHESIS, "designed"),
                   "A Lincos-style primer of true sums and products in 8-bit tokens.")
