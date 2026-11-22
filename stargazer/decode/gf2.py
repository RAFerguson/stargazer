"""Linear algebra over GF(2) for the blind code recogniser.

Representation
--------------
Matrices are stored *bit-packed*: row i of an N x n 0/1 matrix becomes ``ceil(n / 64)`` numpy ``uint64`` words, column
j living in bit ``j % 64`` of word ``j // 64`` (least significant bit first). Gaussian elimination then works one
column at a time and, for that column, XORs the pivot row into *every* other row that has the bit set with a single
vectorised numpy operation. This is the right shape for the recogniser, whose matrices are tall (thousands of blocks)
and narrow (n <= a few hundred columns): the cost is n vectorised passes over N x W words rather than N x n Python
operations. Python integers were rejected for matrices because a row operation applied to many rows at once would
still need a Python loop over rows.

Berlekamp-Massey is the opposite case: a single long vector whose update is "XOR a copy shifted by an arbitrary
number of bits". Python's arbitrary-precision integers do exactly that natively (``C ^= B << shift`` and
``(C & window).bit_count()`` run in C over 30-bit limbs), so the algorithm keeps the connection polynomials as Python
integers. Measured: 10^5 random bits in about 0.2 s, 10^4 bits in about 4 ms, 256 bits in about 50 us.

Batched elimination
-------------------
``batched_null_space`` runs elimination on many small matrices at once (shape B x m single-word rows), vectorised over
the batch. The block-code scan uses it to process every alignment offset and every random subset of one block length
in a handful of numpy passes.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "pack_rows", "unpack_rows", "rank", "row_reduce", "null_space", "berlekamp_massey", "berlekamp_massey_int",
    "linear_complexity_profile", "batched_null_space", "parity64", "pack_vector_int", "int_to_bits",
]

_WORD = 64


def _as_bits_2d(bits_2d: np.ndarray) -> np.ndarray:
    m = np.asarray(bits_2d)
    if m.ndim != 2:
        raise ValueError("expected a two-dimensional 0/1 array")
    return (m & 1).astype(np.uint8) if m.dtype != bool else m.astype(np.uint8)


def pack_rows(bits_2d: np.ndarray) -> np.ndarray:
    """Pack an N x n 0/1 matrix into an N x ceil(n/64) uint64 array (column j -> bit j % 64 of word j // 64)."""
    m = _as_bits_2d(bits_2d)
    n_rows, n = m.shape
    words = max(1, -(-n // _WORD))
    padded = np.zeros((n_rows, words * _WORD), dtype=np.uint8)
    padded[:, :n] = m
    packed_bytes = np.packbits(padded, axis=1, bitorder="little")
    return np.ascontiguousarray(packed_bytes).view("<u8").astype(np.uint64).reshape(n_rows, words)


def unpack_rows(packed: np.ndarray, n: int) -> np.ndarray:
    """Inverse of :func:`pack_rows`: an N x n uint8 0/1 matrix from packed uint64 rows."""
    p = np.ascontiguousarray(np.asarray(packed, dtype=np.uint64))
    if p.ndim == 1:
        p = p[:, None]
    as_bytes = p.astype("<u8").view(np.uint8).reshape(p.shape[0], -1)
    return np.unpackbits(as_bytes, axis=1, bitorder="little")[:, :n].astype(np.uint8)


def parity64(words: np.ndarray) -> np.ndarray:
    """Parity (0 or 1) of the number of set bits in each uint64 word."""
    return (np.bitwise_count(np.asarray(words, dtype=np.uint64)) & 1).astype(np.uint8)


def _eliminate(packed: np.ndarray, n: int, full: bool) -> tuple[np.ndarray, list[int]]:
    """Gauss(-Jordan) elimination in place on packed rows. Returns (rows, pivot columns); pivot rows come first."""
    p = packed.copy()
    n_rows = p.shape[0]
    pivots: list[int] = []
    r = 0
    for c in range(n):
        if r >= n_rows:
            break
        w, b = divmod(c, _WORD)
        shift = np.uint64(b)
        column = (p[r:, w] >> shift) & np.uint64(1)
        nz = np.flatnonzero(column)
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            p[[r, piv]] = p[[piv, r]]
        start = 0 if full else r + 1
        hits = np.flatnonzero((p[start:, w] >> shift) & np.uint64(1)) + start
        hits = hits[hits != r]
        if hits.size:
            p[hits] ^= p[r]
        pivots.append(c)
        r += 1
    return p, pivots


def rank(bits_2d: np.ndarray) -> int:
    """Rank over GF(2) by forward elimination on packed rows."""
    m = _as_bits_2d(bits_2d)
    if m.size == 0:
        return 0
    _, pivots = _eliminate(pack_rows(m), m.shape[1], full=False)
    return len(pivots)


def row_reduce(bits_2d: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Reduced row-echelon form over GF(2) and the list of pivot columns.

    The returned matrix has the same shape as the input; its first ``len(pivots)`` rows are the non-zero rows of the
    RREF (row i has its leading one in column ``pivots[i]`` and zeros in every other pivot column), the rest are zero.
    """
    m = _as_bits_2d(bits_2d)
    n_rows, n = m.shape
    if m.size == 0:
        return m.copy(), []
    reduced, pivots = _eliminate(pack_rows(m), n, full=True)
    return unpack_rows(reduced, n), pivots


def null_space(bits_2d: np.ndarray) -> np.ndarray:
    """A basis (as rows) of {h : M h = 0} over GF(2).

    For an N x n matrix whose rows are codewords, these are the parity checks of the code. The basis is the standard
    one read off the RREF: one vector per free column f, with a 1 at f and, for each pivot row i, the bit
    ``R[i, f]`` at column ``pivots[i]``. Returns an array of shape (n - rank, n).
    """
    m = _as_bits_2d(bits_2d)
    n = m.shape[1]
    if m.shape[0] == 0:
        return np.eye(n, dtype=np.uint8)
    reduced, pivots = row_reduce(m)
    free = [c for c in range(n) if c not in set(pivots)]
    basis = np.zeros((len(free), n), dtype=np.uint8)
    for i, f in enumerate(free):
        basis[i, f] = 1
        for row, pc in enumerate(pivots):
            basis[i, pc] = reduced[row, f]
    return basis


def batched_null_space(rows: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Null spaces of many small matrices at once; every row must fit one word (n <= 64).

    rows has shape (B, m): B matrices of m packed rows each. Gauss-Jordan elimination is vectorised across the batch:
    for each column, every matrix picks its first unused row with that bit set as pivot and XORs it into its other
    rows. Returns ``(basis, valid)`` where basis has shape (B, n) holding, for every column f, the free-variable null
    vector of column f packed into a uint64 word, and valid (B, n) marks the columns that are free (only those
    entries are null-space vectors). The null-space dimension of matrix b is ``valid[b].sum()``.
    """
    if n > _WORD:
        raise ValueError("batched_null_space handles rows of at most 64 bits")
    a = np.array(rows, dtype=np.uint64, copy=True)
    n_batch, m = a.shape
    used = np.zeros((n_batch, m), dtype=bool)
    pivot_row_of_col = np.full((n_batch, n), -1, dtype=np.int64)
    batch_index = np.arange(n_batch)
    one = np.uint64(1)
    bits = np.empty((n_batch, m), dtype=bool)
    candidates = np.empty((n_batch, m), dtype=bool)
    for c in range(n):
        np.not_equal(a & np.uint64(1 << c), 0, out=bits)
        np.greater(bits, used, out=candidates)          # bits & ~used
        p = np.argmax(candidates, axis=1)
        has = candidates[batch_index, p]
        if not has.any():
            continue
        pivot_rows = a[batch_index, p] * has            # zero for matrices without a pivot here
        bits[batch_index, p] = False
        a ^= pivot_rows[:, None] * bits
        hb = batch_index[has]
        used[hb, p[has]] = True
        pivot_row_of_col[hb, c] = p[has]
    valid = pivot_row_of_col < 0
    basis = np.zeros((n_batch, n), dtype=np.uint64)
    f_bits = np.uint64(1) << np.arange(n, dtype=np.uint64)
    basis[:] = f_bits[None, :]
    for c in range(n):
        rows_c = pivot_row_of_col[:, c]
        have = rows_c >= 0
        if not have.any():
            continue
        prow = np.where(have, a[batch_index, np.maximum(rows_c, 0)], np.uint64(0))
        contrib = ((prow[:, None] >> np.arange(n, dtype=np.uint64)[None, :]) & one) << np.uint64(c)
        basis |= contrib
    basis = np.where(valid, basis, np.uint64(0))
    return basis, valid


def pack_vector_int(bits: np.ndarray) -> int:
    """A 0/1 vector as a Python integer with bit i = bits[i]."""
    b = np.asarray(bits, dtype=np.uint8).ravel()
    if b.size == 0:
        return 0
    return int.from_bytes(np.packbits(b, bitorder="little").tobytes(), "little")


def int_to_bits(value: int, length: int) -> np.ndarray:
    """Bits 0..length-1 of a Python integer as a uint8 array."""
    raw = np.frombuffer(value.to_bytes(max(1, -(-length // 8)), "little"), dtype=np.uint8)
    return np.unpackbits(raw, bitorder="little")[:length].astype(np.uint8)


def _bm_core(bits: np.ndarray, want_profile: bool) -> tuple[int, int, np.ndarray | None]:
    """Berlekamp-Massey on Python integers. Returns (L, connection polynomial as int, profile or None)."""
    s = np.asarray(bits, dtype=np.uint8).ravel() & 1
    n = int(s.size)
    # r[t] = s[n-1-t]: the window s[i], s[i-1], ..., s[i-L] is the contiguous slice r[q .. q+L] with q = n-1-i,
    # so the discrepancy is the parity of C & (R >> q).
    reversed_int = pack_vector_int(s[::-1])
    c_poly, b_poly = 1, 1
    big_l, m = 0, -1
    profile = np.empty(n, dtype=np.int64) if want_profile else None
    for i in range(n):
        d = (c_poly & (reversed_int >> (n - 1 - i))).bit_count() & 1
        if d:
            t = c_poly
            c_poly ^= b_poly << (i - m)
            if 2 * big_l <= i:
                big_l = i + 1 - big_l
                b_poly, m = t, i
        if profile is not None:
            profile[i] = big_l
    return big_l, c_poly, profile


def berlekamp_massey(bits: np.ndarray) -> tuple[int, np.ndarray]:
    """Linear complexity and shortest LFSR of a 0/1 sequence (Massey 1969).

    Returns ``(L, c)`` with c a 0/1 array of length L + 1, c[0] = 1, such that for every L <= i < len(bits)
    ``sum_j c[j] * bits[i - j] = 0 (mod 2)``, i.e. bits[i] = c[1] bits[i-1] + ... + c[L] bits[i-L]. (c is the
    *connection* polynomial C(x) = 1 + c1 x + ... + cL x^L, the reciprocal of the characteristic polynomial; c[L] may
    be 0 for sequences that are not purely periodic.) The algorithm is O(n^2) bit operations, done with Python
    integers so that each step costs a few word-parallel operations; 10^5 bits take about 0.2 s.
    """
    big_l, c_poly, _ = _bm_core(bits, want_profile=False)
    return big_l, int_to_bits(c_poly, big_l + 1)


def berlekamp_massey_int(bits: np.ndarray) -> tuple[int, int]:
    """As :func:`berlekamp_massey`, with the connection polynomial as a Python integer (bit j = c[j]).

    Cheaper when many short sequences are compared (the integer is hashable and needs no unpacking).
    """
    big_l, c_poly, _ = _bm_core(bits, want_profile=False)
    return big_l, c_poly


def linear_complexity_profile(bits: np.ndarray) -> np.ndarray:
    """L_i, the linear complexity of bits[:i+1], for every prefix (full resolution, no stride).

    Berlekamp-Massey yields the whole profile at no extra cost, so it is computed at every position even for long
    inputs; memory is one int64 per bit. For a random sequence the profile follows i/2 within a few bits; an LFSR
    of length L gives a plateau at L from about 2L bits on.
    """
    _, _, profile = _bm_core(bits, want_profile=True)
    assert profile is not None
    return profile
