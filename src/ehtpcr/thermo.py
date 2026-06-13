"""Melting temperature (Tm) models.

Default is a **verbatim port** of the legacy nearest-neighbor calculator
(original ``Primer.py`` ``TmCalculator``): Breslauer/Borer-style NN enthalpy/
entropy values with fixed 50 mM salt / 200 nM primer. It is NOT SantaLucia-1998
and must not be relabeled as such. It is the default to preserve the validated
behavior (which primers pass the Tm filter). Alternative models
(``design.tm.model``): the Biopython nearest-neighbor
tables ``breslauer``/``sugimoto``/``santalucia``/``santalucia2004``, plus
``wallace``, ``gc``, and ``primer3``.
"""
from __future__ import annotations

import math

# Nearest-neighbor enthalpy (H) and entropy (S) tables — ported verbatim from
# the legacy TmCalculator. Indexed [first base][second base].
_H = {
    "A": {"A": 9.1, "C": 6.5, "G": 7.8, "T": 8.6},
    "C": {"A": 5.8, "C": 11.0, "G": 11.9, "T": 7.8},
    "G": {"A": 5.6, "C": 11.1, "G": 11.0, "T": 6.5},
    "T": {"A": 6.0, "C": 5.6, "G": 5.8, "T": 9.1},
}
_S = {
    "A": {"A": 24.0, "C": 17.3, "G": 20.8, "T": 23.9},
    "C": {"A": 12.9, "C": 26.6, "G": 27.8, "T": 20.8},
    "G": {"A": 13.5, "C": 26.7, "G": 26.6, "T": 17.3},
    "T": {"A": 16.8, "C": 13.5, "G": 12.9, "T": 24.0},
}
_CONC_SALT = 50.0
_CONC_PRIMER = 200.0


def legacy_breslauer_tm(primer: str) -> float:
    """Verbatim port of the legacy ``TmCalculator`` (default Tm model).

    Assumes A/C/G/T (case-insensitive); other bases are skipped in the NN sum,
    matching the legacy behavior.
    """
    if not primer:
        raise ValueError("legacy_breslauer_tm: empty primer")
    p = primer.upper()
    sum_h = 0.0
    sum_s = 0.0
    bh = ""
    is_first = True
    for ch in p:
        if is_first:
            is_first = False
            bh = ch
            continue
        if ch not in _S or ch not in _H or bh not in _S or bh not in _H:
            continue
        sum_h += _H[bh][ch]
        sum_s += _S[bh][ch]
        bh = ch

    # 5' initiation
    if p[0] == "G" or p[0] == "C":
        sum_h += 0.25 / 2
        sum_s += -0.62 / 2
    else:
        sum_h += 0.7 / 2
        sum_s += 0.62 / 2
    # 3' initiation
    if p[-1] == "G" or p[-1] == "C":
        sum_h += 0.25 / 2
        sum_s += -0.62 / 2
    else:
        sum_h += 0.7 / 2
        sum_s += 0.62 / 2

    return (
        -1000 * sum_h / (-1 * (sum_s + 16.8) + 1.987 * math.log(_CONC_PRIMER / 4000000000))
        - 273.15
        + 16.6 * math.log(_CONC_SALT / 1000) / math.log(10)
    )


# Alternative Tm models. All conditions are pinned to the legacy's (Na⁺ 50 mM,
# primer 200 nM) so the numbers are comparable to the default. Nearest-neighbor
# models map to Biopython's published parameter tables:
_NN_TABLES = {
    "breslauer": "DNA_NN1",        # Breslauer et al. 1986
    "sugimoto": "DNA_NN2",         # Sugimoto et al. 1996
    "santalucia": "DNA_NN3",       # Allawi & SantaLucia 1997
    "santalucia2004": "DNA_NN4",   # SantaLucia & Hicks 2004 (unified)
}


def _biopython_mt():
    try:
        from Bio.SeqUtils import MeltingTemp as mt
    except ImportError as e:
        raise ImportError("this tm.model needs biopython (`pip install biopython`)") from e
    return mt


def tm_nn(primer: str, table: str = "DNA_NN3") -> float:
    """Nearest-neighbor Tm via a Biopython parameter table (Na⁺ 50 mM, 200 nM)."""
    if not primer:
        raise ValueError("tm_nn: empty primer")
    mt = _biopython_mt()
    # saltcorr=5 (Owczarzy 2004) pinned explicitly so the values can't silently
    # move if Biopython changes its default salt-correction method.
    return mt.Tm_NN(primer.upper(), nn_table=getattr(mt, table),
                    Na=50, dnac1=200, dnac2=200, saltcorr=5)


def tm_nn_santalucia(primer: str) -> float:
    """Allawi & SantaLucia 1997 NN Tm (Biopython ``DNA_NN3``) — a few °C cooler
    than the default legacy Breslauer on average. Needs biopython."""
    return tm_nn(primer, "DNA_NN3")


def tm_wallace(primer: str) -> float:
    """Wallace '2(A+T) + 4(G+C)' rule — a rough estimate for short oligos."""
    if not primer:
        raise ValueError("tm_wallace: empty primer")
    return float(_biopython_mt().Tm_Wallace(primer.upper()))


def tm_gc(primer: str) -> float:
    """Empirical GC-fraction Tm (Biopython ``Tm_GC``, Na⁺ 50 mM)."""
    if not primer:
        raise ValueError("tm_gc: empty primer")
    return _biopython_mt().Tm_GC(primer.upper(), Na=50)


def tm_primer3(primer: str) -> float:
    """primer3 thermodynamic Tm (SantaLucia 1998), Na⁺ 50 mM, primer 200 nM.

    The de-facto reference for PCR primer Tm. Needs primer3-py.
    """
    if not primer:
        raise ValueError("tm_primer3: empty primer")
    try:
        import primer3
    except ImportError as e:
        raise ImportError("tm.model 'primer3' needs primer3-py (`pip install primer3-py`)") from e
    # dv_conc/dntp_conc = 0 (Na⁺-only), matching the other models — primer3's
    # defaults are 1.5 mM Mg²⁺ / 0.6 mM dNTP, which would shift the Tm.
    return primer3.calc_tm(primer.upper(), mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=200)


def tm(primer: str, model: str = "legacy") -> float:
    if model == "legacy":
        return legacy_breslauer_tm(primer)
    if model in _NN_TABLES:
        return tm_nn(primer, _NN_TABLES[model])
    if model == "wallace":
        return tm_wallace(primer)
    if model == "gc":
        return tm_gc(primer)
    if model == "primer3":
        return tm_primer3(primer)
    raise NotImplementedError(f"unknown Tm model {model!r}")
