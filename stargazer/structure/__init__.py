"""Sequence structure: does one sequence of discrete symbols carry structure beyond what a simple null explains?

Estimators (bias-corrected entropies, block and conditional entropies, lagged mutual information, compression),
surrogates and permutation p-values, Witten-Bell n-gram models with held-out gain in bits, mutual-information decay
fits, and grammar (Re-Pair) and chunking statistics. All entropies and information values are in bits.
"""
from .estimators import (ENTROPY_METHODS, block_counts, block_entropy, compression_bits, compression_ratio,
                         conditional_entropy, entropy, entropy_chao_shen, entropy_grassberger, entropy_miller_madow,
                         entropy_plugin, entropy_rate_curve, lagged_mutual_information, max_trustworthy_order)
from .grammar import grammar_gain, kt_code_bits, repair, tp_segmentation, tp_segmentation_null
from .mi import fit_decay, mutual_information_decay
from .ngram import NGramModel, benjamini_hochberg, gain_in_bits, held_out_bits_per_symbol, transition_matrix
from .surrogates import block_shuffle, doublet_shuffle, make_surrogate, markov_surrogate, permutation_pvalue, shuffle

__all__ = [
    "ENTROPY_METHODS", "NGramModel", "benjamini_hochberg", "block_counts", "block_entropy", "block_shuffle",
    "compression_bits", "compression_ratio", "conditional_entropy", "doublet_shuffle", "entropy", "entropy_chao_shen",
    "entropy_grassberger", "entropy_miller_madow", "entropy_plugin", "entropy_rate_curve", "fit_decay",
    "gain_in_bits", "grammar_gain", "held_out_bits_per_symbol", "kt_code_bits", "lagged_mutual_information",
    "make_surrogate", "markov_surrogate", "max_trustworthy_order", "mutual_information_decay",
    "permutation_pvalue", "repair", "shuffle", "transition_matrix", "tp_segmentation", "tp_segmentation_null",
]
