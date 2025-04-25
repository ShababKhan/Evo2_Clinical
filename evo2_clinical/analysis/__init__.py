"""
Analysis modules for the Evo2_Clinical package.
"""

from .variant_scoring import (
    run_evo2_variant_scoring,
    analyze_emt_genes,
    predict_gata2_as1_variant_effects,
    analyze_epigenetic_mediators,
    run_swap_snp_pipeline,
    Evo2Runner
)

from .lncrna_scoring import (
    run_evo2_lncrna_scoring,
    analyze_specific_lncrnas,
    lncRNAScorer
)