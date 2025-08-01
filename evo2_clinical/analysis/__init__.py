"""
Analysis module for the Evo2_Clinical package

This module contains functions and classes for analyzing genetic data
related to endothelial cells and lncRNAs.
"""

import logging
logger = logging.getLogger(__name__)

# Initialize components as None in case imports fail
LncRNAScorer = None
variant_scoring = None
lncrna_scoring = None
epigenetic_analysis = None

try:
    from .lncrna_scoring import LncRNAScorer
    logger.info("Successfully imported LncRNAScorer")
except ImportError as e:
    logger.error(f"Error importing LncRNAScorer: {e}")

try:
    from . import variant_scoring
    logger.info("Successfully imported variant_scoring module")
except ImportError as e:
    logger.error(f"Error importing variant_scoring module: {e}")

try:
    from . import lncrna_scoring
    logger.info("Successfully imported lncrna_scoring module")
except ImportError as e:
    logger.error(f"Error importing lncrna_scoring module: {e}")

try:
    from . import epigenetic_analysis
    logger.info("Successfully imported epigenetic_analysis module")
except ImportError as e:
    logger.error(f"Error importing epigenetic_analysis module: {e}")

# Add other analysis module imports here