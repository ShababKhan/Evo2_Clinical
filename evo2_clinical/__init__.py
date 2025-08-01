"""
Evo2_Clinical Package
=====================

A computational framework for investigating endothelial influence on pulmonary disease,
integrating computational insights with molecular biology for advancements in cancer care.

This package implements the Evo2 Pipeline for analyzing genetic variants in the context of
endothelial cell function, lung inflammation, and fibrosis.
"""

__version__ = '0.1.0'

# Set up logging for import errors
import logging
logger = logging.getLogger(__name__)

# Module imports with error handling
modules = {}

try:
    from . import config
    modules['config'] = config
except ImportError as e:
    logger.error(f"Failed to import config module: {e}")

try:
    from . import data
    modules['data'] = data
except ImportError as e:
    logger.error(f"Failed to import data module: {e}")

try:
    from . import analysis
    modules['analysis'] = analysis
except ImportError as e:
    logger.error(f"Failed to import analysis module: {e}")

try:
    from . import database
    modules['database'] = database
except ImportError as e:
    logger.error(f"Failed to import database module: {e}")

try:
    from . import pipeline
    modules['pipeline'] = pipeline
except ImportError as e:
    logger.error(f"Failed to import pipeline module: {e}")

try:
    from . import utils
    modules['utils'] = utils
except ImportError as e:
    logger.error(f"Failed to import utils module: {e}")

try:
    from . import cli
    modules['cli'] = cli
except ImportError as e:
    logger.error(f"Failed to import cli module: {e}")