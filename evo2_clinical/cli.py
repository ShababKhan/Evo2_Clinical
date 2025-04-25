"""
Command-line interface for the Evo2_Clinical package.

This module provides command-line functionality for running the Evo2 pipeline
and related tools.
"""

import os
import sys
import argparse
import logging
from typing import List, Optional

from .config import config
from .pipeline.main import run_evo2_pipeline, Pipeline
from .utils import helpers
from .database.manager import DatabaseManager, create_all_databases

logger = logging.getLogger(__name__)


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Args:
        args: Command line arguments. If None, sys.argv[1:] will be used.
        
    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description='Evo2_Clinical: Computational framework for investigating endothelial influence on pulmonary disease'
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Run pipeline command
    run_parser = subparsers.add_parser('run', help='Run the complete Evo2 pipeline')
    run_parser.add_argument(
        '--config', '-c',
        dest='config_file',
        help='Path to configuration file'
    )
    run_parser.add_argument(
        '--log-level',
        dest='log_level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Logging level'
    )
    
    # Create databases command
    db_parser = subparsers.add_parser('create-db', help='Create and initialize databases')
    db_parser.add_argument(
        '--db-dir',
        dest='db_dir',
        help='Directory to store databases'
    )
    
    # Configure command
    config_parser = subparsers.add_parser('config', help='View or create configuration')
    config_parser.add_argument(
        'action',
        choices=['show', 'create', 'validate'],
        help='Action to perform with configuration'
    )
    config_parser.add_argument(
        '--output', '-o',
        dest='output_file',
        help='Output file for configuration'
    )
    
    # Info command
    info_parser = subparsers.add_parser('info', help='Show information about installed package')
    
    # Parse arguments
    return parser.parse_args(args)


def run_command(args: argparse.Namespace) -> None:
    """
    Run the pipeline based on command line arguments.
    
    Args:
        args: Parsed command line arguments.
    """
    # Set up logging
    log_level = getattr(logging, args.log_level) if hasattr(args, 'log_level') else logging.INFO
    log_dir = os.path.join(os.path.dirname(__file__), '..', 'logs')
    helpers.setup_logging(log_dir, log_level)
    
    try:
        logger.info(f"Running Evo2 pipeline with configuration: {args.config_file}")
        success = run_evo2_pipeline(args.config_file)
        
        if success:
            logger.info("Pipeline completed successfully!")
            sys.exit(0)
        else:
            logger.error("Pipeline failed!")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Error running pipeline: {e}", exc_info=True)
        sys.exit(1)


def create_db_command(args: argparse.Namespace) -> None:
    """
    Create databases based on command line arguments.
    
    Args:
        args: Parsed command line arguments.
    """
    # Set up logging
    log_dir = os.path.join(os.path.dirname(__file__), '..', 'logs')
    helpers.setup_logging(log_dir, logging.INFO)
    
    try:
        db_dir = args.db_dir if hasattr(args, 'db_dir') and args.db_dir else None
        db_manager = DatabaseManager(db_dir)
        
        logger.info(f"Creating databases in: {db_manager.db_dir}")
        success = create_all_databases(db_manager)
        
        if success:
            logger.info("Databases created successfully!")
            sys.exit(0)
        else:
            logger.error("Database creation failed!")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Error creating databases: {e}", exc_info=True)
        sys.exit(1)


def config_command(args: argparse.Namespace) -> None:
    """
    Handle configuration commands.
    
    Args:
        args: Parsed command line arguments.
    """
    if args.action == 'show':
        # Display current configuration
        print("Current Evo2_Clinical configuration:")
        print("\nDATA_PATHS:")
        for key, value in config.DATA_PATHS.items():
            print(f"  {key}: {value}")
        
        print("\nTOOL_PATHS:")
        for key, value in config.TOOL_PATHS.items():
            print(f"  {key}: {value}")
        
        print("\nOUTPUT_DIRS:")
        for key, value in config.OUTPUT_DIRS.items():
            print(f"  {key}: {value}")
    
    elif args.action == 'create':
        # Create a new configuration file
        output_file = args.output_file if hasattr(args, 'output_file') and args.output_file else 'evo2_clinical_config.yaml'
        
        config_dict = {
            'DATA_PATHS': config.DATA_PATHS,
            'TOOL_PATHS': config.TOOL_PATHS,
            'OUTPUT_DIRS': config.OUTPUT_DIRS
        }
        
        try:
            helpers.save_yaml_config(config_dict, output_file)
            print(f"Configuration file created: {output_file}")
        except Exception as e:
            print(f"Error creating configuration file: {e}")
            sys.exit(1)
    
    elif args.action == 'validate':
        # Validate all paths in the config
        missing_data_paths = [path for path, value in config.DATA_PATHS.items() 
                              if value and not os.path.exists(value)]
        
        missing_tool_paths = [path for path, value in config.TOOL_PATHS.items() 
                              if value and not os.path.exists(value)]
        
        if missing_data_paths:
            print("WARNING: The following data paths do not exist:")
            for path in missing_data_paths:
                print(f"  - {path}: {config.DATA_PATHS[path]}")
        
        if missing_tool_paths:
            print("WARNING: The following tool paths do not exist:")
            for path in missing_tool_paths:
                print(f"  - {path}: {config.TOOL_PATHS[path]}")
        
        if not missing_data_paths and not missing_tool_paths:
            print("Configuration is valid. All specified paths exist.")


def info_command(args: argparse.Namespace) -> None:
    """
    Display information about the installed package.
    
    Args:
        args: Parsed command line arguments.
    """
    from . import __version__
    
    print(f"Evo2_Clinical version {__version__}")
    print("\nAbout:")
    print("  A computational framework for investigating endothelial influence on")
    print("  pulmonary disease, integrating computational insights with molecular")
    print("  biology for advancements in cancer care.")
    
    print("\nInsight:")
    print("  This package implements the Evo2 Pipeline for analyzing genetic variants")
    print("  in the context of endothelial cell function, lung inflammation, and fibrosis.")
    
    print("\nAuthors:")
    print("  Jeffrey Man, Aninda Dibya Saha, Shabab Khan")
    
    print("\nContents:")
    import pkgutil
    print("  Modules:")
    for _, name, ispkg in pkgutil.iter_modules([os.path.dirname(__file__)]):
        if ispkg:
            print(f"    - {name} (package)")
        else:
            print(f"    - {name}")


def main(args: Optional[List[str]] = None) -> None:
    """
    Main entry point for the command-line interface.
    
    Args:
        args: Command line arguments. If None, sys.argv[1:] will be used.
    """
    parsed_args = parse_args(args)
    
    if parsed_args.command == 'run':
        run_command(parsed_args)
    elif parsed_args.command == 'create-db':
        create_db_command(parsed_args)
    elif parsed_args.command == 'config':
        config_command(parsed_args)
    elif parsed_args.command == 'info':
        info_command(parsed_args)
    else:
        # No command specified, show help
        parse_args(['--help'])


if __name__ == '__main__':
    main()