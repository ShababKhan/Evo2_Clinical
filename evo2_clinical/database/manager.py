"""
Database management module for the Evo2_Clinical package.

This module provides classes and functions for creating, populating, and querying
databases that store genetic variant information, lncRNA data, and analysis results.
"""

import os
import sqlite3
import pandas as pd
import logging
from typing import Dict, List, Union, Optional, Tuple, Any

from ..config import config

logger = logging.getLogger(__name__)


class DatabaseManager:
    """
    Manages the creation and population of the project databases.
    """
    
    def __init__(self, db_dir: Optional[str] = None):
        """
        Initializes the DatabaseManager.

        Args:
            db_dir: Directory where databases will be stored. If None, uses the
                   directory specified in the config.
        """
        self.db_dir = db_dir or config.OUTPUT_DIRS.get('databases', 'output/databases')
        
        # Convert to absolute path if needed
        if not os.path.isabs(self.db_dir):
            # Get the package root directory (assuming config is in the package root)
            package_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            self.db_dir = os.path.join(package_root, self.db_dir)
        
        # Ensure database directory exists
        os.makedirs(self.db_dir, exist_ok=True)
        logger.info(f"Database directory set to {self.db_dir}")
        
        # Connection cache to avoid reconnecting repeatedly
        self._connections = {}
    
    def create_database(self, db_name: str, schema: Dict[str, List[str]]) -> bool:
        """
        Creates a new SQLite database with the specified schema.

        Args:
            db_name: Name of the database (without .db extension).
            schema: Dictionary defining the database schema.
                   Keys are table names, values are lists of column definitions.

        Returns:
            bool: True if successful, False otherwise.
            
        Example:
            schema = {
                'variants': [
                    'id INTEGER PRIMARY KEY',
                    'chrom TEXT',
                    'pos INTEGER',
                    'ref TEXT',
                    'alt TEXT',
                    'functional_score REAL'
                ]
            }
            create_database('my_database', schema)
        """
        db_path = os.path.join(self.db_dir, f"{db_name}.db")
        logger.info(f"Creating database: {db_path}")
        
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create each table defined in the schema
            for table_name, columns in schema.items():
                create_table_sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join(columns)})"
                logger.debug(f"Executing SQL: {create_table_sql}")
                cursor.execute(create_table_sql)
            
            conn.commit()
            conn.close()
            logger.info(f"Successfully created database {db_name} with tables: {', '.join(schema.keys())}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating database {db_name}: {e}")
            if os.path.exists(db_path):
                logger.warning(f"Removing partially created database file: {db_path}")
                try:
                    os.remove(db_path)
                except:
                    pass
            return False
    
    def populate_database(self, db_name: str, data: pd.DataFrame, table_name: Optional[str] = None) -> bool:
        """
        Populates a database table with data from a DataFrame.

        Args:
            db_name: Name of the database (without .db extension).
            data: DataFrame containing the data to insert.
            table_name: Name of the table to populate. Defaults to db_name if None.

        Returns:
            bool: True if successful, False otherwise.
        """
        table_name = table_name or db_name
        db_path = os.path.join(self.db_dir, f"{db_name}.db")
        logger.info(f"Populating table '{table_name}' in database '{db_name}' with {len(data)} rows")
        
        if len(data) == 0:
            logger.warning("No data to insert, returning")
            return True
        
        try:
            # Get database connection
            conn = self._get_connection(db_name)
            
            # Insert data
            data.to_sql(table_name, conn, if_exists='append', index=False)
            conn.commit()
            
            logger.info(f"Successfully inserted {len(data)} rows into {db_name}.{table_name}")
            return True
            
        except Exception as e:
            logger.error(f"Error populating database {db_name} table {table_name}: {e}")
            return False
    
    def query_database(self, db_name: str, query: str, params: Optional[Tuple] = None) -> pd.DataFrame:
        """
        Executes a query on a database and returns the results as a DataFrame.

        Args:
            db_name: Name of the database (without .db extension).
            query: SQL query string.
            params: Optional tuple of parameters to bind to the query.

        Returns:
            DataFrame containing the query results.
            
        Raises:
            sqlite3.Error: If there's an error executing the query.
            FileNotFoundError: If the database doesn't exist.
        """
        db_path = os.path.join(self.db_dir, f"{db_name}.db")
        
        if not os.path.exists(db_path):
            logger.error(f"Database file not found: {db_path}")
            raise FileNotFoundError(f"Database file not found: {db_path}")
        
        logger.info(f"Executing query on database '{db_name}': {query}")
        
        try:
            # Get database connection
            conn = self._get_connection(db_name)
            
            # Execute query
            if params:
                results_df = pd.read_sql_query(query, conn, params=params)
            else:
                results_df = pd.read_sql_query(query, conn)
            
            logger.info(f"Query returned {len(results_df)} rows")
            return results_df
            
        except sqlite3.Error as e:
            logger.error(f"SQL error executing query on {db_name}: {e}")
            raise
        except Exception as e:
            logger.error(f"Error querying database {db_name}: {e}")
            raise
    
    def get_table_info(self, db_name: str, table_name: str) -> pd.DataFrame:
        """
        Gets information about a table's schema.

        Args:
            db_name: Name of the database (without .db extension).
            table_name: Name of the table.

        Returns:
            DataFrame containing table schema information.
        """
        query = f"PRAGMA table_info({table_name})"
        return self.query_database(db_name, query)
    
    def list_tables(self, db_name: str) -> List[str]:
        """
        Lists all tables in a database.

        Args:
            db_name: Name of the database (without .db extension).

        Returns:
            List of table names.
        """
        query = "SELECT name FROM sqlite_master WHERE type='table'"
        tables_df = self.query_database(db_name, query)
        return tables_df['name'].tolist() if 'name' in tables_df.columns else []
    
    def close_all_connections(self) -> None:
        """Closes all database connections."""
        for db_name, conn in self._connections.items():
            try:
                conn.close()
                logger.debug(f"Closed connection to database {db_name}")
            except:
                pass
        self._connections = {}
    
    def _get_connection(self, db_name: str) -> sqlite3.Connection:
        """Gets a cached connection or creates a new one."""
        db_path = os.path.join(self.db_dir, f"{db_name}.db")
        
        # Use cached connection if available
        if db_name in self._connections:
            conn = self._connections[db_name]
            try:
                # Test if connection is still good
                conn.execute("SELECT 1")
                return conn
            except:
                # Connection is dead, remove from cache
                del self._connections[db_name]
        
        # Create new connection
        conn = sqlite3.connect(db_path)
        self._connections[db_name] = conn
        return conn


def create_endothelial_variants_database(db_manager: DatabaseManager) -> bool:
    """
    Creates the database for endothelial gene variants.

    Args:
        db_manager: Instance of DatabaseManager.

    Returns:
        bool: True if successful, False otherwise.
    """
    schema = {
        'variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'gene TEXT',
            'is_common INTEGER',  # 0 = false, 1 = true
            'is_rare INTEGER',    # 0 = false, 1 = true
            'allele_freq REAL',
            'functional_score REAL'
        ]
    }
    
    return db_manager.create_database('endothelial_variants_db', schema)


def create_lncrna_functionality_database(db_manager: DatabaseManager) -> bool:
    """
    Creates the database for lncRNA functionality.

    Args:
        db_manager: Instance of DatabaseManager.

    Returns:
        bool: True if successful, False otherwise.
    """
    schema = {
        'lncrnas': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'lncrna_name TEXT UNIQUE',
            'functionality_score REAL',
            'structure_score REAL',
            'conservation_score REAL',
            'binding_domains INTEGER'
        ]
    }
    
    return db_manager.create_database('lncrna_functionality_db', schema)


def create_gata2_as1_predictions_database(db_manager: DatabaseManager) -> bool:
    """
    Creates the database for GATA2-AS1 variant effect predictions.

    Args:
        db_manager: Instance of DatabaseManager.

    Returns:
        bool: True if successful, False otherwise.
    """
    schema = {
        'predictions': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'predicted_effect TEXT',
            'confidence_score REAL'
        ]
    }
    
    return db_manager.create_database('gata2_as1_predictions_db', schema)


def create_epigenetic_mediator_scores_database(db_manager: DatabaseManager) -> bool:
    """
    Creates the database for epigenetic mediator variant scores.

    Args:
        db_manager: Instance of DatabaseManager.

    Returns:
        bool: True if successful, False otherwise.
    """
    schema = {
        'epigenetic_variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'gene TEXT',
            'score REAL'
        ]
    }
    
    return db_manager.create_database('epigenetic_mediator_scores_db', schema)


def create_functional_variants_database(db_manager: DatabaseManager) -> bool:
    """
    Creates the database for functional variants in various pathways/diseases.

    Args:
        db_manager: Instance of DatabaseManager.

    Returns:
        bool: True if successful, False otherwise.
    """
    schema = {
        'emt_variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'gene TEXT',
            'score REAL'
        ],
        'cteph_variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'locus TEXT',
            'score REAL'
        ],
        'pah_variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'gene TEXT',
            'score REAL'
        ],
        'mesothelioma_variants': [
            'id INTEGER PRIMARY KEY AUTOINCREMENT',
            'variant_id TEXT',
            'chrom TEXT',
            'pos INTEGER',
            'ref TEXT',
            'alt TEXT',
            'gene TEXT',
            'score REAL'
        ]
    }
    
    return db_manager.create_database('functional_variants_db', schema)


def create_all_databases(db_manager: Optional[DatabaseManager] = None) -> bool:
    """
    Creates all databases needed for the project.

    Args:
        db_manager: Optional instance of DatabaseManager. If None, a new one will be created.

    Returns:
        bool: True if all databases were created successfully, False otherwise.
    """
    if db_manager is None:
        db_manager = DatabaseManager()
    
    success = True
    
    # Create each database
    database_creators = [
        create_endothelial_variants_database,
        create_lncrna_functionality_database,
        create_gata2_as1_predictions_database,
        create_epigenetic_mediator_scores_database,
        create_functional_variants_database
    ]
    
    for create_func in database_creators:
        if not create_func(db_manager):
            logger.error(f"Failed to create database using {create_func.__name__}")
            success = False
    
    return success