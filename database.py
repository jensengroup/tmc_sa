import sqlite3
import logging

logger = logging.getLogger(__name__)


class DatabaseManager:
    def __init__(self, db_path):
        """Initialize the database manager."""
        self.db_path = db_path
        self.initialize_database()

    def initialize_database(self):
        """Initialize the SQLite database and create a table if it does not exist."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        try:
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS familiarity_scores (
                    smiles TEXT PRIMARY KEY,
                    familiarity1 REAL,
                    familiarity2 REAL
                )
                """
            )
            conn.commit()
            logger.info(f"Database initialized at {self.db_path}")
        except sqlite3.Error as e:
            logger.error(f"Error initializing database: {e}")
        finally:
            conn.close()

    def get_existing_familiarity_scores(self, smiles):
        """Retrieve the familiarity scores for a SMILES string if they exist."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        try:
            cursor.execute(
                "SELECT familiarity1, familiarity2 FROM familiarity_scores WHERE smiles = ?",
                (smiles,),
            )
            row = cursor.fetchone()
            if row:
                return row
            return None
        except sqlite3.Error as e:
            logger.error(f"Failed to retrieve familiarity scores: {e}")
            return None
        finally:
            conn.close()

    def store_familiarity_scores(self, smiles, familiarity1, familiarity2):
        """Store the familiarity scores for a SMILES string in the database."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        try:
            cursor.execute(
                """
                INSERT INTO familiarity_scores (smiles, familiarity1, familiarity2)
                VALUES (?, ?, ?)
                ON CONFLICT(smiles) DO UPDATE SET 
                    familiarity1 = excluded.familiarity1,
                    familiarity2 = excluded.familiarity2
                """,
                (smiles, familiarity1, familiarity2),
            )
            conn.commit()
            logger.info(f"Stored familiarity scores for SMILES: {smiles}")
        except sqlite3.Error as e:
            logger.error(f"Failed to store familiarity scores: {e}")
        finally:
            conn.close()
