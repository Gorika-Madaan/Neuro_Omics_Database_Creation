import sqlite3

def initialize_database():
    try:
        # Connect to the SQLite database (or create it if it doesn't exist)
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        print("Connected to the database.")

        # Create or alter the table for genomics data to include additional columns
        print("Creating or ensuring the 'genomics_data' table exists...")
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS genomics_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            disease TEXT NOT NULL,
            gene_id TEXT,
            gene_name TEXT,
            location TEXT,
            mim TEXT,
            ncbi_link TEXT,
            description TEXT
        )
        ''')
        print("Table creation or alteration completed.")

        # Commit the changes
        connection.commit()
        print("Changes committed to the database.")

        # Retrieve and display the schema of the genomics_data table
        print("Retrieving schema for 'genomics_data' table...")
        cursor.execute("PRAGMA table_info(genomics_data);")
        schema = cursor.fetchall()
        if schema:
            print("Schema for 'genomics_data':")
            for column in schema:
                print(f"Column Name: {column[1]}, Data Type: {column[2]}, Not Null: {column[3]}, Default Value: {column[4]}")
        else:
            print("No schema found. The table might not exist or be accessible.")

        # Retrieve and display the data from the genomics_data table
        print("Retrieving data from 'genomics_data' table...")
        cursor.execute('SELECT * FROM genomics_data')
        rows = cursor.fetchall()
        if rows:
            print("Data in 'genomics_data':")
            for row in rows:
                print(row)
        else:
            print("No data found in 'genomics_data'.")

    except sqlite3.Error as e:
        print(f"SQLite error: {e}")

    finally:
        # Close the connection
        if connection:
            connection.close()
            print("Database connection closed.")

if __name__ == "__main__":
    initialize_database()
