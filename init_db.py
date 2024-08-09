import sqlite3

try:
    # Connect to the SQLite database (or create it if it doesn't exist)
    connection = sqlite3.connect('database.db')

    # Create a cursor object to interact with the database
    cursor = connection.cursor()

    # Create a table for genomics data if it doesn't exist
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS genomics_data (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        disease TEXT NOT NULL,
        gene TEXT NOT NULL,
        description TEXT
    )
    ''')

    # Insert sample data (optional)
    sample_data = [
        ('Parkinson\'s Disease', 'Gene1', 'Description of Gene1'),
        ('Parkinson\'s Disease', 'Gene2', 'Description of Gene2')
    ]

    # Check if data already exists to avoid duplicates
    cursor.executemany('''
    INSERT INTO genomics_data (disease, gene, description)
    SELECT ?, ?, ?
    WHERE NOT EXISTS (
        SELECT 1 FROM genomics_data WHERE disease = ? AND gene = ?
    )
    ''', [(disease, gene, desc, disease, gene) for disease, gene, desc in sample_data])

    # Commit the changes
    connection.commit()

    # Retrieve and display the data
    cursor.execute('SELECT * FROM genomics_data')
    rows = cursor.fetchall()
    for row in rows:
        print(row)

except sqlite3.Error as e:
    print(f"SQLite error: {e}")

finally:
    # Close the connection
    if connection:
        connection.close()

print("Database initialized with sample data.")
