import sqlite3

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

cursor.executemany('''
INSERT INTO genomics_data (disease, gene, description) VALUES (?, ?, ?)
''', sample_data)

# Commit the changes and close the connection
connection.commit()
connection.close()

print("Database initialized with sample data.")
