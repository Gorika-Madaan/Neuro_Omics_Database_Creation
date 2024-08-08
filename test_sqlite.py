import sqlite3

def test_db():
    try:
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        # Check existing tables
        cursor.execute('SELECT name FROM sqlite_master WHERE type="table";')
        tables = cursor.fetchall()
        print("Tables in database:", tables)

        # Check schema of genomics_data
        cursor.execute('PRAGMA table_info(genomics_data);')
        schema = cursor.fetchall()
        print("Schema of genomics_data:", schema)

        # Fetch data for Parkinson's Disease
        cursor.execute('''
        SELECT gene, description FROM genomics_data WHERE disease = 'Parkinson\'s Disease'
        ''')
        data = cursor.fetchall()
        print("Data fetched:", data)

        connection.close()
    except sqlite3.OperationalError as e:
        print("SQLite OperationalError:", e)
    except sqlite3.DatabaseError as e:
        print("SQLite DatabaseError:", e)
    except Exception as e:
        print("Error:", e)

if __name__ == "__main__":
    test_db()
