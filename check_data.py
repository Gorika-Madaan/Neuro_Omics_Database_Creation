import sqlite3

def check_data():
    try:
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        cursor.execute('SELECT * FROM genomics_data WHERE disease = "Parkinson\'s Disease"')
        rows = cursor.fetchall()

        print("Data in genomics_data for Parkinson's Disease:")
        for row in rows:
            print(row)

    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
    finally:
        if connection:
            connection.close()

check_data()
