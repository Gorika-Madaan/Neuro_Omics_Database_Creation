import pandas as pd
import sqlite3

# Load your data from Excel
data = pd.read_excel('ALZHEIMER.xlsx')

# Connect to SQLite database
conn = sqlite3.connect('database_2.db')
cursor = conn.cursor()

data.columns = data.columns.str.strip()
# Insert data into the database
data.to_sql('ALZHEIMER', conn, if_exists='append', index= False)

# Commit and close connection
conn.commit()
conn.close()
