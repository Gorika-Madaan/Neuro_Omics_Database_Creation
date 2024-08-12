import time
import requests
import sqlite3
from flask import Flask, render_template, g

app = Flask(__name__)


@app.route("/")
def home():
    return render_template('home.html')


@app.route("/genomics")
def genomics():
    return render_template('genomics.html')


@app.route("/proteomics")
def proteomics():
    return render_template('proteomics.html')


@app.route("/transcriptomics")
def transcriptomics():
    return render_template('transcriptomics.html')


@app.route("/metabolomics")
def metabolomics():
    return render_template('metabolomics.html')

DATABASES = {
    'Genomics_db' : 'database_2.db'  # Ensure the database path is correct
}

def get_db(db_name):
    db_key = f'_database_{db_name}'
    db = getattr(g, db_key, None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASES[db_name])
        g.__setattr__(db_key, db)
    return db

@app.teardown_appcontext
def close_connection(exception):
    for db_key in g.__dict__.keys():
        if db_key.startswith('_database_'):
            db = getattr(g, db_key, None)
            if db is not None:
                db.close()

@app.route('/ALZ_genomics')
def ALZ_genomics():
    db = get_db('Genomics_db')
    cursor = db.cursor()
    cursor.execute("SELECT gene_symbol, gene_name, location, description FROM ALZHEIMER")
    genes = cursor.fetchall()
    return render_template('ALZ_genomics.html', genes=genes)




if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)
