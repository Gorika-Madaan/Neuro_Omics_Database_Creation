import time
import requests
import sqlite3
from flask import Flask, render_template

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

def fetch_gene_ids(disease_keyword):
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "gene",
        "term": disease_keyword,
        "retmode": "json",
        "retmax": 50  # Adjust as needed
    }
    response = requests.get(search_url, params=params)
    response.raise_for_status()
    search_results = response.json()
    gene_ids = search_results['esearchresult']['idlist']
    return gene_ids

def fetch_gene_details(gene_id):
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }
    retry_attempts = 5
    delay = 1
    for attempt in range(retry_attempts):
        try:
            response = requests.get(summary_url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                print(f"Rate limit exceeded. Retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
            else:
                raise e
    raise Exception("Max retry attempts exceeded") 
    
def store_gene_data(connection, disease, gene_id, gene_data):
    if 'result' in gene_data and gene_id in gene_data['result']:
        gene_summary = gene_data['result'][gene_id]

        gene_name = str(gene_summary.get('name', 'Unknown'))
        description = str(gene_summary.get('description', 'No description available'))
        location = str(gene_summary.get('maploc', 'Unknown location'))
        mim = str(gene_summary.get('mim', 'N/A'))
        ncbi_link = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"

        try:
            with connection:
                connection.execute(
                    """
                    INSERT INTO genomics_data (disease, gene_id, gene_name, location, mim, ncbi_link, description)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (disease, gene_id, gene_name, location, mim, ncbi_link, description)
                )
        except sqlite3.Error as e:
            print(f"SQLite error: {e}")

@app.route('/parkinsons_genomics')
def parkinsons_genomics():
    disease_keyword = "Parkinson's Disease"
    conn = sqlite3.connect('database.db')

    # Fetch gene IDs related to the disease keyword
    gene_ids = fetch_gene_ids(disease_keyword)
    print(f"Fetched {len(gene_ids)} gene IDs for {disease_keyword}.")

    # Fetch details for each gene ID and store them in the database
    for gene_id in gene_ids:
        gene_data = fetch_gene_details(gene_id)
        if gene_data:
            store_gene_data(conn, disease_keyword, gene_id, gene_data)

    # Retrieve and display the stored data
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM genomics_data WHERE disease = ?', (disease_keyword,))
    rows = cursor.fetchall()

    conn.close()

    # Render the data in an HTML template
    return render_template('parkinsons_genomics.html', genes=rows)

if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)
