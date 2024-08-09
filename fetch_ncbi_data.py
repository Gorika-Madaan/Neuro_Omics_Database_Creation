import requests
import sqlite3

def initialize_database():
    try:
        # Connect to the SQLite database (or create it if it doesn't exist)
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        print("Connected to the database.")

        # Create the table for genomics data if it doesn't exist
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

    except sqlite3.Error as e:
        print(f"SQLite error: {e}")

    finally:
        # Close the connection
        if connection:
            connection.close()
            print("Database connection closed.")

def fetch_gene_ids(disease_keyword):
    """
    Search NCBI for gene IDs related to the specified disease keyword.
    """
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "gene",
        "term": disease_keyword,
        "retmode": "json",
        "retmax": 100  # Set the number of results you want to fetch
    }

    response = requests.get(search_url, params=params)
    response.raise_for_status()
    search_results = response.json()
    gene_ids = search_results['esearchresult']['idlist']
    return gene_ids

def fetch_gene_details(gene_id):
    """
    Fetch detailed information about a gene from NCBI using its gene ID.
    """
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }

    response = requests.get(summary_url, params=params)
    response.raise_for_status()
    return response.json()

def store_gene_data(connection, disease, gene_id, gene_data):
    """
    Store fetched gene data into the SQLite database.
    """
    if 'result' in gene_data and gene_id in gene_data['result']:
        gene_summary = gene_data['result'][gene_id]

        # Debug output to verify the structure of gene_summary
        print(f"Gene Summary for {gene_id}: {gene_summary}")

        try:
            gene_name = str(gene_summary.get('name', 'Unknown'))
            description = str(gene_summary.get('description', 'No description available'))
            location = str(gene_summary.get('maploc', 'Unknown location'))
            mim = str(gene_summary.get('mim', 'N/A'))
            ncbi_link = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"

            # Print types to verify they are not lists
            print(f"Inserting data: disease={disease}, gene_id={gene_id}, gene_name={gene_name} (type: {type(gene_name)}), location={location} (type: {type(location)}), mim={mim} (type: {type(mim)}), ncbi_link={ncbi_link} (type: {type(ncbi_link)}), description={description} (type: {type(description)})")

            with connection:
                connection.execute(
                    """
                    INSERT INTO genomics_data (disease, gene_id, gene_name, location, mim, ncbi_link, description)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (disease, gene_id, gene_name, location, mim, ncbi_link, description)
                )
        except Exception as e:
            print(f"Error processing gene data for gene ID {gene_id}: {e}")
    else:
        print(f"Gene data not found for gene ID: {gene_id}")

def main():
    # Initialize the database
    initialize_database()

    # Connect to the SQLite database
    conn = sqlite3.connect('database.db')

    # Specify the disease keyword
    disease_keyword = "Parkinson's Disease"

    # Fetch gene IDs related to the disease keyword
    gene_ids = fetch_gene_ids(disease_keyword)
    print(f"Fetched {len(gene_ids)} gene IDs for {disease_keyword}.")

    # Fetch details for each gene ID and store them in the database
    for gene_id in gene_ids:
        gene_data = fetch_gene_details(gene_id)
        if gene_data:
            # Print the entire gene_data for debugging
            print(f"Gene Data for {gene_id}: {gene_data}")
            store_gene_data(conn, disease_keyword, gene_id, gene_data)

    # Display the stored data
    print("Displaying stored data from 'genomics_data' table:")
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM genomics_data')
    rows = cursor.fetchall()
    for row in rows:
        print(row)

    # Close the database connection
    conn.close()

if __name__ == "__main__":
    main()
