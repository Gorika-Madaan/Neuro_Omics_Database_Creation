from Bio import Entrez
import sqlite3

# Initialize Entrez email (required for NCBI API access)
Entrez.email = "allaboutbioinformatics21@gmail.com"  # Replace with your email address

def fetch_gene_data(disease_name, gene_symbol):
    """Fetch gene data from NCBI for a specific disease and gene symbol."""
    try:
        # Search for gene information in NCBI
        handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[sym] AND {disease_name}")
        record = Entrez.read(handle)
        handle.close()

        # If we have results, fetch the summary
        if record["IdList"]:
            gene_id = record["IdList"][0]
            handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(handle)
            handle.close()

            gene_summary = summary[0]["Summary"]
            return gene_summary
        else:
            return "No data found for the specified gene and disease."
    except Exception as e:
        print(f"Error fetching data: {e}")
        return None

def insert_gene_data():
    """Fetch gene data from NCBI and insert it into the SQLite database."""
    # Connect to the SQLite database
    connection = sqlite3.connect('database.db')
    cursor = connection.cursor()

    # List of genes and diseases to fetch data for
    genes_and_diseases = [
        ('Parkinson\'s Disease', 'LRRK2'),
        ('Alzheimer\'s Disease', 'APP')
    ]

    # Fetch and insert data into the database
    for disease, gene in genes_and_diseases:
        description = fetch_gene_data(disease, gene)
        if description:
            cursor.execute('''
            INSERT INTO genomics_data (disease, gene, description)
            VALUES (?, ?, ?)
            ''', (disease, gene, description))

    # Commit changes and close the connection
    connection.commit()
    connection.close()

    print("Gene data inserted successfully.")

if __name__ == "__main__":
    insert_gene_data()
