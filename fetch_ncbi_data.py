from Bio import Entrez
import sqlite3

# Initialize Entrez email (required for NCBI API access)
Entrez.email = "allaboutbioinformatics21@gmail.com"  # Replace with your email address

def fetch_gene_data(disease_name):
    """Fetch gene data from NCBI for a specific disease name."""
    try:
        print(f"Fetching data for disease: {disease_name}")
        # Search for gene information in NCBI
        handle = Entrez.esearch(db="gene", term=f"{disease_name}[All Fields]", retmax=375)  # Fetch up to 375 results
        record = Entrez.read(handle)
        handle.close()

        print("Record:", record)

        gene_data = []
        # Fetch data for each gene ID
        for gene_id in record["IdList"]:
            print(f"Found Gene ID: {gene_id}")
            handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(handle)
            handle.close()

            # Extracting details from the summary
            gene_name = summary[0].get("Name", "No name available")
            location_info = summary[0].get("GenomicInfo", [{}])[0]
            chromosome = location_info.get("Chromosome", "No chromosome available")
            location = f"{chromosome} {location_info.get('MapLocation', 'No location available')}"
            mim = summary[0].get("MIM", "No MIM available")
            ncbi_link = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"

            description = summary[0].get("Summary", "No summary available")
            gene_data.append((gene_id, gene_name, location, mim, ncbi_link, description))

        if gene_data:
            return gene_data
        else:
            print("No data found for the specified disease.")
            return [("No gene data found", "", "", "", "", "No data found for the specified disease.")]
    except Exception as e:
        print(f"Error fetching data: {e}")
        return [(None, None, None, None, None, None)]

def insert_gene_data():
    """Fetch gene data from NCBI and insert it into the SQLite database."""
    try:
        # Connect to the SQLite database
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        # Ensure the table exists with the additional columns
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
        print("Table creation ensured or already exists.")

        # Fetch and insert data for Parkinson's Disease
        disease = 'Parkinson\'s Disease'
        gene_data = fetch_gene_data(disease)
        if gene_data:
            for gene_id, gene_name, location, mim, ncbi_link, description in gene_data:
                if description and gene_id:
                    print(f"Inserting: {gene_id}, {gene_name}, {location}, {mim}, {ncbi_link}, {description}")  # Debug print
                    cursor.execute('''
                    INSERT INTO genomics_data (disease, gene_id, gene_name, location, mim, ncbi_link, description)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', (disease, gene_id, gene_name, location, mim, ncbi_link, description))
                else:
                    print(f"Skipping insertion for gene_id: {gene_id} due to missing description or ID")  # Debug print
        else:
            print("No gene data to insert.")  # Debug print

        # Commit changes and close the connection
        connection.commit()
        print("Gene data committed to the database.")
    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
    finally:
        if connection:
            connection.close()
            print("Database connection closed.")

if __name__ == "__main__":
    insert_gene_data()
