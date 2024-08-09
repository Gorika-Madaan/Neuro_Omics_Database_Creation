from flask import Flask, render_template
import sqlite3

app = Flask(__name__)


@app.route("/") 
#empty route for homepage
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

@app.route("/parkinsons_genomics")
def parkinsons_genomics():
    try:
        # Connect to the SQLite database
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        # Use parameterized query to avoid syntax errors
        query = '''
        SELECT gene, description FROM genomics_data WHERE disease = ?
        '''
        cursor.execute(query, ('Parkinson\'s Disease',))
        data = cursor.fetchall()  # Fetch all rows

        # Close the database connection
        connection.close()

        print("Data fetched successfully:", data)  # Debugging statement

        # Render the HTML template with the data
        return render_template('parkinsons_genomics.html', data=data)
    except Exception as e:
        # Print and return the error if any
        print("Error:", e)
        return f"An error occurred while accessing the database: {e}"


print(__name__)
if __name__ == "__main__":
    app.run(host = '0.0.0.0', debug=True) 