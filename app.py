from flask import Flask, render_template
app = Flask(__name__)
import sqlite3

@app.route("/") 
#empty route for homepage
def home():
    return render_template('home.html')

@app.route("/genomics")
def genomics():
    return render_template('genomics.html')

@app.route("/Parkinsons_genomics")
def Parkinsons_genomics():
    return render_template('Parkinsons_genomics.html')

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
        connection = sqlite3.connect('database.db')
        cursor = connection.cursor()

        cursor.execute('''
        SELECT gene, description FROM genomics_data WHERE disease = 'Parkinson\'s Disease'
        ''')
        data = cursor.fetchall()
        connection.close()

        print("Data fetched successfully:", data)  # Debugging statement
        return render_template('parkinsons_genomics.html', data=data)
    except Exception as e:
        print("Error:", e)  # Print the error
        return "An error occurred while accessing the database."


print(__name__)
if __name__ == "__main__":
    app.run(host = '0.0.0.0', debug=True) 