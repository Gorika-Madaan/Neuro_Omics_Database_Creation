from flask import Flask, render_template
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

print(__name__)
if __name__ == "__main__":
    app.run(host = '0.0.0.0', debug=True) 