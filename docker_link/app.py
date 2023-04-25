from flask import Flask
import presence_absence
app = Flask(__name__)

@app.route('/')
def initialise():
    title = "Mutation-Gene KO combinations with significant p values <br />"
    output = presence_absence.run_example().to_html()
    return title + output
