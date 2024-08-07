import json
import re
import time
import os
import ast
import logging
from flask import (
    Flask,
    render_template,
    request,
    jsonify,
    Response,
    stream_with_context,
)
import requests
from requests.adapters import HTTPAdapter, Retry

# Import process from process.py
from process import process

app = Flask(__name__)


def stream_data(data):
    return f"data: {json.dumps(data)}\n\n"


@app.route("/")
def index():
    return render_template("index2.html")


@app.route("/process", methods=["GET", "POST"])
def process_route():
    gene_name = request.args.get("gene_name")
    residue1 = request.args.get("residue1")
    position = int(request.args.get("position"))
    residue2 = request.args.get("residue2")

    # Validate input
    if not gene_name:
        return "Gene name is required", 400
    if not residue1:
        return "Residue 1 is required", 400
    if not position:
        return "Position is required", 400
    if not residue2:
        return "Residue 2 is required", 400

    def generate():
        yield stream_data({"message": "Processing request."})
        for data in process(gene_name, residue1, position, residue2):
            print(data)
            yield stream_data(data)
        yield stream_data({"type": "done", "message": "Done processing."})

    return Response(stream_with_context(generate()), content_type="text/event-stream")
