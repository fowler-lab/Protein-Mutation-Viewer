from flask import Flask, render_template
import os.path



import pandas as pd
from collections import defaultdict


REFERENCE_DIR = "../../data/reference-proteins/"
MUTATIONS_DIR = "../../data/mutations/"
DATASET = pd.read_pickle("../../data/MUTATIONS-Spike.pkl.gz")

def get_muatation_counts(mutations):
    assert len(mutations) > 0
    mutation_counts = defaultdict(int)
    #Checking for SNPs first
    for mutation in mutations:
        if "del" in mutation or "ins" in mutation or "indel" in mutation:
            mutation_counts[int(mutation.split("_")[0])] += 1
        else:
            mutation_counts[int(mutation[1:-1])] += 1

    #Normalise the counts to fall between 0 and 1
    max_count = max(mutation_counts.values())

    def normalise(count, max_count):
        return round(0.75* (count / max_count)**2 + 0.25, 2)

    mutation_counts = {index: normalise(mutation_counts[index], max_count) for index in mutation_counts.keys()}
    return mutation_counts

def write_new_pdb(mutation_counts, filename="out.pdb", reference=REFERENCE_DIR+"6vxx.pdb"):
    filename = MUTATIONS_DIR+filename
    #Read the reference pdb file
    with open(reference) as f:
        ref = [line for line in f]

    def trim(s):
        '''Takes in a string and returns it without trailing or leading whitespace

        Args:
            s (str): String
        '''
        if s[0] == " ":
            return trim(s[1:])
        if s[-1] == " ":
            return trim(s[:-1])
        return s    

    #Write new file
    with open(filename, "w") as f:
        for line in ref:
            if line[:5] == "ATOM ":
                line = list(line)
                #Match the ATOM to the amino acid number from the dataset
                aa_index = int(trim(''.join(line[22:26])))
                if aa_index in mutation_counts.keys():
                    #Make the occupancy the normalised count if mutated at this atom
                    line[56:60] = list(str(mutation_counts[int(trim(''.join(line[22:26])))]))
                else:
                    #Default to 0 if there were no mutations
                    line[56:60] = list("0.00")
                line = ''.join(line)
            f.write(line)

def run():
    app = Flask(__name__)

    @app.route("/")
    def index():
        return render_template("index.html")
    
    @app.route("/viewer/pango/<lineage>")
    def view_lineage(lineage):
        #Check if the pdb has already been made, and generate as required
        if not os.path.isfile(MUTATIONS_DIR+lineage+".pdb"):
            try:
                write_new_pdb(get_muatation_counts(DATASET[DATASET["pango_lineage"] == lineage]["mutation"]), filename=lineage+".pdb")
            except AssertionError:
                #This mutation has no data
                return "<h3>No data found for lineage "+lineage+"</h3>"
        return render_template("viewer.html", pdb=lineage)
    
    @app.route("/list/pango")
    def list_lineages():
        lineages = sorted(list(set(DATASET["pango_lineage"])))
        return render_template("pango_list.html", lineages=lineages)
    
    app.run(debug=True, port=4000)

if __name__ == "__main__":
    run()