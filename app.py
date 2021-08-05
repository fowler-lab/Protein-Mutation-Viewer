import os.path
from collections import Counter, defaultdict
from os import mkdir

import pandas as pd
from flask import Flask, redirect, render_template, url_for, request

REFERENCE_DIR = "data/reference-proteins/"
MUTATIONS_DIR = "data/mutations/"
DATASET = pd.read_pickle("data/MUTATIONS-Spike.pkl.gz")

def get_muatation_counts(mutations, threshold=5):
    '''Function to get the normalised number of mutations at each amino acid residue

    Args:
        mutations (list): List of mutations in GARC
        threshold (int, optional): Threshold to consider a mutation. Defaults to 5

    Returns:
        (dict, dict): Tuple of dictionaries mapping mutation_index->normalised_count and mutation_index->reference_amino_acid
    '''    
    assert len(mutations) > 0
    mutation_counts = defaultdict(int)
    references = {}
    mutants = defaultdict(list)
    #Checking for SNPs first
    for mutation in mutations:
        if "del" in mutation or "ins" in mutation or "indel" in mutation:
            mutation_counts[int(mutation.split("_")[0])] += 1
            mutants[int(mutation.split("_")[0])].append("_".join(mutation.split("_")[1:]))
        else:
            references[int(mutation[1:-1])] = str(mutation[0])
            mutation_counts[int(mutation[1:-1])] += 1
            mutants[int(mutation[1:-1])].append(mutation[-1])

    #Normalise the counts to fall between 0 and 1
    max_count = max(mutation_counts.values())

    def normalise(count, max_count):
        return round(count/max_count, 2)
        # return round(0.75* (count / max_count)**2 + 0.25, 2)

    mutation_counts = {index: normalise(mutation_counts[index], max_count) for index in mutation_counts.keys() if mutation_counts[index] >= threshold}
    #Ignore all values which have normalised counts of 0
    mutation_counts = {index : mutation_counts[index] for index in mutation_counts.keys() if mutation_counts[index] > 0}
    mutants = {index: Counter(mutants[index]) for index in mutants.keys() if index in mutation_counts.keys()}
    return mutation_counts, references, mutants

def get_counts_by_lineage(lineage="pango"):
    '''Get the number of mutations by lineage

    Args:
        lineage (str, optional): The name of the lineage type. One of ['pango', 'scorpio']
    Returns:
        dict: A counter dict to count the number of occurances of each lineage
    '''    
    if lineage == "scorpio":
        lineage = "scorpio_call"
    elif lineage == "pango":
        lineage = "pango_lineage"
    else:
        return {}
    return Counter([i for i in DATASET[lineage] if not pd.isnull(i)])

def write_new_pdb(mutation_counts, filename="out.pdb", reference=REFERENCE_DIR+"6vxx.pdb"):
    '''Writes a new pdb file with the occupancy field altered appropriately for the mutation_counts

    Args:
        mutation_counts (dict): Dictionary mapping mutation_index->normalised_count
        filename (str, optional): Output filename for the resultant pdb file. Defaults to "out.pdb".
        reference (str, optional): Filename of the reference pdb file. Defaults to REFERENCE_DIR+"6vxx.pdb".
    '''    
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
                    continue
                    line[56:60] = list("0.00")
                line = ''.join(line)
            f.write(line)


def write_reference_pdb(filename="out.pdb", reference=REFERENCE_DIR+"6vxx.pdb"):
    '''Writes a new pdb file with the occupancy field altered to be 0 consistently - 
    creating a blank cartoon representation

    Args:
        filename (str, optional): Output filename for the resultant pdb file. Defaults to "out.pdb".
        reference (str, optional): Filename of the reference pdb file. Defaults to REFERENCE_DIR+"6vxx.pdb".
    '''    
    filename = REFERENCE_DIR+filename
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
                line[56:60] = list("0.00")
                line = ''.join(line)
            f.write(line)

def run():
    '''Factory method for the Flask server
    '''    
    app = Flask(__name__, static_url_path='', static_folder='', template_folder='templates')

    @app.errorhandler(404)
    def page_not_found(e):
        return render_template('404.html')

    @app.route("/")
    def index():
        '''Render the index page
        '''  
        return render_template("index.html")
    
    @app.route("/viewer/covid/spike")
    def viewer_home():
        return render_template("viewer.html")
    
    @app.route("/viewer/covid/spike/<lin_type>")
    def viewer_lineage_home(lin_type):
        if lin_type in ["pango", "scorpio"]:
            return render_template('viewer.html', type=lin_type)
        else:
            return render_template('404.html'), 404
    
    @app.route("/viewer/covid/spike/<lin_type>/<lineage>")
    def view_lineage(lin_type, lineage):
        '''Page for viewing mutations in the covid spike protein

        Args:
            lineage (str): Name of the pango lineage
        '''        
        if lin_type == "pango":
            pass
            # lineage = lineage.upper()
        elif lin_type != "scorpio":
            #If the lineage type is not pango or scorpio, 404
            return render_template('404.html'), 404
        #Check if the pdb has already been made, and generate as required
        # if not os.path.isfile(MUTATIONS_DIR+lineage+".pdb"):
        try:
            if lin_type == "pango":
                mutation_counts, references, mutations = get_muatation_counts(DATASET[DATASET["pango_lineage"] == lineage]["mutation"])
            else:
                mutation_counts, references, mutations = get_muatation_counts(DATASET[DATASET["scorpio_call"] == lineage]["mutation"])
            write_new_pdb(mutation_counts, filename=lineage+".pdb")
        except AssertionError:
            #This mutation has no data
            return render_template('viewer.html', lineage=lineage)
        if not os.path.isfile(REFERENCE_DIR+"6vxx-blank.pdb"):
            write_reference_pdb(filename="6vxx-blank.pdb")
        return render_template("viewer.html", pdb=lineage, mutation_counts=mutation_counts, references=references, mutations=mutations, type=lin_type)
    
    @app.route("/list/<lin_type>")
    def list_lineages(lin_type):
        '''Render a page containing a list of lineages
        '''        
        if lin_type == "pango":
            lineages = sorted(list(set(DATASET["pango_lineage"])))
        elif lin_type == "scorpio":
            lineages = sorted(list(set([i for i in DATASET["scorpio_call"] if not pd.isnull(i)])))
        else:
            return render_template('404.html'), 404
        counts = get_counts_by_lineage(lin_type)
        return render_template("lineage_list.html", lineages=lineages, counts=counts, type=lin_type)
    
    @app.route("/search")
    def search():
        '''Render a page of results for a lineage search
        '''        
        query = request.args.get("query")
        type_restriction = request.args.get("type")
        if query:
            scorpio = list(set([i for i in DATASET["scorpio_call"] if not pd.isnull(i)]))
            pango = list(set(DATASET["pango_lineage"]))
            if type_restriction != "pango":
                scorpio_results = [i for i in scorpio if query.lower() in i.lower()]
            else:
                scorpio_results = []
            if type_restriction != "scorpio":
                pango_results = [i for i in pango if query.lower() in i.lower()]
            else:
                pango_results = []
            return render_template('search.html', query=query, scorpio_results=scorpio_results, pango_results=pango_results, type=type_restriction)
        else:
            return render_template('viewer.html', type="scorpio")




    app.run(debug=True, port=4000)

if __name__ == "__main__":
    #Ensure that the MUTATIONS_DIR and REFERENCE_DIR exist
    if not os.path.isdir(REFERENCE_DIR):
        mkdir(REFERENCE_DIR)
    if not os.path.isdir(MUTATIONS_DIR):
        mkdir(MUTATIONS_DIR)
    run()
