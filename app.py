import os.path
from collections import Counter, defaultdict
from os import mkdir
import re

import pandas as pd
from flask import Flask, redirect, render_template, url_for, request, make_response

#Global definitions for important paths
REFERENCE_DIR = "data/reference-proteins/"
MUTATIONS_DIR = "data/mutations/"
DATASET = pd.read_pickle("data/MUTATIONS-Spike-Fixed.pkl.gz")

class MissingDataException(Exception):
    '''Exception class used to show when there was no data available for a given lineage
    '''    
    pass

def get_muatation_counts(mutations, threshold=5):
    '''Function to get the normalised number of mutations at each amino acid residue

    Args:
        mutations (list): List of mutations in GARC
        threshold (int, optional): Threshold to consider a mutation. Defaults to 5

    Returns:
        (dict, dict, dict): Tuple of dictionaries mapping mutation_index->normalised_count and mutation_index->reference_amino_acid, mutation_index->{mutation->frequency}
    '''    
    if len(mutations) == 0:
        raise MissingDataException("No data given")
    mutation_counts = defaultdict(int)
    references = {}
    mutants = defaultdict(list)
    #Checking for SNPs first
    for mutation in mutations:
        if "del" in mutation or "ins" in mutation or "indel" in mutation:
            mutation_counts[(int(mutation.split("_")[0])-1)//3 + 1] += 1
            mutants[(int(mutation.split("_")[0])-1)//3 + 1].append("_".join(mutation.split("_")[1:]))
        else:
            references[int(mutation[1:-1])] = str(mutation[0])
            mutation_counts[int(mutation[1:-1])] += 1
            mutants[int(mutation[1:-1])].append(mutation[-1])

    #Normalise the counts to fall between 0 and 1
    max_count = max(mutation_counts.values())

    def normalise(count, max_count):
        #Simple normalisation with rounding to 2 d.p
        return round(count/max_count, 2)

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

def trim(s):
    '''Takes in a string and returns it without trailing or leading whitespace

    Args:
        s (str): String
    '''
    if len(s) == 0:
        return s
    if s[0] == " ":
        return trim(s[1:])
    if s[-1] == " ":
        return trim(s[:-1])
    return s

def write_new_pdb(mutation_counts, filename="out.pdb", reference=REFERENCE_DIR+"6vxx.pdb", find_missing=False):
    '''Writes a new pdb file with the occupancy field altered appropriately for the mutation_counts

    Args:
        mutation_counts (dict): Dictionary mapping mutation_index->normalised_count
        filename (str, optional): Output filename for the resultant pdb file. Defaults to "out.pdb".
        reference (str, optional): Filename of the reference pdb file. Defaults to REFERENCE_DIR+"6vxx.pdb".
        find_missing (bool, optional): Boolean to determine whether this should return a list of 
                                        amino acid indices which were not added to the file (due to the reference 
                                        missing them)
    Returns:
        list: List of amino acid indices where the reference was not updated. Only returns if find_missing is set to True.
    '''    
    filename = MUTATIONS_DIR+filename
    #Read the reference pdb file
    with open(reference) as f:
        ref = [line for line in f]

    seen = set()
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
                    seen.add(int(trim(''.join(line[22:26]))))
                else:
                    #Default to 0 if there were no mutations
                    continue
                line = ''.join(line)
            f.write(line)
    if find_missing:
        return sorted(list(set(mutation_counts.keys()).difference(seen)))
    else:
        return None
        
def compare_mutations(mutations1, mutations2):
    '''Merge two dictionaries of mutations together using the cantor pairing function to give a single value.
    Cantor pairing function: (x, y)-> z -> (x, y) where x, y, z are all natural numbers.
        f(x, y) = 1/2 * (x + y) * (x + y + 1) + y = z
        f'(z) = (w - z + t, z - t) = (x, y) where
            w = int(((8 * z + 1)**0.5 - 1) / 2)
            t = w * (w + 1) / 2
    This is adapted here to convert the decimal values to natural numbers and vice versa for this purpose.
    This involves truncating the values to [0.0, 0.1, .. 0.9, 1.0] and then mutliplying by 10 to lie in [0, 1, ... 9, 10]
    This results in the cantor pairing function having a domain of [0, 220], 
        which can be divided by 100 to lie within [0.00, 9.99] to fit within the occupancy field -
        allowing reconstitution of the x and y values within the viewer's code

    Args:
        mutations1 (dict): Dictionary mapping index->mutation_count
        mutations2 (dict): Dictionary mapping index->mutation_count
    
    Returns:
        dict: Dictionary mapping index->value where value is the product of the modified cantor pairing function of the corresponding pair
    '''
    keys = set(mutations1.keys()).union(set(mutations2.keys()))
    merged = {}
    for key in keys:
        #Adjust the values to be natural numbers in range 0-10 (this looses some precision, but allows reduction of values)
        x = round(mutations1.get(key, 0), 1)*10
        y = round(mutations2.get(key, 0), 1)*10
        #The output of the cantor pairing function lies between 220-0, so divide by 100 to fit within 0.00-9.99
        val = round((0.5*(x+y)*(x+y+1)+y)/100, 2)
        merged[key] = val
    return merged

def merge_mutations(mutations1, mutations2):
    '''Merge together two amino acid mutations 

    Args:
        mutations1 (dict): Dictionary of index->{amino_acid->frequency}
        mutation2 (dict): Dictionary of index->{amino_acid->frequency}
    Returns:
        dict: Dictionary mapping index->{amino_acid->[freq1, freq2]}
    '''    
    keys = set(mutations1.keys()).union(mutations2.keys())
    merged = {}
    for key in keys:
        x = mutations1.get(key, dict())
        y = mutations2.get(key, dict())
        merged_aa = {}
        aa_keys = set(x.keys()).union(set(y.keys()))
        for aa_key in aa_keys:
            aa_x = x.get(aa_key, 0)
            aa_y = y.get(aa_key, 0)
            merged_aa[aa_key] = [aa_x, aa_y]
        merged[key] = merged_aa
    return merged

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
        '''Render the index page. Also used to deal with cookie logic.
            Cookie logic:
                There are 2 possible cookies.
                    1. `cookie_consent`: Cookie to store whether a user has consented to non-essential cookies.
                    2. `colour`: Cookie used to store the current colour scheme between viewers.
        '''  
        #Check for colour setting
        colour = request.args.get("colour")
        consent = request.args.get("consent")
        return_path = request.args.get("return_path")
        if colour or consent:
            resp = make_response(render_template("set_cookie.html", return_path=return_path))
            if colour:
                #Get the value of the colour
                if colour == "1":
                    value = "default"
                elif colour in ["default", "red", "green", "blue"]:
                    value = colour
                else:
                    #Give an errro if cookie is not in approved values
                    return render_template('400.html'), 400
                #Check for consent
                if request.cookies.get("cookie_consent") == "yes":
                    resp.set_cookie("colour", value=value, samesite="strict", path="/", max_age=60*60*24*365)
                else:
                    #Consent has not been given, so add a GET parameter for the return path showing colour
                    #This will not persist between pages, but will recolour the current viewer
                    if "?" in return_path:
                        #There are already GET parameters, so append
                        return_path = return_path + "&colour=" + value
                    else:
                        #No GET parameters, so add appropriately
                        return_path = return_path + "?colour=" + value
                    return render_template("set_cookie.html", return_path=return_path)
            if consent:
                if consent == "1":
                    value = "yes"
                elif consent == "0":
                    value = "no"
                else:
                    render_template('400.html'), 400
                resp.set_cookie("cookie_consent", value=value, samesite="strict", path="/", max_age=60*60*24*365)
            return resp
        else:
            return render_template("index.html")
    @app.route("/cookies")
    def cookies():
        '''Render an info page about the uses of cookies.
        Contains a table of values for each cookie, detailing if they are set and their value.
        '''        
        cookie_consent = request.cookies.get("cookie_consent")
        colour = request.cookies.get("colour")
        return render_template("cookie_policy.html", cookie_consent=cookie_consent, colour=colour)
    
    @app.route("/viewer/covid/spike")
    def viewer_home():
        '''Home page for the viewer. Shows a viewer of all mutations in the dataset, and a search bar for searching
            by lineage.
        '''
        lin_type1 = request.args.get("lin_type1")
        lineage1 = request.args.get("lineage1")
        if lin_type1 is not None and lineage1 is not None:
            return render_template("viewer.html",  lin_type1=lin_type1, lineage1=lineage1)
        else:
            if not os.path.exists("data/mutations/all.pdb"):
                mutation_counts, references, mutations = get_muatation_counts(DATASET["mutation"])
                write_new_pdb(mutation_counts, filename="all.pdb", find_missing=True)
            return render_template("viewer.html", pdb="all")
    
    @app.route("/viewer/covid/spike/<lin_type>/<lineage>")
    def view_lineage(lin_type, lineage):
        '''Page for viewing mutations in the covid spike protein. Displays a viewer, as well as a table of the 
            mutations and their possible values.

        Args:
            lin_type (str): Name of the type of lineage. Must be either `scorpio` or `pango`
            lineage (str): Name of the lineage
        '''        
        if lin_type not in  ["scorpio", "pango"]:
            #If the lineage type is not pango or scorpio, 404
            return render_template('404.html'), 404
        #Check if the pdb has already been made, and generate as required
        # if not os.path.isfile(MUTATIONS_DIR+lineage+".pdb"):
        try:
            if lin_type == "pango":
                mutation_counts, references, mutations = get_muatation_counts(DATASET[DATASET["pango_lineage"] == lineage]["mutation"])
            else:
                mutation_counts, references, mutations = get_muatation_counts(DATASET[DATASET["scorpio_call"] == lineage]["mutation"])
            missing = write_new_pdb(mutation_counts, filename=lineage+".pdb", find_missing=True)
        except MissingDataException:
            #This mutation has no data
            return render_template('viewer.html', lineage=lineage)
        if not os.path.isfile(REFERENCE_DIR+"6vxx-blank.pdb"):
            write_reference_pdb(filename="6vxx-blank.pdb")
        return render_template("viewer.html", pdb=lineage, mutation_counts=mutation_counts, references=references, 
                                mutations=mutations, type=lin_type, missing=missing)
    
    @app.route("/viewer/compare/covid/spike/")
    def comparison():
        '''Render a page for viewing a comparison of 2 lineages. Displays a viewer showing relative mutation for both
            lineages, as well as a table of relative mutations. Uses GET request arguments.

        Args:
            lin_type1 (str): Name of the first lineage type.
            lin_type2 (str): Name of the second lineage type.
            lineage1 (str): Name of the first lineage.
            lineage2 (str): Name of the second lineage. 
        '''        
        #Get the lineages to compare
        lin_type1 = request.args.get("lin_type1")
        lin_type2 = request.args.get("lin_type2")
        lineage1 = request.args.get("lineage1")
        lineage2 = request.args.get("lineage2")

        #Get the mutations for both
        if lin_type1 and lin_type2 and lineage1 and lineage2:
            unknown = []
            try:
                if lin_type1 == "pango":
                    mutation_counts1, references1, mutations1 = get_muatation_counts(DATASET[DATASET["pango_lineage"] == lineage1]["mutation"])
                else:
                    mutation_counts1, references1, mutations1 = get_muatation_counts(DATASET[DATASET["scorpio_call"] == lineage1]["mutation"])
            except MissingDataException:
                unknown.append((lin_type1, lineage1))
            try:
                if lin_type2 == "pango":
                    mutation_counts2, references2, mutations2 = get_muatation_counts(DATASET[DATASET["pango_lineage"] == lineage2]["mutation"])
                else:
                    mutation_counts2, references2, mutations2 = get_muatation_counts(DATASET[DATASET["scorpio_call"] == lineage2]["mutation"])
            except MissingDataException:
                unknown.append((lin_type2, lineage2))
            if len(unknown) > 0:
                return render_template('compare_viewer.html', unknown=unknown)
            #Merge them
            mutation_counts = compare_mutations(mutation_counts1, mutation_counts2)
            mutations = merge_mutations(mutations1, mutations2)
            missing = write_new_pdb(mutation_counts, filename=f"{lineage1}_{lineage2}.pdb", find_missing=True)
            return render_template('compare_viewer.html', mutation_counts=mutation_counts, lin_type1=lin_type1, lineage1=lineage1,
                                    lin_type2=lin_type2, lineage2=lineage2, references1=references1, references2=references2,
                                    mutation_counts1=mutation_counts1, mutation_counts2=mutation_counts2, mutations=mutations,
                                    missing=missing)
        else:
            return render_template('comparison.html')
        
    
    @app.route("/list/<lin_type>")
    def list_lineages(lin_type):
        '''Render a page containing a list of lineages filtered by lineage type.

        Args:
            lineage_type (str): Name of the lineage type. Must be either `scorpio` or `pango`
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
        '''Render a page for search results. Has several uses:
            1. A query is given -> render results page
            2. A query is given with a lineage type restriction -> render results page with only that lineage type
            3. A comparison lineage is given -> render a search page to find a lineage to compare against
            4. A comparison lineage is given and a query is given -> render a results page to find a lineage to compare against
            Args are given through GET request arguments.
        
        Args:
            query (str): A lineage query
            type (str): Name of the lineage type to restrict the search to
            lin_type1 (str): Name of the lineage type to compare against
            lineage1 (str): Name of the lineage to compare against
            compare (int): Flag to show whether a comparison is occuring. 
        '''        
        query = request.args.get("query")
        type_restriction = request.args.get("type")
        lin_type1 = request.args.get("lin_type1")
        lineage1 = request.args.get("lineage1")
        compare = bool(request.args.get("compare"))
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
            if lin_type1 and lineage1:
                return render_template('search.html', query=query, scorpio_results=scorpio_results, pango_results=pango_results, 
                                        type=type_restriction, lin_type1=lin_type1, lineage1=lineage1, compare=compare)
            else:
                return render_template('search.html', query=query, scorpio_results=scorpio_results, pango_results=pango_results, 
                                        type=type_restriction, compare=compare)

        else:
            if lin_type1 and lineage1:
                return render_template('viewer.html', lin_type1=lin_type1, lineage1=lineage1, compare=compare)
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
