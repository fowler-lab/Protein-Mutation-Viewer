# Protein Mutation Viewer
A flask server to allow visualisation of protein mutations.
Uses a customised version of Mol* to provide a viewer.

## Installation
```
git clone --recurse-submodules git@github.com:JeremyWesthead/Protein-Mutation-Viewer.git
cd Protein-Mutation-Viewer

#Build the viewer
cd molstar-mutation
npm install
npm run build

#Install requirements for flask
cd ..
python -m virtualenv env
source env/bin/activate
pip install -r requirements.txt
```

## Run
This will start the webserver on port 4000
```
python app.py
```
## Documentation
### Python
Documentation for Python logic for website flow, as well as creating mutation PDB files can be found with `python -m pydoc -b app.py` and navigating to `Protein-Mutation-Viewer` module.
This skips a lot of the documentation of the flask endpoints, but this can be found in the associated docstrings to all of the endpoints.

### Molstar Mutation
#### Viewer
In order to use the viewer, a `div` must be created with a known ID, and styling to give the size of the div. In order to enable resizing of the viewer, extra logic has to be added. In the majority of use cases in this project, the logic to perform this is below:
```
<style>
    #app {
        position:relative;
        width: 60rem;
        height: 750px;
        float: right;
        resize: both;
        overflow: auto;
    }
</style>
<div id="app" onmouseup="resize_check();" onmouseover="being_resized=false;" onmousemove="being_resized=false;"></div>
<script type="text/javascript">
    var being_resized;
    var viewer;
    async function load_structure() {
        viewer.plugin.colour = window.viewer_colour;
        viewer.plugin.is_reference = false;
        await viewer.loadStructureFromUrl("/data/mutations/{{pdb}}.pdb", "pdb", false);
        viewer.plugin.is_reference = true;
        await viewer.loadStructureFromUrl("/data/reference-proteins/6vxx-blank.pdb", "pdb", false);
    }
    window.onload = async function(){
        viewer = new molstar.Viewer('app', {
            layoutIsExpanded: false,
            layoutShowRemoteState: false,
            layoutShowSequence: true,
            layoutShowLog: false,
            layoutShowLeftPanel: false,
            viewportShowExpand: false,
            viewportShowSelectionMode: false,
            viewportShowAnimation: false,
            layoutShowControls: true,
            viewportShowExpand: false,
            collapseLeftPanel: true,
            pdbProvider: 'pdbe',
            emdbProvider: 'pdbe',
        });
        await load_structure();
        const resize_ob = new ResizeObserver(function(entries) {
            // Set the bool to show the app has been resized
            being_resized = true;
        });
        resize_ob.observe(document.querySelector('#app'));
        being_resized = false;
    }
    async function resize_check(){
        if (being_resized === true){
            //The element has just been resized so redraw it
            viewer.handleResize();
            being_resized = false;
        }
    }
</script>
```
This creates a `div` with the id of `app`, which has event listeners for `mouseup`, `mouseover`, and `mousemove` to detect resizes (which is enabled for the `div` using the `resize: both;` line of css).

When the window loads, a new viewer is instanciated, passing the `div` id, as well as a variety of options.

`load_structure()` is an asynchronous function for the purpose of loading the structure. Initially, it sets the colour scheme by setting `viewer.plugin.colour`, and then denoting that the first structure is not a reference protein with `viewer.plugin.is_reference`. The viewer then loads the mutation structure (which only has the points of mutation). The viewer is then updated to say that the next structure is a reference with `viewer.plugin.is_reference = true;`, and then the blank structure is loaded. This is the reference `6VXX` PDB, but with the occupancy field set to 0 for all values. This results in the mutation points being loaded and displayed as spacefill with occupancy based colouring, and the reference structure being loaded as an all white protein in cartoon representation. 

#### Comparisons
Comparisons can be handled, using the cantor pairing function to create a single PDB file with 2 sets of mutations stored within it. To setup a viewer to show a comparison, some flags must be set:
```
viewer.plugin.is_comparison = true //Sets colouring to use inverse cantor pairing function
window.is_comparison = true //Updates onhover labels to detail both mutations
```