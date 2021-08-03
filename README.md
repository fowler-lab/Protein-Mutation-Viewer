# Protein-Mutation-Viewer
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