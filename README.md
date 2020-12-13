To install the required environment run:
> conda env create -f environment.yml
> conda activate pdga

Modify the run.py file and run it giving a seed as argument 
> python run.py 0

In this implementation you can use the MAP4 fingerprint or the RDKit AP fingerprint

Each run will create a folder (named according to the query name and the used seed) and will write out two files:
- one containing all sequences generated during each generation, their average distance from the query and the minimum distance from the query
- another containing all sequences below the defined threshold.

For more results run the algorithm with different random seeds.

To read the output you can use the jupyter notebook in this folder.

