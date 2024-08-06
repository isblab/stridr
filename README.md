[Add PubMed link]: [![PubMed](https://salilab.org/imp-systems/static/images/pubmed.png)](https://pubmed.ncbi.nlm.nih.gov/36040254/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11044598.svg)](https://zenodo.org/doi/10.5281/zenodo.11044598)


# StrIDR: Structures of Intrinsically Disordered Regions

StrIDR is a database of IDRs, confirmed via experimental or homology-based evidence, that are resolved in experimentally determined structures.  
Existing databases provide extensive information for IDRs at the sequence level. However, only a tiny fraction of these IDRs are associated with an experimentally determined protein structure. Moreover, the disordered region of interest could be unresolved even if a structure exists.  
StrIDR is expected to be useful for gaining insights into the dynamics, folding, and interactions of IDRs. It may provide structural data for computational methods to studying IDRs including molecular dynamics (MD) simulations and machine learning (ML).  

![Main_fig]()


## Directory structure
1. [Database](Database/) : contains the input files for database creation and the database.
2. [scripts](scripts/) : contains all the scripts used for creation of the database.


## Creating StrIDR
### Obtain data files from disorder databases and datasets
1. DIBS : https://dibs.enzim.ttk.mta.hu/downloads/
2. MFIB : https://mfib.pbrg.hu/downloads.php
3. FuzDB : https://fuzdb.org/api/entries?&sort_field=entry_id&sort_value=asc&format=tsv
4. PDBtot and PDBcdr : https://doi.org/10.1016/j.jmb.2020.02.017
5. DisProt : https://disprot.org/api/search?release=2023_12&show_ambiguous=true&show_obsolete=false&format=tsv&namespace=all&get_consensus=false
6. IDEAL : https://www.ideal-db.org/download/current/IDEAL.xml.gz  

Add the downloaded data files to the `Database/raw/` directory. The file names must match those specified in 1_disobind_database.py.


### Creating input files for StrIDR construction
```
python 1_download_databases.py -c CORES
```
CORES - number of cores to parallelize on.
This script generates the following files that are used by the downstream script:  
<ul>
	<li>Merged_Uniprot_IDs.txt : file containing UniProt accessions obtained from the aforementioned databases and datasets.</li>
	<li>Merged_PDB_IDs.txt : txt file containing PDB IDs for all UniProt accessions.</li>
</ul>


### Creating StrIDR database
```
python 2_create_database_dataset_files.py -c CORES
```


### Move files to the website directory
```
python move_to_web_dir.py
```


### Populate the SQL database and make migrations
```
cd ../website/dbsite/
python populate.py
python manage.py migrate
python manage.py makemigrations database
python manage.py migrate
```

