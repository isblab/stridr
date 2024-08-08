import subprocess
import time

tic = time.time()
# Move StrIDR JSON file to website/db_site/database/static/
print( "Moving StrIDR JSON file..." )
subprocess.call( 
	["cp", 
	"../Database/StrIDR_database/StrIDR_database.json",
	"../website/db_site/database/static/"]
 	)

print( "Moving Uniprot sequences JSON file..." )
# Move StrIDR UniProts JSON file to website/db_site/database/static/
subprocess.call( 
	["cp", 
	"../Database/StrIDR_database/Uniprot_seq.json",
	"../website/db_site/database/static/"]
 	)

# Move all CIF files to website/db_site/database/static/
print( "Moving CIF files to web dir..." )
subprocess.call( 
	["cp", "-r", 
	"../Database/Combined_PDBs/",
	"../website/db_site/database/static/Combined_PDBs/"]
 	)

subprocess.call( 
	["zip", 
	"../website/db_site/database/static/Combined_PDBs.zip",
	"../website/db_site/database/static/Combined_PDBs/"]
 	)

# Move all SIFTS mapping files to website/db_site/database/static/
print( "Moving SIFTS files to web dir..." )
subprocess.call( 
	["cp", "-r", 
	"../Database/Mapped_PDBs/",
	"../website/db_site/database/static/Mapped_PDBs/"]
 	)

subprocess.call( 
	["zip", 
	"../website/db_site/database/static/Mapped_PDBs.zip",
	"../website/db_site/database/static/Mapped_PDBs/"]
 	)


toc = time.time()
print( f"Total time taken = {(toc-tic)/60} minutes" )
print( "May the Force serve you well..." )

