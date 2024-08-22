import os
from omegaconf import OmegaConf
os.environ.setdefault( "DJANGO_SETTINGS_MODULE", "db_site.settings" )

import django
django.setup()

from database.models import PDB, Uniprot, UniprotPDB, ChainID, DisorderResidues, CrossRefs

def populate( db, N = 200 ):
	problematic = {}
	total = len( db.keys() )
	count = 0
	for idx, pdb_id in enumerate( db.keys() ):
		print( f"\nEntry {idx}/{total}: {pdb_id}" )
		pdb_obj = PDB( pdb_id = pdb_id )
		pdb_obj.save()

		# if pdb_id != "7qo6":
		# 	continue
		# else:
		# 	print( f"\nEntry {idx}/{total}: {pdb_id}" )

		for uni_id in db[pdb_id].keys():
			# print( db[pdb_id][uni_id]["Name"], "  ", uni_id )
			if db[pdb_id][uni_id]["Name"] == []:
				count += 1
				continue

			name = db[pdb_id][uni_id]["Name"][0]

			uniprot_obj, _ = Uniprot.objects.get_or_create( 
												uni_id = uni_id,
												name = name
												 )
			uniprot_obj.save()

			uni_pdb_obj, _ = UniprotPDB.objects.get_or_create( 
												pdb = pdb_obj,
												uniprot = uniprot_obj
												 )
			uni_pdb_obj.save()
			
			# Add PDB chains info.
			asym_id, auth_asym_id = [], []
			for ind in range( len( db[pdb_id][uni_id]["Asym_ID"] ) ):
				chain_id_obj, _ = ChainID.objects.get_or_create(  
												asym_id = db[pdb_id][uni_id]["Asym_ID"][ind],
												auth_asym_id = db[pdb_id][uni_id]["Auth_Asym_ID"][ind]
												)
				uni_pdb_obj.chain_ids.add( chain_id_obj )
				
				# Add disorder residue info.
				disorder_in_seq, disorder_in_struct = [], []
				for pos in db[pdb_id][uni_id]["Disorder_in_seq"][ind]:
					if pos[0] != pos[1]:
						disorder_in_seq.append( f"{pos[0]}-{pos[1]}" )
					else:
						disorder_in_seq.append( f"{pos[0]}" )
				disorder_in_seq = ", ".join( disorder_in_seq )

				for pos in db[pdb_id][uni_id]["Disorder_in_struct"][ind]:
					if pos[0] != pos[1]:
						disorder_in_struct.append( f"{pos[0]}-{pos[1]}" )
					else:
						disorder_in_struct.append( f"{pos[0]}" )
				disorder_in_struct = ", ".join( disorder_in_struct )

				disorder_res_obj, _ = DisorderResidues.objects.get_or_create( 
															disorder_in_seq = disorder_in_seq,
															disorder_in_struct = disorder_in_struct
															 )
				uni_pdb_obj.disorder_residues.add( disorder_res_obj )


			
			# Add cross references.

			disprot = list( set( db[pdb_id][uni_id]["cross_refs"]["DisProt"] ) )
			if "DP00303" in disprot:
				print( pdb_id )
				print( db[pdb_id] )
			ideal = list( set( db[pdb_id][uni_id]["cross_refs"]["IDEAL"] ) )
			mobidb = list( set( db[pdb_id][uni_id]["cross_refs"]["MobiDB"] ) )

			max_elements = max( [len( disprot ), len( ideal ), len( mobidb )] )
			len_disprot = len( disprot )
			len_ideal = len( ideal )
			len_mobidb = len( mobidb )

			# Making all 3 cross_refs of the same size.
			if len_disprot < max_elements:
				[disprot.append( "" ) for i in range( ( max_elements - len_disprot ) )]
			if len_ideal < max_elements:
				[ideal.append( "" ) for i in range( ( max_elements - len_ideal ) )]
			if len_mobidb < max_elements:
				[mobidb.append( "" ) for i in range( ( max_elements - len_mobidb ) )]

			for ind in range( max_elements ):
				cross_ref_obj, _ = CrossRefs.objects.get_or_create( 
															disprot = disprot[ind],
															ideal = ideal[ind],
															mobidb = mobidb[ind]
														 )
				uni_pdb_obj.cross_refs.add( cross_ref_obj )

			# exit()
			# Add Uniprot info to the PDB.
			pdb_obj.uniprot.add( uniprot_obj )
	print( count )


if __name__ == "__main__":
	import json
	db_path = "/data2/kartik/Disorder_Proteins/disobind/Database/StrIDR_database/"
	with open( f"{db_path}StrIDR_database.json" ) as f:
		db = json.load( f )
	
	populate( db )
	print( "May the Force be with you..." )
