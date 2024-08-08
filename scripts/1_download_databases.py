######### Parse sequence and structure databases to get  #########
######### lists of PDB IDs containing IDRs in complex      ##########
######### ------>"May the Force serve u well..." <------##########
##################################################################

############# One above all #############
##-------------------------------------##

import numpy as np
import pandas as pd
import random
from datetime import date
import argparse
import tqdm

import json
from functools import partial
from multiprocessing import Pool, get_context
import tqdm

from from_APIs_with_love import *
import xml.etree.ElementTree as ET

import warnings
warnings.filterwarnings( "ignore" )
random.seed( 11 )


##########################################################################
#------------------------------------------------------------------------#
class TheIlluminati():
	"""
	Get all UniProt IDs containing IDRs from the following sources:
		DIBS, MFIB, FuzDB, PDBtot, PDBcdr,
		DisProt, IDEAL, MobiDB
	"""

	def __init__( self, dir_name, cores ):
		"""
		Constructor. 
		The raw files (self.dibs, self.mfib, self.fuzzdb) etc 
		are downloaded as is from the web servers or papers of the databases. 
		The output_file containts a list of PDBs and associated information parsed from the raw files.
		"""

		self.base_path = dir_name # "./input_files/"
		self.dibs_path    = "./raw/DIBS_complete_19Jul24.xml"                 # DIBS XML file  path.
		self.mfib_path    = "./raw/MFIB_complete_xml_19Jul24/xml/"            # MFIB XML files dir path.
		self.fuzdb_path  = "./raw/browse_fuzdb_19Jul24.tsv"           # Fuzdb tsv file path.
		self.pdb_tot_path = "./raw/FP_pdbtot_modified.xlsx"           # PDB TOT xlsx file path.
		self.pdb_cdr_path = "./raw/FP_pdbcdr_modified.xlsx"           # PDB CDR xlsx file path.
		self.disprot_path = "./raw/DisProt release_2024_06 with_ambiguous_evidences.tsv"      # DisProt tsv file  path.
		self.ideal_path = "./raw/IDEAL_19Jul24.xml"                   # IDEAL XML file  path.


		# Output files contaiing relevant info. from each database.
		self.dibs_output_file = f"{self.base_path}DIBS.txt"
		self.mfib_output_file = f"{self.base_path}MFIB.txt"
		self.fuzdb_output_file = f"{self.base_path}FuzDB.txt"
		self.pdb_tot_cdr_output_file = f"{self.base_path}PDB_tot_cdr.txt"
		self.disprot_output_file = f"{self.base_path}DisProt.csv"
		self.ideal_output_file = f"{self.base_path}IDEAL.csv"
		self.mobidb_output_file = f"{self.base_path}MobiDB.csv"
		self.merged_uniprot_ids = f"{self.base_path}Merged_Uniprot_IDs.txt"
		
		self.mobidb_batch_size = 1000      # No. of entries per batch to be downloaded.
		self.max_cores = cores
		self.max_trials = 25               # max no. of attempts to download.
		self.wait_time = 15
		

		self.pdb_tot_cdr_csv = f"{self.base_path}PDB_tot_cdr.csv"
		self.output_file = f"{self.base_path}Merged_PDB_IDs.txt"

		self.logger_file = f"{self.base_path}Logs_database"
		self.logger = {key:{} for key in ["counts", "time_taken"]}


	def forward( self ):
		self.tim = time.time()
		"""
		Extract PDB IDs and associated info from the five datasets 
		on structures of IDR complexes. Combine them into one dataframe/CSV file. 
		""" 
		self.run_checks()
		if os.path.exists( f"{self.logger_file}.json" ):
			with open( f"{self.logger_file}.json", "r" ) as f:
				self.logger = json.load( f )
		
		if os.path.exists( f"{self.output_file}.json" ):
			with open( f"{self.output_file}.json", "r" ) as f:
				self.logger = json.load( f )

		unique_uni_ids = self.process_databases_datasets()
		self.get_all_unique_pdbs( unique_uni_ids )

		self.save_logs()



	def run_checks( self ):
		"""
		Check if the raw files for all databases and datasets are 
			present in the required directory.

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		None
		"""
		if not os.path.exists( self.dibs_path ):
			raise Exception( "DIBS database file not found in Database/raw/..." )
		
		if not os.path.exists( self.mfib_path ):
			raise Exception( "MFIB database file not found in Database/raw/..." )
		
		if not os.path.exists( self.fuzdb_path ):
			raise Exception( "FuzzDB database file not found in Database/raw/..." )
		
		if not os.path.exists( self.pdb_tot_path ):
			raise Exception( "PDBtot dataset file not found in Database/raw/..." )
		
		if not os.path.exists( self.pdb_cdr_path ):
			raise Exception( "PDBcdr dataset file not found in Database/raw/..." )

		if not os.path.exists( self.disprot_path ):
			raise Exception( "DisProt database file not found in Database/raw/..." )
		
		if not os.path.exists( self.ideal_path ):
			raise Exception( "IDEAL dataset file not found in Database/raw/..." )


	def process_databases_datasets( self ):
		"""
		Obtain IDR UniProt IDs from the following databases and datasets:
			DIBS
			MFIB
			FuzDB
			DisProt
			IDEAL
			MobiDB
			PDBtot and PDBcdr

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		unique_uni_ids -- pd.Series of all unique UniProt IDs obtained.
		"""
		t_start = time.time()
		##------------------------------------
		print( "\nDIBS --------------" )
		if os.path.exists( self.dibs_output_file ):
			with open( self.dibs_output_file, "r" ) as f:
				dibs_uni_ids = f.readlines()[0].split( "," )
		else:
			dibs_uni_ids = self.parse_dibs()
		print( f"UniProt IDs from DIBS = {self.logger['counts']['dibs']}" )
		print( f"Time taken = {self.logger['time_taken']['dibs']} seconds" )

		##------------------------------------
		print( "\nMFIB --------------" )
		if os.path.exists( self.mfib_output_file ):
			with open( self.mfib_output_file, "r" ) as f:
				mfib_uni_ids = f.readlines()[0].split( "," )
		else:
			mfib_uni_ids = self.parse_mfib()
		print( f"UniProt IDs from MFIB = {self.logger['counts']['mfib']}" )
		print( f"Time taken = {self.logger['time_taken']['mfib']} seconds" )
		
		##------------------------------------
		print( "\nFuzDB --------------" )
		if os.path.exists( self.fuzdb_output_file ):
			with open( self.fuzdb_output_file, "r" ) as f:
				fuzdb_uni_ids = f.readlines()[0].split( "," )
		else:
			fuzdb_uni_ids = self.parse_fuzdb()
		print( f"UniProt IDs from FuzDB = {self.logger['counts']['fuzdb']}" )
		print( f"Time taken = {self.logger['time_taken']['fuzdb']} seconds" )

		##------------------------------------
		print( "\nPDBtot-PDBcdr --------------" )
		if os.path.exists( self.pdb_tot_cdr_output_file ):
			with open( self.pdb_tot_cdr_output_file, "r" ) as f:
				# print( f.readlines() )

				tot_uni_ids, cdr_uni_ids = f.readlines()
				tot_uni_ids = tot_uni_ids.strip().split( "," )
				cdr_uni_ids = cdr_uni_ids.strip().split( "," )

		else:
			tot_uni_ids, cdr_uni_ids = self.parse_pdb_tot_cdr()
		print( f"UniProt IDs from: PDBtot = {self.logger['counts']['pdb_tot']} \t PDBcdr = {self.logger['counts']['pdb_cdr']}" )
		print( f"Time taken = {self.logger['time_taken']['tot_cdr']} seconds" )

		##------------------------------------
		print( "\nDisProt --------------" )
		if os.path.exists( self.disprot_output_file ):
			disprot = pd.read_csv( self.disprot_output_file )
		else:
			disprot = self.parse_disprot()
		print( f"Total entries obtained from from DisProt = {self.logger['counts']['disprot']}" )
		print( f"Unique entries obtained from from DisProt = {self.logger['counts']['disprot_unique']}" )
		print( f"Total Uniprot IDs obtained from from DisProt = {self.logger['counts']['disprot_uniIDs']}" )
		print( f"Unique Uniprot IDs obtained from from DisProt = {self.logger['counts']['disprot_unique_uniIDs']}" )
		print( f"Time taken = {self.logger['time_taken']['disprot']} seconds" )

		##------------------------------------
		print( "\nIDEAL --------------" )
		if os.path.exists( self.ideal_output_file ):
			ideal = pd.read_csv( self.ideal_output_file )
		else:
			ideal = self.parse_ideal()

			with open( f"{self.logger_file}.json", "w" ) as w:
				json.dump( self.logger, w )

		print( f"Total entries obtained from from IDEAL = {self.logger['counts']['ideal']}" )
		print( f"Total Uniprot IDs obtained from from IDEAL = {self.logger['counts']['ideal_uniIDs']}" )
		print( f"Unique Uniprot IDs obtained from from IDEAL = {self.logger['counts']['ideal_unique_uniIDs']}" )
		print( f"Time taken = {self.logger['time_taken']['ideal']} seconds" )

		##------------------------------------
		print( "\nMobiDB --------------" )
		if os.path.exists( self.mobidb_output_file ):
			mobidb = pd.read_csv( self.mobidb_output_file )
		else:
			mobidb = self.parse_mobidb()

			with open( f"{self.logger_file}.json", "w" ) as w:
				json.dump( self.logger, w )

		print( f"Total Uniprot IDs obtained from MobiDB = {self.logger['counts']['mobidb']}" )
		print( f"Unique Uniprot IDs obtained from MobiDB = {self.logger['counts']['mobidb_unique']}" )
		print( f"Time taken = {self.logger['time_taken']['mobidb']} hours\n" )


		##------------------------------------##
		all_uni_ids = dibs_uni_ids + mfib_uni_ids + fuzdb_uni_ids + tot_uni_ids + cdr_uni_ids
		all_uni_ids += list( disprot["Uniprot ID"] )
		all_uni_ids += list( ideal["Uniprot ID"] )
		all_uni_ids += list( mobidb["Uniprot ID"] )
		# Just to separate some entries with comma-separated Uni IDs.
		all_uni_ids = ",".join( all_uni_ids ).split( "," )

		unique_uni_ids = pd.unique( all_uni_ids )

		t_end = time.time()
		if not os.path.exists( self.merged_uniprot_ids ):
			self.logger["counts"]["total_uni_ids"] = len( all_uni_ids )
			self.logger["counts"]["unique_uni_ids"] = len( unique_uni_ids )
			self.logger["time_taken"]["unique_uni_ids"] = ( t_end - t_start )

			print( f"Total UniProt IDs obtained: {len( all_uni_ids )}" )
			print( f"Unique UniProt IDs obtained: {len( unique_uni_ids )}" )
			print( f"Total time taken to get unique Uniprot IDs = {( t_end - t_start )/3600}" )

			with open( self.merged_uniprot_ids, "w" ) as w:
				w.writelines( ",".join( unique_uni_ids ) )

			with open( f"{self.logger_file}.json", "w" ) as w:
				json.dump( self.logger, w )
		
		else:
			print( "Merged Uniprto IDs already downloaded..." )
			with open( self.merged_uniprot_ids, "r" ) as f:
				unique_uni_ids = f.readlines()[0].split( "," )

			print( f"Total UniProt IDs obtained: {self.logger['counts']['total_uni_ids']}" )
			print( f"Unique UniProt IDs obtained: {self.logger['counts']['unique_uni_ids']}" )
			print( f"Total time taken to get unique Uniprot IDs = {self.logger['time_taken']['unique_uni_ids']/3600}" )


		return unique_uni_ids



	def parse_dibs( self ):
		"""
		Gather relevant info. from MFIB XML files.
		Store all PDB IDs in DIBS on disk.

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		dibs_pdbs --> list of all UniProt IDs in DIBS.
		"""
		tic = time.time()
		dibs_uni_ids = []
		root = ET.parse( self.dibs_path ).getroot()

		for i, parent in enumerate( root ):
			if parent.tag == "entry":
				for child in parent:
					if child.tag == "macromolecules":
						for subchild in child:
							if subchild.tag == "chain":
								if subchild.find( "type" ).text == "Disordered":
									dibs_uni_ids.append( subchild.find( "uniprot/id" ).text )

		with open( self.dibs_output_file, "w" ) as w:
			w.writelines( ",".join( dibs_uni_ids ) )

		toc = time.time()

		self.logger["counts"]["dibs"] = len( dibs_uni_ids )
		self.logger["time_taken"]["dibs"] = ( toc - tic )

		return dibs_uni_ids


	def parse_mfib( self ):
		"""
		Gather relevant info. from MFIB XML files.
		Store all PDB IDs in MFIB on disk.

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		mfib_pdbs --> list of all PDB IDs in MFIB.
		""" 
		tic = time.time()
		mfib_uni_ids = []
		for file in os.listdir( self.mfib_path ):
			root = ET.parse( f"{self.mfib_path}{file}" ).getroot()

			for i, parent in enumerate( root ):
				if parent.tag == "macromolecules":
					for child in parent:
						if child.tag == "chain":
							for subchild in child:
								if subchild.tag == "uniprot":
									mfib_uni_ids.append( subchild.find( "id" ).text )

		with open( self.mfib_output_file, "w" ) as w:
			w.writelines( ",".join( mfib_uni_ids ) )

		toc = time.time()

		self.logger["counts"]["mfib"] = len( mfib_uni_ids )
		self.logger["time_taken"]["mfib"] = ( toc - tic )

		return mfib_uni_ids



	def parse_fuzdb( self ):
		"""
		Gather relevant info. from MFIB XML files.
		Store all PDB IDs in DIBS on disk.

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		fuzdb_uni_ids --> list of all UniProt IDs in FuzDB.
		""" 
		tic = time.time()
		fuzdb_uni_ids = []
		
		data = pd.read_csv( self.fuzdb_path, sep="\t" )
		print( "Total entries in Fuzzdb = ", len( data ) )

		# Loop over all the entries in the file
		for i in range( len( data ) ):
			fuzdb_uni_ids.append( data.iloc[i,1] )

		with open( self.fuzdb_output_file, "w" ) as w:
			w.writelines( ",".join( fuzdb_uni_ids ) )

		toc = time.time()

		self.logger["counts"]["fuzdb"] = len( fuzdb_uni_ids )
		self.logger["time_taken"]["fuzdb"] = ( toc - tic )

		return fuzdb_uni_ids


	def parse_pdb_tot_cdr( self ):
		"""
		Gather relevant info. from PDBtot and PDBcdr.

		Input:
		----------
		Does not take any argument.

		Returns:
		----------
		tot_cdr_uni_ids --> list of all UniProt IDs in PDBtot and PDBcdr.
		""" 
		tic = time.time()
		tot_uni_ids, cdr_uni_ids = [], []

		pdbtot = pd.read_excel( self.pdb_tot_path )

		for i in range( len( pdbtot ) ):
			tot_uni_ids.append( pdbtot["Uniprot ID"][i] )

		pdbcdr = pd.read_excel( self.pdb_cdr_path )
		for i in range( len( pdbcdr ) ):
			cdr_uni_ids.append( pdbcdr["Uniprot ID"][i] )

		with open( self.pdb_tot_cdr_output_file, "w" ) as w:
			w.writelines( ",".join( tot_uni_ids ) + "\n" )
			w.writelines( ",".join( cdr_uni_ids ) )

		toc = time.time()

		self.logger["counts"]["pdb_tot"] = len( tot_uni_ids )
		self.logger["counts"]["pdb_cdr"] = len( cdr_uni_ids )
		self.logger["time_taken"]["tot_cdr"] = ( toc - tic )

		return tot_uni_ids, cdr_uni_ids



	def parse_disprot( self ):
		""" 
		DisProt Parser
		Parse the DisProt database tsv file to extract relevant info.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		df_new --> pd.DataFrame containing DisProt ID, UniProt ID, 
			and Disorder regions for all entries.
		"""
		 
		df = pd.read_csv( self.disprot_path, sep = "\t" )

		tic = time.time()
		disprot_holocron = { "Disprot ID":[], "Uniprot ID":[], "Disorder regions":[] }
		for i, id_ in enumerate( df["disprot_id"] ):
			disprot_holocron["Disprot ID"].append( id_ )
			disprot_holocron["Uniprot ID"].append( df["acc"][i] )

			start = df["start"][i]
			end = df["end"][i]	
			disprot_holocron["Disorder regions"].append( f"{start}-{end}" )

		df_new = pd.DataFrame()
		for key in disprot_holocron.keys():
			df_new[key] = disprot_holocron[key]
		df_new.to_csv( self.disprot_output_file, index = False )

		toc = time.time()

		self.logger["counts"]["disprot"] = len( df_new["Disprot ID"] )
		self.logger["counts"]["disprot_unique"] = len( set( df_new["Disprot ID"] ) )
		self.logger["counts"]["disprot_uniIDs"] = len( ",".join( df_new["Uniprot ID"] ).split( "," ) )
		self.logger["counts"]["disprot_unique_uniIDs"] = len( set( ",".join( df_new["Uniprot ID"] ).split( "," ) ) )
		self.logger["time_taken"]["disprot"] = ( toc - tic )

		return df_new



	def get_obsolete_uni_id( self, uni_ids ):
		"""  
		Check if the uniprot ID is obsolete or not.

		Input:
		----------
		uni_ids --> comma separated string of uniprot IDs.

		Returns:
		----------
		2 tuples, one of active IDs and another of obsolete IDs.  
		"""    
	 
		uni_ids = uni_ids.split( "," )
		obsolete, active = [], []
		num_uni_ids = len( uni_ids )

		# cores = max_cores if num_uni_ids > self.max_cores else num_uni_ids
		with Pool( self.max_cores ) as p:
			for result in p.imap_unordered( 
								partial( get_uniprot_seq, 
									max_trials = self.max_trials, 
									wait_time = self.wait_time, 
									return_id = True ), 
								uni_ids ):
				id_, seq = result
				if seq == []:
					obsolete.append( id_ )
				else:
					active.append( id_ )

		return ",".join( active ), ",".join( obsolete )


	def parse_ideal( self ):
		""" 
		IDEAL Parser
		Extract Uniprot ID and other relevant info from IDEAL xml file.
		Considering all disordered regions for verified ProS.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		df_new --> pd.DataFrame containing IDEAL ID, UniProt ID, 
			and Disorder regions for all selected entries.
		"""

		root = ET.parse( self.ideal_path ).getroot()
		
		tic = time.time()
		ideal_holocron = {i:[] for i in ["IDP ID", "Uniprot ID", "Disorder regions", "Uniprot ID (obsolete)"]}
		
		for i, tag in enumerate( root.findall( f"IDEAL_entry" ) ):
			verified = False
			idp_id = tag.find(f"idp_id").text

			# if idp_id == "IID00165":
			# Get all uniprot IDs for the entry.
			uni_ids = ",".join( [uni_id.text for uni_id in tag.findall(f"General/uniprot") ] )

			# Get all verified ProS disorder regions.
			pros_type = [t_.text for t_ in tag.findall( "Function/pros_type" )]
			start = [ds.text for ds in tag.findall( "Function/disorder_location/disorder_region_start" )]
			end = [de.text for de in tag.findall( "Function/disorder_location/disorder_region_end" )] 

			tmp_dr = []

			for idx, type_ in enumerate( pros_type ):
				if type_ == "verified" and start[idx] != None and end[idx] != None:
					# print( f"IDP ID {idp_id} contains a verified ProS..." )
					verified = True
					tmp_dr.append( f"{start[idx]}-{end[idx]}" )

			# Selecting only disorder regions for verified ProS.
			if verified:
				ideal_holocron["IDP ID"].append( idp_id )
				
				# Segregate the obsolete and active uniprot IDs.
				active, obsolete = self.get_obsolete_uni_id( uni_ids )

				# active, obsolete = uni_ids, ""
				ideal_holocron["Uniprot ID"].append( active )
				ideal_holocron["Uniprot ID (obsolete)"].append( obsolete )		
				
				ideal_holocron["Disorder regions"].append( ",".join( tmp_dr ) )
			print( f"{i} --> Done for IDP ID: {idp_id} --------------------------" )

		df = pd.DataFrame()
		for key in ideal_holocron.keys():
			df[key] = ideal_holocron[key]
		df.to_csv( self.ideal_output_file, index = False )

		toc = time.time()
		self.logger["counts"]["ideal"] = len( pd.unique( df["IDP ID"] ) )
		self.logger["counts"]["ideal_uniIDs"] = len( ",".join( df["Uniprot ID"] ).split( "," ) )
		self.logger["counts"]["ideal_unique_uniIDs"] = len( set( ",".join( df["Uniprot ID"] ).split( "," ) ) )
		self.logger["time_taken"]["ideal"] = ( toc - tic )

		return df


	def mobidb_parallel( self, partition, category ):
		"""
		Download a subset (partition) of MobiDB entries.
		Input:
		----------
		partition --> contains the batch to be downloaded; contains the limit and skip values.
		category --> specifies the category of data to be downloaded from MobiDB.

		Returns:
		----------
		dict containing Uniprot ID, disordered regions, and their annotation. 
		"""

		skip, limit = partition
		current = f"{skip}-{skip+limit}"
		# print( f"\nPartition {current}..." )

		# Randomly delay downloading different partitions by rn seconds; rn ~ U(a,b).
		rn = random.uniform( 1, 2 )
		time.sleep( rn )
		
		# Considering only disorder regions with a priority consensus.
		fields_included = [
							"curated-disorder-priority", "curated-lip-priority", 
							"homology-disorder-priority", "homology-lip-priority",
							"derived-lip-priority", "derived-mobile-th_90",
							"derived-mobile_context_dependent-th_90",
							"derived-binding_mode_context_dependent-priority",
							"curated-binding_mode_disorder_to_disorder-priority",
							"derived-binding_mode_disorder_to_disorder-priority",
							"homology-binding_mode_disorder_to_disorder-priority",
							"derived-binding_mode_disorder_to_order-priority",
							]
		
		download_url = f"https://mobidb.org/api/download?format=json&limit={limit}&skip={skip}&{category}=exists"
		
		"""
		Download the batch and check if the first entry is loadable or not (low pass test).
		Loop over all entries and load the JSON dict.
			Obtain the Uniprot ID and Uniref ID. 
			Also obtain the priority disorder, lip regions from curated and homology evidence.
		If it fails at any point, redownload the batch.
		For each redownload, start from the last unsuccessfully loaded entry instead of completly restarting again.
		"""
		last_outpost = 0
		data_dict = {}
		for trial in range( self.max_trials ):
			response = send_request( download_url, _format = None, max_trials = self.max_trials )

			try:
				data = json.loads( next( response.iter_lines() ) )
				
				all_entries = [entry for entry in response.iter_lines()]
				jump_start = last_outpost
				jump_end = len( all_entries )

				for i,entry in enumerate( all_entries[jump_start:jump_end] ):

					data = json.loads( entry )

					i += jump_start  
					data_dict[f"{skip+i}"] = {key:[] for key in ["Uniprot ID", "Disorder regions", "Annotation"]}

					data_dict[f"{skip+i}"]["Uniprot ID"].append( data["acc"] )

					for key in fields_included:
						if key in data.keys():
							data_dict[f"{skip+i}"]["Disorder regions"].extend( [f"{x[0]}-{x[1]}" for x in data[key]["regions"]] )
							data_dict[f"{skip+i}"]["Annotation"].append( key )
					last_outpost += 1
				break

			except Exception as e:
				if trial != self.max_trials-1:
					# print( f"Handling the error in \t Error in partition {current}...\n" )
					time.sleep( trial*1 )
					continue
				else:
					print( "Nooooooooooooooooooo....\n" )
					exit()

		return data_dict	


	def parse_mobidb( self ):
		""" 
		MobiDB Parser
		Obtain Uniprot IDs, disorder region and annotation from MobiDB
		This calls mobidb_parallel to download MobiDB entries in batches in parallel. 

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		df_new --> pd.DataFrame containing UniProt ID, Disorder regions, 
			and Annotation osurce for all entries.
		"""

		tim = time.time()
		mobidb = {i:[] for i in ["Uniprot ID", "Disorder regions", "Annotation"]}
		categories_included = [ "derived-lip-merge", "derived-mobile-th_90", 
								"derived-binding_mode_disorder_to_disorder-mobi", "derived-binding_mode_disorder_to_order-mobi",
								"curated-disorder-merge", "curated-lip-merge", "homology-disorder-merge", "homology-lip-merge"]
		
		for category in categories_included:
			tic = time.time()
			print( f"Downloading data for category: {category}" )
			count_url = f"https://mobidb.org/api/count?format=json&{category}=exists"

			response = send_request( count_url, _format = "json" )
			total = response["n"]

			print( f"Total entries in Category - {category} = {total}" )

			if total < self.mobidb_batch_size:
				skip, limit = [0], [total]
			else:
				skip = np.arange( 0, total, self.mobidb_batch_size )
				limit = [self.mobidb_batch_size for i in range( len( skip )-1 )]
				last = ( total - skip[-1] )
				limit = np.append( limit, last )

			partitions = list( zip( skip, limit ) )


			with Pool( 10 ) as p:
				results = tqdm.tqdm( p.imap_unordered( 
											partial( self.mobidb_parallel, category = category ), partitions ), 
										total = len( partitions ) )

				for result in results:
					for key in result.keys():
						if result[key]["Disorder regions"] != []:
							mobidb["Uniprot ID"].append( ",".join( result[key]["Uniprot ID"] ) )
							mobidb["Disorder regions"].append( ",".join( result[key]["Disorder regions"] ) )
							mobidb["Annotation"].append( ",".join( result[key]["Annotation"] ) )

			toc = time.time()

			print( f"Total entries so far = {len( mobidb['Uniprot ID'] )}\n" )
			print( f"Time taken = {( toc - tic )/60} minutes" )

		df = pd.DataFrame()
		for key in mobidb:
			print( f"{key} \t {len( mobidb[key] )}" )
			df[key] = mobidb[key]

		df.to_csv( self.mobidb_output_file, index = False )
		tom = time.time()

		self.logger["counts"]["mobidb"] = len( df["Uniprot ID"] )
		self.logger["counts"]["mobidb_unique"] = len( set( df["Uniprot ID"] ) )
		self.logger["time_taken"]["mobidb"] = ( tom - tim )

		return df


	def get_PDBs_for_uniprot_id( self, uni_id ):
		"""
		Obtain PDB IDs corresponding to a set of Uniprot IDs. 
		All PDB IDs related to each Uniprot ID are obtained using the PDB REST API.
		
		Input:
		----------
		uni_id --> UniProt accession ID.

		Returns:
		----------
		[Uniprot ID, list of PDB IDs]
		"""
		pdbs = get_PDB_from_Uniprot_pdb_api( uni_id, 
											max_trials = self.max_trials, 
											wait_time = self.wait_time )

		return [uni_id, pdbs]



	def get_all_unique_pdbs( self, all_uniprot_ids ):
		"""
		Obtain all PDB IDs corresponding to all the Uniprot IDs from DisProt/IDEAL/MobiDB and save as a csv file.
		"""

		if os.path.exists( self.output_file ):
			print( f"{self.output_file} already exists..." )
			print( f"Total PDBs obtained = {self.logger['counts']['total_pdbs']}" )
			print( f"Total unique PDBs obtained so far = {self.logger['counts']['unique_pdbs']}\n" )
			print( f"Time taken to get all unique PDB IDs = {self.logger['time_taken']['unique_pdbs']/3600} hours" )
		
		else:
			tic = time.time()
			self.logger["no_PDBs"] = []
			pdb_ids = []
			
			counter = 0
			with Pool( self.max_cores ) as p:
				for result in tqdm.tqdm( 
										p.imap_unordered( 
														self.get_PDBs_for_uniprot_id, 
														all_uniprot_ids ), 
										total = len( all_uniprot_ids ) ):

					if len( result[1] ) == 0:
						self.logger["no_PDBs"].append( result[0] )
					else:
						# Multiple PDBs can exist for a Uniprot ID.
						pdb_ids.extend( result[1] )

					if counter != 0 and counter%50000 == 0:
						print( f"Total PDBs obtained so far = {len( pdb_ids )} \t {len( pd.unique( pdb_ids ) )}" )

					counter += 1

			unique_pdbs = list( pd.unique( pdb_ids ) )

			self.logger["counts"]["total_pdbs"] = len( pdb_ids )
			self.logger["counts"]["unique_pdbs"] = len( unique_pdbs )
			with open( self.output_file, "w" ) as w:
				w.writelines( ",".join( unique_pdbs ) )

			toc = time.time()
			self.logger["time_taken"]["unique_pdbs"] = ( toc - tic )
			
			print( f"Total PDBs obtained = {len( pdb_ids )}" )
			print( f"Unique PDBs obtained so far = {len( unique_pdbs )}\n" )
			print( f"Time taken to get all unique PDB IDs = {self.logger['time_taken']['unique_pdbs']/3600} hours" )



	def save_logs( self ):
		"""
		Save logs on disk.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		proc = subprocess.Popen( "hostname", shell = True, stdout = subprocess.PIPE )
		system = proc.communicate()[0]
		proc = subprocess.Popen( "date", shell = True, stdout = subprocess.PIPE )
		sys_date = proc.communicate()[0]

		self.tom = time.time()
		self.logger["time_taken"]["total"] = ( self.tom - self.tim )
		with open( f"{self.logger_file}.json", "w" ) as w:
			json.dump( self.logger, w )

		with open( f"{self.logger_file}.txt", "w" ) as w:
			w.writelines( f"Created on: System = {system} \t Date = {sys_date}\n" )

			w.writelines( "\n\n-----------------------------------------------------------\n" )
			w.writelines( "\t\t Counts ------------\n" )
			w.writelines( f"Total Uniprot IDs obtained from DIBS: {self.logger['counts']['dibs']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from MFIB: {self.logger['counts']['mfib']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from FuzDB: {self.logger['counts']['fuzdb']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from PDBtot: {self.logger['counts']['pdb_tot']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from PDBcdr: {self.logger['counts']['pdb_cdr']}\n" )
			w.writelines( f"Total entries obtained from DisProt: {self.logger['counts']['disprot']}\n" )
			w.writelines( f"Unique entries obtained from DisProt: {self.logger['counts']['disprot_unique']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from DisProt: {self.logger['counts']['disprot_uniIDs']}\n" )
			w.writelines( f"Unique Uniprot IDs obtained from DisProt: {self.logger['counts']['disprot_unique_uniIDs']}\n" )
			w.writelines( f"Total entries obtained from IDEAL: {self.logger['counts']['ideal']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from IDEAL: {self.logger['counts']['ideal_uniIDs']}\n" )
			w.writelines( f"Unique Uniprot IDs obtained from IDEAL: {self.logger['counts']['ideal_unique_uniIDs']}\n" )
			w.writelines( f"Total Uniprot IDs obtained from MobiDB: {self.logger['counts']['mobidb']}\n" )
			w.writelines( f"Unique Uniprot IDs obtained from MobiDB: {self.logger['counts']['mobidb_unique']}\n" )
			w.writelines( f"Total Uniprot IDs obtained: {self.logger['counts']['total_uni_ids']}\n" )
			w.writelines( f"Unique Uniprot IDs obtained: {self.logger['counts']['unique_uni_ids']}\n" )
			w.writelines( f"Total PDB IDs obtained: {self.logger['counts']['total_pdbs']}\n" )
			w.writelines( f"Unique PDB IDs obtained: {self.logger['counts']['unique_pdbs']}\n" )

			w.writelines( "\n\n-----------------------------------------------------------\n" )
			w.writelines( "\t\t Time taken ------------\n" )
			w.writelines( f"Time taken - DIBS: {self.logger['time_taken']['dibs']} seconds\n" )
			w.writelines( f"Time taken - MFIB: {self.logger['time_taken']['mfib']} seconds\n" )
			w.writelines( f"Time taken - FuzDB: {self.logger['time_taken']['fuzdb']} seconds\n" )
			w.writelines( f"Time taken - PABtot and PDBcdr: {self.logger['time_taken']['tot_cdr']} seconds\n" )
			w.writelines( f"Time taken - DisProt: {self.logger['time_taken']['disprot']} seconds\n" )
			w.writelines( f"Time taken - IDEAL: {self.logger['time_taken']['ideal']/60} minutes\n" )
			w.writelines( f"Time taken - MobiDB: {self.logger['time_taken']['mobidb']/3600} hours\n" )
			w.writelines( f"Time taken to get Unique Uniprot IDs: {self.logger['time_taken']['unique_uni_ids']/3600} hours\n" )
			w.writelines( f"Time taken to get PDB IDs: {self.logger['time_taken']['unique_pdbs']/3600} hours\n" )
			w.writelines( f"Total Time taken: {self.logger['time_taken']['total']/3600} hours" )

			w.writelines( "\n\n-----------------------------------------------------------\n" )
			w.writelines( f"UniProt IDs with no PDB IDs = {len( self.logger['no_PDBs'] )}\n" )



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Parse a) sequence and b) structure databases to get two corresponding "+
								  "lists of PDBs (Merged_*.csv) containing IDRs in complex.") 

	parser.add_argument( '--max_cores', '-c', dest="c", help = "No. of cores to be used.", type = int, required=False, default=10 )
	
	database_path = f"../Database/"
	dir_name = "input_files/"

	if not os.path.exists( f"{database_path}{dir_name}" ):
		os.makedirs( f"{database_path}{dir_name}" )

	os.chdir( database_path )

	cores = parser.parse_args().c

	# Get the PDB IDs of structures of complexes with IDRs from disorder protein databases.

	TheIlluminati( dir_name, cores ).forward()

	print( "May the Force be with you..." )

