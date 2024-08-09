########## StrIDR (database) mode: Get the PDB structures, sequences, and associated info for the database of IDRs with structures. ##########
########## Disobind (dataset) mode: Download PDBs of IDRs in complexes and parse to obtain binary complexes of IDRs with a partner protein. ##########

########## ------>"May the Force serve u well..." <------###########
####################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import sys
import time
from omegaconf import OmegaConf
import os
import json
import subprocess
import argparse
import tqdm
from multiprocessing import Pool
import random
import h5py

from from_APIs_with_love import ( sifts_map_shell_command, 
								download_SIFTS_Uni_PDB_mapping,
								get_sifts_mapping, 
								from_pdb_rest_api_with_love, pdb_valid,
								download_pdb, get_superseding_pdb_id,
								get_uniprot_entry_name, get_uniprot_seq )

from utility import ( change_basis2, remove_nulls2, 
					load_disorder_dbs,
					find_disorder_regions, 
					get_overlap, check_for_overlap,
					merged_seq_exceeds_maxlen,
					ranges, sort_by_residue_positions 
					)

import warnings
warnings.filterwarnings("ignore")

class parse_pdbs_for_idrs():
	def __init__( self, cores, create_dataset ):
		"""
		Constructor.
		"""
		# Set the seeds for PRNG.
		self.global_seed = 11
		self.seed_worker()
		
		self.cores = cores
		self.max_trials = 25
		self.wait_time = 20
		self.batch_size = 25000
		self.create_dataset = create_dataset

		self.base_path = "/data2/kartik/Disorder_Proteins/disobind/Database/"
		# self.base_path = "../Database/"
		self.merged_PDBs = f"../input_files/Merged_PDB_IDs.txt"
		# Dir for PDB/CIF files.
		self.PDB_path = f"../Combined_PDBs/"
		# Dir for json files of PDB entry and entity info from PDB REST API. 
		self.PDB_api_path = f"../PDB_api/"
		#Dir for SIFTS output files.
		self.mapped_PDB_path = f"../Mapped_PDBs/"
		self.disprot_path = f"../input_files/DisProt.csv"
		self.ideal_path = f"../input_files/IDEAL.csv"
		self.mobidb_path = f"../input_files/MobiDB.csv"

		if self.create_dataset:
			print( "------------------------------- Creating files for Disobind dataset...\n" )
			self.uni_max_len = 10000 if self.create_dataset else None
			self.max_len = 100
			self.min_len = 20
			self.max_num_chains = 50
			self.min_disorder_percent = 0.2
			self.fragment_exceeding_maxlen = True
			self.prot2_frags_limit = None
			self.acceptance_threshold = 0.0

			# PDBs that were downloaded after parallelize_PDB_download stage.
			self.downloaded_pdbs_file = "../Downloaded_PDBs_disobind.txt" 

			# Complexes obtained after running dataset_creation prior to 
			# 	creating all vs all binary complexes.
			self.idr_complexes_dir = "./IDR_PDBs/"
			self.idr_pdbs_file = "./IDR_PDBs.json"
			self.idr_pdbs_selected = {}

			# Directory to store binary complexes for each Uniprot ID pair.
			self.binary_complexes = {}
			self.binary_complexes_dir = f"./Binary_complexes_{self.prot2_frags_limit}/"
			# File to store the Uniprot ID pair keys post creating binary complexes.
			self.binary_complexes_file = f"./Binary_complexes_{self.prot2_frags_limit}.txt" #TODOshru what is this
			
			# File to store the Uniprot ID pair keys.
			self.uniprot_pairs_file = "./Uniprot_pairs.txt"
			self.uniprot_pairs = []

			self.output_csv_path = f"../csv_files_{self.max_len}_{self.min_len}_{self.min_disorder_percent}/"
			self.output_file_name = f"Disobind_dataset"
			self.output_dir_name = f"{self.base_path}Disobind_dataset_{self.max_len}_{self.min_len}_{self.min_disorder_percent}/"

		else:
			print( "------------------------------- Creating files for StrIDR database...\n" )
			self.uni_max_len = None
			self.min_len = 1

			# PDBs that were downloaded after parallelize_PDB_download stage
			self.downloaded_pdbs_file = "../Downloaded_PDBs_stridr.txt" 

			
			self.output_file_name = f"StrIDR_database"
			self.output_dir_name = f"{self.base_path}{self.output_file_name}/"
			self.processed_entry_dir = f"{self.output_dir_name}Processed_entry/"
			self.database_dict = {}

		# File to store Uniprot sequences.
		# Uniprot seq that remain post removing invalid/obsolete ones.
		self.uniprot_seq_file = "./Uniprot_seq.json"
		# Uniprot seq post selecting IDR PDBs.
		self.sel_uniprot_seq_file = "./Disobind_Uniprot_seq.json"
		self.uniprot_seq = {}

		# File to store PDB info obtained after downloading the PDBs.
		self.pdb_info_file = f"Selected_PDBs_info.h5"  # File corresponding to stored pdb_info_df below.
		# Table on each PDB, stores info about PDB ID, auth asym ID, asym ID, entity IDs, Uniprot IDs.
		self.pdb_info_df = [] 

		# File to store PDBs leftover post removing obsolete/invalid Uni IDs.
		self.filtered_pdb_file = f"Filtered_PDB_entries.h5"

		self.logger_file = f"Logs_{self.output_file_name}.json"
		self.logger = {"time_taken": {}, "counts": {}}


	def seed_worker( self ):
		"""
		Set the seed for PRNG.
		"""
		random.seed( self.global_seed )
		np.random.seed( self.global_seed )


	def forward( self ):
		"""
		Disobind: Parse the PDBs containing IDRs to store information on binary complexes of IDRs and a partner protein. 
		StrIDR: Download PDB and Uniprot files for IDRs.
		"""

		"""
		Create the required directories if they don't exist.
		"""
		self.tim = time.time()
		if not os.path.exists( f"{self.output_dir_name}" ):
			os.makedirs( f"{self.output_dir_name}" )
		
		if not os.path.exists( f"{self.output_dir_name}parse_sifts.py" ):
			subprocess.call( ["cp", "./parse_sifts.py", f"{self.output_dir_name}parse_sifts.py"] )
		
		os.chdir( f"{self.output_dir_name}" )
		
		if not os.path.exists( f"{self.mapped_PDB_path}" ):
			os.makedirs( f"{self.mapped_PDB_path}" )
		if not os.path.exists( f"{self.PDB_path}" ):
			os.makedirs( f"{self.PDB_path}" )
		if not os.path.exists( f"{self.PDB_api_path}" ):
			os.makedirs( f"{self.PDB_api_path}" )


		# Load the disordered protein databases.
		self.disprot, self.ideal, self.mobidb = load_disorder_dbs( 
																	self.disprot_path,
																	self.ideal_path,
																	self.mobidb_path )
		# <=======================================================> #
		"""
		Filter PDBs as specified in self.parallelize_PDB_download().
		If complexes already downloaded, use the existing ones else download anew.
		Download only complexes PDB for Disobind dataset.
		Download both monomers and complexes PDB for StrIDR database.
		"""

		if os.path.exists( self.pdb_info_file ):
			# Flag to clear logger keys required downstream.
			print( "Loading state post PDB download...\n" )
			self.pdb_info_df = pd.read_hdf( self.pdb_info_file )
			print( "No. of complexes obtained: ", len( pd.unique( self.pdb_info_df["PDB ID"] ) ) )

			if os.path.exists( self.logger_file ):
				with open( self.logger_file, "r" ) as f:
					self.logger = json.load( f )
				 
		else:
			t1 = time.time()
			# Create keys for logging relevant info to global logger.
			for key in ["defectors",
						"pdb_not_exist", 
						"not_obtained_from_rest_api",
						"deprecated_pdb_id", "chimeric", "np_entity",
						"no_chains_left_in_PDB", "monomer",
						"no_SIFTS_mapping", "pdb_not_downloaded"]:
				self.logger[key] = [[], 0]
			
			self.logger["counts"]["all_chains_in_pdb"] = 0

			self.parallelize_PDB_download()
	
			t2 = time.time()
			self.logger["time_taken"]["PDB_download"] = t2-t1
			print( f"Time taken for filtering PDBs: {( t2-t1 )/60} minutes", "\n" )

			self.logger["counts"]["PDB_download"] = len( set( self.pdb_info_df["PDB ID"].values ) )

			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )

		print( "\n------------------------------------------------------------------\n" )
		# <=======================================================> #
		"""
		Download Uniprot sequences for the selected PDB complexes.
		"""
		if os.path.exists( self.uniprot_seq_file ):
			# Flag to clear all logger keys required downstream of parallelize_PDB_download()
			print( "Loading state post uniprot seq download...\n" )
			with open( self.uniprot_seq_file, "r" ) as f:
				self.uniprot_seq = json.load( f )

			print( "Total Uniprot seq obtained: %s"%len( self.uniprot_seq ) )

			if os.path.exists( self.logger_file ):
				with open( self.logger_file, "r" ) as f:
					self.logger = json.load( f )
		else:
			t1 = time.time()
			for key in ["obsolete_uni_id", f">{self.uni_max_len}"]:
				self.logger[key] = [[], 0]
			self.parallelize_uniprot_seq_download()
	
			t2 = time.time()
			self.logger["time_taken"]["Uniprot_seq_download"] = t2-t1
			print( f"Time taken for downloading Uniprot seq: {( t2-t1 )/60} minutes", "\n" )

			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )

		
		print( "\n------------------------------------------------------------------\n" )
		# <=======================================================> #
		"""
		Remove rows which contain excluded Uniprot IDs.
		"""
		print( "Removing invalid/obsolete Uniprot IDs from PDB entries..." )
		if os.path.exists( self.filtered_pdb_file ):
			self.pdb_info_df = pd.read_hdf( self.filtered_pdb_file )
			x = len( set( self.pdb_info_df["PDB ID"].values ) )
			print( "PDB IDs left post filtering by Uniprot IDs: %s\n"%( x ) )
		
		else:
			self.exclude_problematic_uniprots_from_df()
			
			self.pdb_info_df.to_hdf( self.filtered_pdb_file, key = "data", mode = "w" )
			
			x = len( set( self.pdb_info_df["PDB ID"].values ) )
			print( "PDB IDs left post post filtering by Uniprot IDs: %s\n"%( x ) )
			self.logger["counts"]["PDBs_after_filtering_by_uniprot_ids"] = x

			# Save the state ( logs.
			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )
		print( "\n------------------------------------------------------------------\n" )

		"""
		Here we create either of two things:
		1. PDB and Uniprot files for the StrIDR Database.
		2. Binary complexes of IDRs and a partner protein for Disobind.  
		"""
		# <=======================================================> #
		if self.create_dataset:
			self.dataset_module()
		else:
			self.database_module()


	def dataset_module( self ):
		"""
		Consists of 3 stages:
			1. Identify disordered and ordered chains in PDB - self.parallelize_dataset_creation().
			2. Split non-binary complexes into binary complexes - create_binary_complexes().
			3. Segregate all binary complexes into non-overlapping sets and save - 
												parallelize_nonoverlapping_sets_creation().

		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""

		# <=======================================================> #
		# Load the unsplit complexes if already existing.
		if os.path.exists( self.idr_pdbs_file ):
			print( "Loading saved extracted complexes...\n" )
			# with open( self.dataset_dict_file, "r" ) as f:
			# 	self.dataset_dict = json.load( f )
			with open( self.idr_pdbs_file, "r" ) as f:
				self.idr_pdbs_selected = json.load( f )

			with open( self.sel_uniprot_seq_file, "r" ) as f:
				self.uniprot_seq = json.load( f )

			print( "PDB IDs obtained post obtaining IDR and ordered chains: %s"%( len( self.idr_pdbs_selected ) ) )
			print( "Unique Uniprot IDs left: %s\n"%( len( self.uniprot_seq ) ) )

		else:
			tic = time.time()
			# If re-running post filtration, clear all downstream keys.
			for key in [f">{self.max_num_chains}_chains_in_pdb", "chain_not_mapped", 
						f"#residues<{self.min_len}",
						"mismatch_in_num_uni/pdb_frags", 
						"mismatch_in_uni/pdb_frag_lengths", "aberrant_unIID/mapping",
						 "no_IDRs_in_PDB"]:
				self.logger[key] = [[], 0]
			
			if not os.path.exists( self.idr_complexes_dir ):
				os.makedirs( self.idr_complexes_dir )
			self.parallelize_dataset_creation()
			self.logger["counts"]["IDR_PDBs_obtained"] = len( self.idr_pdbs_selected )
			self.logger["counts"]["Uniprot_IDs_remaining"] = len( self.uniprot_seq )
			
			print( "PDB IDs obtained post obtaining IDR and ordered chains: %s"%( len( self.idr_pdbs_selected ) ) )
			print( "Unique Uniprot IDs left: %s\n"%( len( self.uniprot_seq ) ) )
			
			toc = time.time()

			self.logger["time_taken"]["IDR_PDBs_obtained"] = toc-tic
			print( f"Time taken for extracting IDR entries = {( toc-tic )/60} minutes...")

			# Save the state ( logs and extracted complexes ) on disk.
			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )

		del self.pdb_info_df
		print( "\n------------------------------------------------------------------\n" )
		# <=======================================================> #
		# Load the split complexes if already existing.
		if os.path.exists( self.binary_complexes_file ):
			print( "Using split binary complexes...\n" )
			with open( self.binary_complexes_file, "r" ) as f:
				# binary_complexes = json.load( f )
				binary_complexes = f.readlines()[0].split( "," )
			print( "Total Uniprot ID pairs: ", self.logger["counts"]["uniprot_ID_pairs"] )
			print( "Total binary complexes: ", self.logger["counts"]["total_binary_complexes"] )

		else:
			tic = time.time()
			# Split non-binary complexes into binary complexes.
			print( "Splitting non-binary complexes to binary complexes...\n" )
			if not os.path.exists( self.binary_complexes_dir ):
				os.makedirs( self.binary_complexes_dir )
			self.create_binary_complexes()

			self.logger["counts"]["uniprot_ID_pairs"] = len( self.binary_complexes )
			print( "Total unique Uniprot ID pairs: ", len( self.binary_complexes ) )
			print( "Total binary complexes: ", self.logger["counts"]["total_binary_complexes"] )
			
			del self.binary_complexes
			toc = time.time()
			self.logger["time_taken"]["binary_complexes"] = toc-tic
			print( f"Time taken for creating binary complexes = {( toc-tic )/60} minutes...\n")

			# Save the state ( logs and extracted complexes ) on disk.
			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )

		print( "\n------------------------------------------------------------------\n" )

		self.save_logs()


	def database_module( self ):
		"""
		Consists of 2 stages:
			1. Identify disordered chains in PDB - self.database_creation_module().
			2. Save collected data - self.save_data().

		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		if os.path.exists( f"{self.output_file_name}.json" ):
			print( "Loading saved extracted complexes...\n" )
			with open( f"{self.output_file_name}.json", "r" ) as f:
				self.database_dict = json.load( f )

			print( "Total PDBs included: %s\n"%( len( self.database_dict ) ) )
		
		else:
			tic = time.time()
			for key in ["chain_not_mapped", "no_disordered_regions", "entities_IDR_in_seq", 
						"chains_IDR_in_seq", "entities_IDR_in_struct", "chains_IDR_in_struct", 
						"no_IDRs_in_PDB"]:
				self.logger[key] = [[], 0]
			self.logger["counts"]["pdbs_IDR_in_seq"] = 0
			self.logger["counts"]["all_chains_post_download"] = 0
			self.logger["counts"]["all_entities_post_download"] = 0
			self.logger["counts"]["uni_ids_with_IDR_in_struct"] = 0

			self.parallelize_database_creation()
			
			toc = time.time()
			self.logger["time_taken"]["IDR_PDBs_obtained"] = toc-tic

			print( f"Time taken for creating the database = {( toc-tic )/60} minutes...")

			self.logger["counts"]["IDR_PDBs_obtained"] = len( self.database_dict )

			with open( self.logger_file, "w" ) as w:
				json.dump( self.logger, w )

		print( "\nCreating plots..." )
		self.get_disorder_count_for_entries()
		self.plot_expt_methods_count()
	
		self.save_logs()


	def exclude_problematic_uniprots_from_df( self ):
		"""
		Remove rows from the self.pdb_info_df which contain:
			obsolete uni id.
			invalid uni id.
			uni_id exceeding max Uniprot seq length.

		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		self.pdb_info_df = self.pdb_info_df.reset_index( drop = True )
		uni_ids_to_remove = []
		uni_ids = self.pdb_info_df["Uniprot ID"].values
		uni_ids_to_remove = self.logger["obsolete_uni_id"][0] + self.logger[f">{self.uni_max_len}"][0]
		# self.logger["invalid_uni_id"][0]
		# uni_ids_to_remove += self.logger[f">{self.uni_max_len}"][0]
		# uni_ids_to_remove = [uni_id for uni_id in uni_ids if uni_id not in self.uniprot_seq.keys()]
			
		for i in self.pdb_info_df.index:
			uni_ids_to_keep = list( set( self.pdb_info_df.loc[i, "Uniprot ID"].split( "," )  ) - set( uni_ids_to_remove ) )
			
			# Drop the row if all Uniprot IDs are obsolete.
			if len( uni_ids_to_keep ) == 0:
				self.pdb_info_df = self.pdb_info_df.drop( i )
			
			# Keep only non-obsolete Uniprot IDs.
			else:
				self.pdb_info_df["Uniprot ID"][i] = ",".join( uni_ids_to_keep )



	def remove_missing_residues( self, mapped_uni_pos, mapped_pdb_pos ):
		"""
		Missing residues do not have a structure in the PDB file. 
		We remove them because:
			May pose a problem while merging (for considering heterogeneity) downstream.
				They may contribute to FPs and FNs while training the model.
		In SIFTS mapping, missing residues appear as "nulls" or nan.

		e.g. 
		Indices:: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] --> [[0, 1, 2, 3], [6, 7, 8, 9]]
		PDB pos:: [1, 2, 3, 4, "null", "null", 7, 8, 9, 10] --> [[1, 2, 3, 4], [7, 8, 9, 10]]
		Use index list to remove corresponding residue positions from mapped_uni_pos.
		
		Input:
		----------
		mapped_uni_pos --> list of Uniprot residue positions from SIFTS PDB-Uniprot mapping.
		mapped_pdb_pos --> list of PDB residue positions from SIFTS PDB-Uniprot mapping.

		Returns:
		----------
		uni_pos_frags --> list of lists (Uniprot residue positions) with missing residue positions removed.
		pdb_pos_frags --> list of lists (PDB residue positions) with missing residue positions removed.
		"""
		indices = np.arange( 0, len( mapped_pdb_pos ), 1 )
		pdb_pos_no_nulls, indices = remove_nulls2( mapped_pdb_pos, indices )

		mapped_pdb_pos_uni_basis = []
		uni_pos_frags, pdb_pos_frags = [], []
		for idx in range( len( pdb_pos_no_nulls ) ):
			if len( pdb_pos_no_nulls[idx] ) < self.min_len:
				uni_pos_frags, pdb_pos_frags = [], []
			
			else:
				# Consider residue in Uniprot mapping only if they are in PDB residues also - sanity check.
				uni_pos_frags.append( [mapped_uni_pos[ind] for ind in indices[idx]] )
				uni_pos_frags[-1] = list( map( int,  uni_pos_frags[-1] ) )

				pdb_pos_frags.append( list( map( int,  pdb_pos_no_nulls[idx] ) ) )


		return uni_pos_frags, pdb_pos_frags


	def split_frag( self, fragment ):
		"""
		Split fragments into half recursively, until all are < self.max_len.
		
		Input:
		----------
		fragment --> list of residue positions.

		Returns:
		----------
		fragment --> list of residue positions.
		"""
		if len( fragment ) <= self.max_len:
			return [fragment]
		else:
			shatter_point = len( fragment )// 2
			chunk1 = fragment[:shatter_point]
			chunk2 = fragment[shatter_point:]
			return self.split_frag( chunk1 ) + self.split_frag( chunk2 )


	def fragment_chain( self, plist ):
		"""
		Loop over all the fragements and further break into half until
				length of each fragment < self.max_len.
		Output a new_list containing fragments < self.max_len.
		plist --> list of fragments.
			Each fragment is a list of residue positions.

		Input:
		----------
		plist --> list of lists. Each fragment is a list of residue positions.

		Returns:
		----------
		new_list --> list of lists with all fragments of length < self.max_len.
		"""
		new_list = []
		for region in plist:
			new_list.extend( self.split_frag( region ) )
		return new_list


	def create_key( self, uni_id1, uni_id2 ):
		key = f"{uni_id1}--{uni_id2}"
		
		# print( key, "\t", uni_id1, "\t", uni_id2 )
		if key not in self.binary_complexes.keys():
			# print( key, "\t", uni_id1, "\t", uni_id2, "\t------------------------" )
			self.binary_complexes[key] = {k:[] for k in ["PDB ID", "Auth Asym ID1", "Auth Asym ID2",
															"Asym ID1", "Asym ID2", 
															"Uniprot positions1", "Uniprot positions2",
															"PDB positions1", "PDB positions2"
															] }
		return key


	def download_pdb_and_sift( self, pdb ):
		"""
		Obtain all entity IDs, asym IDs, auth asym IDs, Uniprot IDs from PDB REST API.
		Replace a PDB ID with a superseeding ID, in case a superseding ID exists.

		Remove an entry if:
			PDB does not exist.
			Monomer - monomers only for StrIDR dataset creation.
			No SIFTS mapping.
		
		Filter chains that are not protein and chimeric chains. 
		Download and save the SIFTS mapping as a .tsv file and PDB as a .pdb or .cif file.
		For errors encountered during downloading, try until max_trials are exhausted.
		
		Input:
		----------
		pdb --> PDB ID.

		Returns:
		----------
		logs_dict --> dict for logging exceptions (See parallelize_pdb_download).
		pdb --> PDB ID (same as input).
		pdb_df --> dataframe containing the following info for each PDB chain
					PDB ID
					Entity IDs
					Asym IDs
					Auth Asym IDs
					UniProt IDs
		"""

		# For uniformity, converting all PDB IDs to lowercase.
		pdb_id = pdb.lower().strip()
		
		# For logs.
		# not_exist, deprecated_pdb, not_mapped, monomer = None, None, None, None
		logs_dict = {key:[] for key in ["pdb_not_exist", "deprecated_pdb_id", 
										"not_obtained_from_rest_api",
										"chimeric", "np_entity", "no_SIFTS_mapping", "monomer",
										"no_chains_left_in_PDB", "pdb_not_downloaded"]}
		logs_dict["all_chains_in_pdb"] = 0
		logs_dict["all_uni_ids_in_pdb"] = []


		for trial in range( self.max_trials ):
			# Randomly delay each trial by a few seconds.
			rn = np.random.uniform( 1, self.wait_time )
			time.sleep( rn )
			
			try:
				# Check if the PDB ID has been superseded by a new one.
				# Use the new PDB ID if True.
				pdb = get_superseding_pdb_id( pdb_id, max_trials = self.max_trials, wait_time = self.wait_time )
				pdb = pdb.lower()

				if pdb != pdb_id:
					# deprecated_pdb = f"{pdb_id}--{pdb}"
					logs_dict["deprecated_pdb_id"].append( f"{pdb_id}--{pdb}" )
				
				# If PDB data dict is already saved, load and read it.
				if os.path.exists( f"{self.PDB_api_path}{pdb}_entry.json" ):
					with open( f"{self.PDB_api_path}{pdb}_entry.json" ) as f:
						entry_data = json.load( f )
					with open( f"{self.PDB_api_path}{pdb}_entity.json" ) as f:
						entity_data = json.load( f )
				else:
					entry_data, entity_data = None, None
				
				# Get entity and entry info from pdb. 
				pdb_info = from_pdb_rest_api_with_love( pdb, 
														max_trials = self.max_trials, 
														wait_time =  self.wait_time,
														custom = [entry_data, entity_data] )


				# PDB ID does not exist.
				if pdb_info == None: #TODO above function returns a tuples of Nones. Check accordingly --Done(look from_pdb_api func)
					if trial != self.max_trials:
						continue
					else:
						logs_dict["pdb_not_exist"].append( pdb )
						pdb, pdb_df = None, None

				else:
					pdb_df, entry_data, entity_data, chimeric, np_entity, total_chains, all_uni_ids = pdb_info

					if entry_data == None or entity_data == None:
						if os.path.exists( f"{self.PDB_api_path}{pdb}_entry.json" ):
							subprocess.call( ["rm", f"{self.PDB_api_path}{pdb}_entry.json"] )
						elif os.path.exists( f"{self.PDB_api_path}{pdb}_entity.json" ):
							subprocess.call( ["rm", f"{self.PDB_api_path}{pdb}_entity.json"] )
						
						if trial == self.max_trials - 1:
							print( "Couldn't fetch info from PDB REST API: ", pdb )
							logs_dict["not_obtained_from_rest_api"] = pdb
							pdb, pdb_df = None, None

						else:
							continue

					# Not considering chimeric entities.
					logs_dict["chimeric"].extend( chimeric )
					# Non protein (np) entity.
					logs_dict["np_entity"].extend( np_entity )
					# All protein and non-protein chains in the PDB entry.
					logs_dict["all_chains_in_pdb"] = total_chains
					logs_dict["all_uni_ids_in_pdb"].extend( all_uni_ids )

					if len( pdb_df ) == 0:
						logs_dict["no_chains_left_in_PDB"].append( pdb )
						pdb, pdb_df = None, None

					else:					
						# Check if PDB is a monomer or not.
						# 	Monomer PDB has a single entity with a single chain.
						# 		only 1 row and 1 chain present.
						# 	We also count PDBs with only protein chain left as a monomer.
						# Removing monomers only for the dataset creation.
						if len( pdb_df ) == 1 and self.create_dataset:
							logs_dict["monomer"].append( pdb )
							pdb, pdb_df = None, None

						else:
							mapping_success = download_SIFTS_Uni_PDB_mapping( mapped_PDB_path = self.mapped_PDB_path, 
																				pdb = pdb, 
																				max_trials = self.max_trials, 
																				wait_time = self.wait_time )
							# TODO where are the cases? discrepancy in 50294 - 48755 --Done: They were not downloaded.
							if mapping_success:
								if os.path.exists( f"{self.PDB_path}{pdb}.pdb" ):
									if pdb_valid( f"{self.PDB_path}{pdb}.pdb", "pdb" ):
										pass
									else:
										subprocess.call( ["rm", f"{self.PDB_path}{pdb}.pdb"] )
										continue

								elif os.path.exists( f"{self.PDB_path}{pdb}.cif" ):
									if pdb_valid( f"{self.PDB_path}{pdb}.pdb", "cif" ):
										pass
									else:
										subprocess.call( ["rm", f"{self.PDB_path}{pdb}.cif"] )
										continue
								
								else:
									# Download .pdb or .cif for the PDB id not already existing.
									ext = download_pdb( pdb, max_trials = self.max_trials, 
														wait_time = self.wait_time, 
														return_id = False )
									
									if ext == None:
										logs_dict["pdb_not_downloaded"].append( pdb )
										print( "Not downloaded: ", pdb )
										pdb, pdb_df = None, None
									
									else:
										os.rename( f"./{pdb}.{ext}", f"{self.PDB_path}{pdb}.{ext}" )

										if not os.path.exists( f"{self.mapped_PDB_path}{pdb}.tsv" ):
											os.rename( f"./{pdb}.tsv", f"{self.mapped_PDB_path}{pdb}.tsv" )
									
										pdb_entry_path = f"{self.PDB_api_path}{pdb}_entry.json"
										pdb_entity_path = f"{self.PDB_api_path}{pdb}_entity.json"

										if not os.path.exists( pdb_entry_path ):
											out_file1 = open( pdb_entry_path, "w" )
											json.dump( entry_data, out_file1 )
											out_file1.close()

										if not os.path.exists( pdb_entity_path ):
											out_file2 = open( pdb_entity_path, "w" )
											json.dump( entity_data, out_file2 )
											out_file2.close()

							# If no SIFTS mapping obtained.
							else:
								logs_dict["no_SIFTS_mapping"].append( pdb )
								pdb, pdb_df = None, None

				output = [logs_dict, pdb, pdb_df]
				break

			except Exception as e:
				if trial == 0:
					print( f"Trial {trial}: Exception {e} \t --> {pdb}" )
				
				time.sleep( self.wait_time )
				if trial != self.max_trials-1:
					continue
				else:
					print( f"Trial {trial}: Exception {e} \t --> {pdb}" )
					output = [0, pdb, None]
		
		return output


	def fetch_uniprot_seq_and_name( self, uni_id ):
		"""
		Download Uniprot sequence for the given Uniprot ID.
		Log exceptions (invalid, obsolete, too long IDs). 

		Input:
		----------
		uni_id --> Uniprot accession ID.

		Returns:
		----------
		uni_id --> Uniprot accession ID (same as input).
		uni_seq --> sequence for the given uni_id.
		logs_dict --> dict for logging exceptions (See parallelize_uniprot_seq_download).
		"""	
		
		logs_dict = {key:[] for key in ["obsolete_uni_id", f">{self.uni_max_len}"]}

		uni_seq = get_uniprot_seq( uni_id = uni_id, 
									max_trials = self.max_trials, 
									wait_time = self.wait_time, 
									return_id = False )

		# Uniprot ID is obsolete.
		if uni_seq == []:
			logs_dict["obsolete_uni_id"].append( uni_id )
			name = None
		else:
			# # Fetch the name for the Uniprot ID if it is not obsolete.
			name = get_uniprot_entry_name( uni_id = uni_id, 
											max_trials = self.max_trials, 
											wait_time = self.wait_time )


		# Remove the chain if the Uniprot seq for the corresponding Uniprot ID exceeds max len.
		if self.uni_max_len != None:
			if len( uni_seq ) > self.uni_max_len:
				logs_dict[f">{self.uni_max_len}"].append( uni_id )
		
		return [uni_id, uni_seq, name, logs_dict]



	def the_forge( self, data ):
		"""
		Obtain info. for creating the database.
		Data obtained from PDB - entity IDs, asym IDs, auth asym IDs, Uniprot IDs.
		Remove entity or chain (see below) if:
				Chain not mapped.
				Mismatch in the no. of mapped PDB and Uniprot fragments.
				Mismatch in the length of mapped PDB and Uniprot fragments.
		Identify IDRs for each chain in each entity in the PDB.
			Obtain all disordered residue using mapped uniprot positions.
			Obtain disordered residues present in structure after removing missing residues.
		Not considering chains with -
			No disordered residues.
			No disordered residues in the structure.
		Remove the entry if no disordered residues in any chain.
		
		Input:
		----------
		data --> pandas groupby object (grouped by PDB ID).

		Returns:
		----------
		logs_dict --> dict for logging exceptions (see parallelize_database_creation()).
		entry_dict --> nested dict.
						{PDB:
							"Uniprot ID":{
									"Name"
									"all Uniprot ID"
									"Sequence"
									"Asym ID"
									"Auth Asym ID"
									"Disordered residues in seq"
									"Disordered residues in struct"
							}

						}
		"""
		logs_dict = { i:[] for i in ["no_disordered_regions", "chain_not_mapped",
										"entities_IDR_in_seq", "chains_IDR_in_seq",
										"entities_IDR_in_struct", "chains_IDR_in_struct"
									] }
		logs_dict["pdbs_IDR_in_seq"] = 0
		logs_dict["all_chains_post_download"] = 0
		logs_dict["all_entities_post_download"] = 0
		logs_dict["uni_ids_with_IDR_in_struct"] = []

		pdb = data[0]
		all_entity_ids = data[1]["Entity ID"].values  
		all_asym_ids = data[1]["Asym ID"].values
		all_auth_asym_ids = data[1]["Auth Asym ID"].values
		all_uniprot_ids = data[1]["Uniprot ID"].values

		if not os.path.exists( f"{self.processed_entry_dir}{pdb}.json" ):
			logs_dict["all_chains_post_download"] += len( all_asym_ids )
			logs_dict["all_entities_post_download"] += len( set( all_entity_ids ) )

			entry_dict = {}
			entry_dict[pdb] = {}

			# <=======================================================> #
			# Get the mapping for the PDB.
			mapped_pdb_file = f"{self.mapped_PDB_path}{pdb.lower()}.tsv"
			mapping = pd.read_csv( mapped_pdb_file, sep = "\t", header = None )
			
			# For all entity IDs.
			# for entity_id in entity_ids:
			for i in range( len( all_entity_ids ) ):
				entity_id = all_entity_ids[i]
				asym_id = all_asym_ids[i]
				auth_asym_id = all_auth_asym_ids[i]
				uniprot_ids = list( set( all_uniprot_ids[i].split( "," ) ) )
				tmp = [id_ for id_ in uniprot_ids if id_ in self.uniprot_seq.keys()] 
				uni_id = tmp[0]

				# <=======================================================> #
				#### Get the mapped Uniprot and PDB seq for the chain.
				mapping_dict, _ = get_sifts_mapping( mapping = mapping, 
													chain1 = auth_asym_id, 
													uni_id1 = uniprot_ids )

				# Remove the chain if the mapping does not contain the apt chain ID.
				# 	We only consider mapping for the PDB specified Uniprot ID -- at SIFTS level.
				if mapping_dict["uni_pos"] == []:
					logs_dict["chain_not_mapped"].append( f"{pdb}_{auth_asym_id}" )
					# not_mapped.append( f"{pdb}_{auth_asym_ids[entity_id][j]}" )
					continue
				#TODO just a comment. Not deleted from the row in pdb_info_df if SIFTS mapping does not exist 

				mapped_uni_pos = mapping_dict["uni_pos"]
				mapped_pdb_pos = mapping_dict["pdb_pos"]

				# <=======================================================> #
				uni_pos_frags, pdb_pos_frags = self.remove_missing_residues( mapped_uni_pos, mapped_pdb_pos )

				# If <= self.min_len residues left in the fragment.
				if pdb_pos_frags == []:
					continue

				"""
				# Entries like '5aps' have been mapped to the same Uniprot sequence along the length twice.
				"""
				# <=======================================================> #
				## Get all disordered regions from DisProt/IDEAL.
				disorder_regions, ref = find_disorder_regions( disprot = self.disprot, 
																ideal = self.ideal, 
																mobidb = self.mobidb, 
																uni_ids = uniprot_ids, 
																min_len = 1,
																return_ids = True )
				disorder_in_seq = []
				disorder_in_struct = []
				
				if disorder_regions == []:
					logs_dict["no_disordered_regions"].append( uni_id )
					continue
				
				overlap_residues = get_overlap( mapped_uni_pos, disorder_regions )
				for overlap in overlap_residues:
					disorder_in_seq += overlap
				
				# <=======================================================> #
				# Disorder residues present in Uniprot sequence.
				if disorder_in_seq != []:
					logs_dict["pdbs_IDR_in_seq"] = 1  
					logs_dict["chains_IDR_in_seq"].append( asym_id )
					logs_dict["uni_ids_with_IDR_in_struct"].extend( uniprot_ids )
					
					if entity_id not in logs_dict["entities_IDR_in_seq"]:
						logs_dict["entities_IDR_in_seq"].append( entity_id )
					for frag in uni_pos_frags:
						overlap_residues = get_overlap( frag, disorder_regions )
						for overlap in overlap_residues:
							disorder_in_struct += overlap

					# <=======================================================> #
					# Disorder residues present in structure.
					if disorder_in_struct != []:
						logs_dict["chains_IDR_in_struct"].append( asym_id )
						if entity_id not in logs_dict["entities_IDR_in_struct"]:
							logs_dict["entities_IDR_in_struct"].append( entity_id )

						if uni_id not in entry_dict[pdb].keys():
							entry_dict[pdb][uni_id] = {
									key: [] for key in ["Name", "all_Uniprot_IDs", 
														"Sequence", "Asym_ID", "Auth_Asym_ID", 
														"Disorder_in_seq", "Disorder_in_struct",
														"Count_disorder_residues"]
							}

						entry_dict[pdb][uni_id]["Count_disorder_residues"].append( len( disorder_in_struct ) )
						disorder_in_seq = ranges( disorder_in_seq )
						disorder_in_struct = ranges( disorder_in_struct )

						entry_dict[pdb][uni_id]["Name"].append( self.uniprot_seq[uni_id][1] )
						entry_dict[pdb][uni_id]["all_Uniprot_IDs"].append( uniprot_ids )
						entry_dict[pdb][uni_id]["Sequence"].append( self.uniprot_seq[uni_id][0] )
						entry_dict[pdb][uni_id]["Asym_ID"].append( asym_id )
						entry_dict[pdb][uni_id]["Auth_Asym_ID"].append( auth_asym_id )
						entry_dict[pdb][uni_id]["Disorder_in_seq"].append( disorder_in_seq )
						entry_dict[pdb][uni_id]["Disorder_in_struct"].append( disorder_in_struct )

						entry_dict[pdb][uni_id]["cross_refs"] = {}
						
						disprot_ids, ideal_ids, mobi_ids = ref.split( "--" )
						for dbname, db in zip( ["DisProt", "IDEAL", "MobiDB"], [disprot_ids, ideal_ids, mobi_ids] ):
							entry_dict[pdb][uni_id]["cross_refs"][dbname] = db.split( "," )

			processed_entry = {"entry": entry_dict, "logs": logs_dict}
			with open( f"{self.processed_entry_dir}{pdb}.json", "w" ) as w:
				json.dump( processed_entry, w )

		else:
			with open( f"{self.processed_entry_dir}{pdb}.json", "r" ) as f:
				processed_entry = json.load( f )
				entry_dict = processed_entry["entry"]
				logs_dict = processed_entry["logs"]

		return [logs_dict, entry_dict]



	def where_the_magic_happens( self, data ):
		"""
		Data obtained from PDB - entity IDs, asym IDs, auth asym IDs, Uniprot IDs.
		Remove entity or chain (see below) if:
				Chain not mapped.
				Mismatch in the no. of mapped PDB and Uniprot fragments.
				Mismatch in the length of mapped PDB and Uniprot fragments.
		Identify IDRs for each chain in each entity in the PDB.
			Remove missing residues form Uniprot and PDB mappings --> this creates fragments.
			Remove fragments < min len.
			If lengtn of fragments > 250(max len), further break them into halves.
			Identify IDRs in each fragment by overlapping the fragemnt with disordered regions 
				from DisProt/IDEAL/MobiDB --> prot1(disordered).
				A fragment is considered as IDR, if it has >= self.min_disorder_percent disordered residues.
					Else considered as prot2(not disordered).
			Fragments with no overlap are considered as prot2(not disordered).
		Remove the entry if no IDRs found.
		
		We consider all Uniprot IDs for all chains in an entity.
		Some entity may have >1 Uniprot ID.
			These correspond to the same protein, since we already removed chimeras.
			We assume the 1st one to be representative, while we use all of these 
					to obtain mapping and disordered regions.
		Input:
		----------
		data --> pandas groupby object (grouped by PDB ID).

		Returns:
		----------
		logs_dict --> dict for logging exceptions (see parallelize_dataset_creation()).
		cargo --> list of lists.
				[pdb, prot1, prot2]
					pdb --> PDB ID
					prot1/prot2 --> 
									Asym/Auth_asym IDs.
									Uniprot IDs.
									Uniprot positions for all fragments per PDB.
									PDB positions for all fragments per PDB.
		"""
		logs_dict = { i:[] for i in [f">{self.max_num_chains}_chains_in_pdb", "chain_not_mapped", f"#residues<{self.min_len}", 
								"mismatch_in_num_uni/pdb_frags", 
								"mismatch_in_uni/pdb_frag_lengths", "aberrant_unIID/mapping"] }
		
		pdb = data[0]
		all_entity_ids = data[1]["Entity ID"].values
		all_asym_ids = data[1]["Asym ID"].values
		all_auth_asym_ids = data[1]["Auth Asym ID"].values
		all_uniprot_ids = data[1]["Uniprot ID"].values

		
		# Process the entry if not already processed. An entry has been processed if the corresponding JSON file exists.
		if not os.path.exists( f"{self.idr_complexes_dir}{pdb}.json" ):
			entry_dict = {pdb:{}}
			for key1 in ["prot1","prot2"]:
				entry_dict[pdb][key1] = {}
				for key2 in ["Asym ID", "Auth Asym ID", "Uniprot ID", "Uniprot positions", "PDB positions"]:
					entry_dict[pdb][key1][key2] = []

			# Ignore PDBs with no. of chains >= self.max_num_chains.
			if len( all_asym_ids ) >= self.max_num_chains:
				logs_dict[f">{self.max_num_chains}_chains_in_pdb"].append( pdb )

			else:
				# <=======================================================> #
				# Get the mapping for the PDB.
				mapped_pdb_file = f"{self.mapped_PDB_path}{pdb.lower()}.tsv"
				mapping = pd.read_csv( mapped_pdb_file, sep = "\t", header = None )

				chain_ids_p1, chain_ids_p2 = [[], []], [[], []]
				uni_ids_p1, uni_ids_p2 = [], []
				uni_pos_p1, uni_pos_p2 = [], []
				pdb_pos_p1, pdb_pos_p2 = [], []

				# <=======================================================> #
				# For all entity IDs.
				for i in range( len( all_entity_ids ) ):
					s1 = time.time()
					asym_id = all_asym_ids[i]
					auth_asym_id = all_auth_asym_ids[i]

					# Representative Uni ID must not be obsolete/invalid.
					uniprot_ids = set( all_uniprot_ids[i].split( "," ) )
					tmp = [id_ for id_ in uniprot_ids if id_ in self.uniprot_seq.keys()]
					uni_id = tmp[0]

					# <=======================================================> #
					# Get the mapped Uniprot and PDB seq for the chain.
					mapping_dict, _ = get_sifts_mapping( mapping = mapping, 
														chain1 = auth_asym_id, 
														uni_id1 = uniprot_ids )

					# Remove the chain if the mapping does not contain the apt chain ID.
					# 	We only consider mapping for the PDB specified Uniprot ID -- at SIFTS level.
					if mapping_dict["uni_pos"] == []:
						logs_dict["chain_not_mapped"].append( f"{pdb}_{auth_asym_id}__{uni_id}" )
						continue

					mapped_uni_pos = mapping_dict["uni_pos"]
					mapped_pdb_pos = mapping_dict["pdb_pos"]

					# <=======================================================> #
					## Remove "nulls"/nan.
					## Get the mapped Uniprot and PDB fragments.
					uni_pos_frags, pdb_pos_frags = self.remove_missing_residues( mapped_uni_pos, mapped_pdb_pos )

					# If < self.min_len residues left in the fragment.
					if len( pdb_pos_frags ) == 0:
						self.logger[f"#residues<{self.min_len}"].append( f"{pdb}_{asym_id}" )
						continue

					if self.fragment_exceeding_maxlen:
						# If fragment length > max_len, break it into 2 halves.
						uni_pos_frags = self.fragment_chain( uni_pos_frags )
						pdb_pos_frags = self.fragment_chain( pdb_pos_frags )

					# Ignore if the no. of uni fragments does not match no. of pdb fragments.
					if len( uni_pos_frags ) != len( pdb_pos_frags ):
						if f"{pdb}_{auth_asym_id}" not in logs_dict["mismatch_in_num_uni/pdb_frags"]:
							# mismatch_frags.append( pdb )
							logs_dict["mismatch_in_num_uni/pdb_frags"].append( 
										f"{pdb}_{auth_asym_id}" 
										)
						continue

					# Remove fragments if the length of uniprot and pdb frags do not match.
					# This may have occured due to problems in mapping e.g. 5nrl, 7jij.
					uni_pos_frags_, pdb_pos_frags_ = [], []
					mismatch_found = False
					
					for uni_frag, pdb_frag in zip( uni_pos_frags, pdb_pos_frags ):
						if len( uni_frag ) != len( pdb_frag ):
							mismatch_found = True
							continue
						else:
							uni_pos_frags_.append( uni_frag )
							pdb_pos_frags_.append( pdb_frag )
					
					if mismatch_found:
						uni_pos_frags, pdb_pos_frags = uni_pos_frags_, pdb_pos_frags_
						del uni_pos_frags_
						del pdb_pos_frags_
						# mismatch_frag_lengths.append( pdb )
						logs_dict["mismatch_in_uni/pdb_frag_lengths"].append( 
										f"{pdb}_{auth_asym_id}"
										 )

					# <=======================================================> #
					## Get all disordered regions from DisProt/IDEAL.					
					disorder_regions = find_disorder_regions( disprot = self.disprot, 
																ideal = self.ideal, 
																mobidb = self.mobidb, 
																uni_ids = uniprot_ids, 
																min_len = 1,
																return_ids = False )
					
					# For all the fragments.
					for idx in range( len( uni_pos_frags ) ):
						# Overlap mapped_uni_pos with disorder regions on DisProt/IDEAl.
						# Consider fragments as disordered if atleast min_disorder_percent
						# 		disordered residues are present.
						if len( uni_pos_frags[idx] ) > self.max_len:
							raise Exception( "Fragment exceeds max length..." )

						# Remove fragments if the Uniprot residues are not continous.
						uni_pos = ranges( uni_pos_frags[idx] )
						if len( uni_pos ) > 1:
							continue
						# Remove fragments if there is duplication of residue e.g. 6ap9_B_102,102.
						if len( set( uni_pos_frags[idx] ) ) != len( uni_pos_frags[idx] ):
							continue

						# Remove fragments if Uniprot seq has less residues than in Uniprot mapping.
						# e.g.: G1SJB4
						# 		Uniprot residues present in mapping do not exist in seq.
						# 		Seq length as per Uniprot 317.
						# 		Uniprot residue positions in mapping = 117-424.
						# 		PDB (7jqb) residue positions in mapping = 2-314.
						if uni_pos_frags[idx][-1] > len( self.uniprot_seq[uni_id] ):
							label = f"{pdb}_{auth_asym_id}_{uni_id}"
							if label not in logs_dict["aberrant_unIID/mapping"]:
								logs_dict["aberrant_unIID/mapping"].append( label )
							continue
						
						total_overlap_residues = 0
						for overlap in get_overlap( uni_pos_frags[idx], disorder_regions ):
							total_overlap_residues += len( overlap )

						percent_disorder = total_overlap_residues/len( uni_pos_frags[idx] )
						
						if percent_disorder >= self.min_disorder_percent:
							entry_dict[pdb]["prot1"]["Asym ID"].append( asym_id )
							entry_dict[pdb]["prot1"]["Auth Asym ID"].append( auth_asym_id )
							entry_dict[pdb]["prot1"]["Uniprot ID"].append( uni_id )
							entry_dict[pdb]["prot1"]["Uniprot positions"].append( uni_pos_frags[idx] )
							entry_dict[pdb]["prot1"]["PDB positions"].append( pdb_pos_frags[idx] )

						else:
							entry_dict[pdb]["prot2"]["Asym ID"].append( asym_id )
							entry_dict[pdb]["prot2"]["Auth Asym ID"].append( auth_asym_id )
							entry_dict[pdb]["prot2"]["Uniprot ID"].append( uni_id )
							entry_dict[pdb]["prot2"]["Uniprot positions"].append( uni_pos_frags[idx] )
							entry_dict[pdb]["prot2"]["PDB positions"].append( pdb_pos_frags[idx] )

			processed_entry = {"entry": entry_dict, "logs": logs_dict}
			with open( f"{self.idr_complexes_dir}{pdb}.json", "w" ) as w:
				json.dump( processed_entry, w )

		else:
			with open( f"{self.idr_complexes_dir}{pdb}.json", "r" ) as f:
				processed_entry = json.load( f )
				entry_dict = processed_entry["entry"]
				logs_dict = processed_entry["logs"]


		return [logs_dict, entry_dict]


	def get_values( self, entry ):
		"""
		Obtain the following for the given entry:
			chain IDs, Uniprot IDs, residue positions.

		Input:
		----------
		entry --> dict containing:
					Asym ID
					Auth Asym ID
					Uniprot IDs
					Uniprot positions
					PDB positions

		Returns:
		----------
		asym_ids --> Asym ID for the PDB entry.
		auth_asym_ids --> Auth Asym ID for the PDB entry.
		uni_ids --> Uniprot IDs for the PDB entry.
		uni_pos --> Uniprot positions for the PDB entry.
		pdb_pos --> PDB positions for the PDB entry.
		"""
		asym_ids = entry["Asym ID"]
		auth_asym_ids = entry["Auth Asym ID"]
		uni_ids = entry["Uniprot ID"]
		uni_pos = entry["Uniprot positions"]
		pdb_pos = entry["PDB positions"]

		return asym_ids, auth_asym_ids, uni_ids, uni_pos, pdb_pos


	def create_binary_complexes( self ):
		"""
		Split non-binary complexes or homomers into binary complexes.
		Creating binary complexes for:
			IDR-IDR chains.
			IDR-ordered chains.
		Not considering intrachain frgaments.
		Load the dict containing selected IDR PDBs.
			Loop over and load each PDB entry one at a time.
		eg: A:B --> AB
			A:BC --> AB, AC
			AB:CD --> AC, AD, BC, BD, AB
		
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		for pdb in self.idr_pdbs_selected.keys():
			with open( f"{self.idr_complexes_dir}{pdb}.json", "r" ) as f:
				entry_dict = json.load( f )
				entry_dict = entry_dict["entry"]

			asym_ids1, auth_asym_ids1, uni_ids1, uni_pos1, pdb_pos1 = self.get_values( entry_dict[pdb]["prot1"] )
			asym_ids2, auth_asym_ids2, uni_ids2, uni_pos2, pdb_pos2 = self.get_values( entry_dict[pdb]["prot2"] )

			# Split the IDR chains into binary pairs if prot2 is empty ('', []).
			# ABC: --> AB, AC, BC.			
			if entry_dict[pdb]["prot2"]["Asym ID"] == "":
				for i in range( 1, len( asym_ids1 ) ):

					# Intrachain fragments are not considered.
					if asym_ids1[i-1] != asym_ids1[i]:
						key = self.create_key( uni_ids1[i-1], uni_ids1[i] )

						self.binary_complexes[key]["PDB ID"            ].append( pdb )
						self.binary_complexes[key]["Auth Asym ID1"     ].append( auth_asym_ids1[i-1] )
						self.binary_complexes[key]["Auth Asym ID2"     ].append( auth_asym_ids1[i] )
						self.binary_complexes[key]["Asym ID1"          ].append( asym_ids1[i-1] )
						self.binary_complexes[key]["Asym ID2"          ].append( asym_ids1[i] )
						self.binary_complexes[key]["Uniprot positions1" ].append( uni_pos1[i-1] )
						self.binary_complexes[key]["Uniprot positions2" ].append( uni_pos1[i] )
						self.binary_complexes[key]["PDB positions1" ].append( pdb_pos1[i-1] )
						self.binary_complexes[key]["PDB positions2" ].append( pdb_pos1[i] )

			else:
				for i in range( len( asym_ids1 ) ):
					# Separate the multi-chain IDR (prot1) into binary complexes(ABC:DEF).
					if i != len( asym_ids1 )-1:
						# Following combinations will be produced ABC --> AB, AC, BC
						for j in range( i, len( asym_ids1 ) ):
							if asym_ids1[i] != asym_ids1[j]:

								key = self.create_key( uni_ids1[i], uni_ids1[j] )

								self.binary_complexes[key]["PDB ID"            ].append( pdb )
								self.binary_complexes[key]["Auth Asym ID1"     ].append( auth_asym_ids1[i] )
								self.binary_complexes[key]["Auth Asym ID2"     ].append( auth_asym_ids1[j] )
								self.binary_complexes[key]["Asym ID1"     ].append( asym_ids1[i] )
								self.binary_complexes[key]["Asym ID2"     ].append( asym_ids1[j] )
								self.binary_complexes[key]["Uniprot positions1"].append( uni_pos1[i] )
								self.binary_complexes[key]["Uniprot positions2"].append( uni_pos1[j] )
								self.binary_complexes[key]["PDB positions1"].append( pdb_pos1[i] )
								self.binary_complexes[key]["PDB positions2"].append( pdb_pos1[j] )

					p2_frags = len( asym_ids2 )
					num_chains2 = len( set( auth_asym_ids2 ) )
					counter = {key:0 for key in set( asym_ids2 )}
					# By default accept all fragments.
					acceptance_prob = 1.0
					for k in range( len( asym_ids2 ) ):
						"""
						Separate the multi-chain prot2 entries into separate entries.
						Intrachain fragments are not considered.
						Combinations of prot2 chains not being made.
						Following combinations will be produced A-D,E,F; B-D,E,F; C-D,E,F.
						The Catch:
							(No. of frags explode as frag size decreases or the PDB contains a lot of entities).
							If prot2_frags_limit is specified, don't consider all prot2 fragments.
							Randomly select some fragments from each chain.
						"""
						if asym_ids1[i] != asym_ids2[k]:
							if self.prot2_frags_limit == None:
								pass

							# Randomly choose a subset of prot2 fragments (up to self.prot2_frag_limit) 
							# 	if no. of chains > self.prot2_frag_limit.
							elif p2_frags > self.prot2_frags_limit and num_chains2 > self.prot2_frags_limit:
								# Randomly sample fragments from each chain upto the limit.
								if counter[asym_ids2[k]] > self.prot2_frags_limit:
									continue

								# Randomly sample acceptance probability from a uniform distribution.
								acceptance_prob = np.random.uniform( 0.0, 1.0 )

							# Consider a fragment if its acceptance probability is higher than threshold.
							if acceptance_prob > self.acceptance_threshold:
								key = self.create_key( uni_ids1[i], uni_ids2[k] )

								self.binary_complexes[key]["PDB ID"            ].append( pdb )
								self.binary_complexes[key]["Auth Asym ID1"     ].append( auth_asym_ids1[i] )
								self.binary_complexes[key]["Auth Asym ID2"     ].append( auth_asym_ids2[k] )
								self.binary_complexes[key]["Asym ID1"     ].append( asym_ids1[i] )
								self.binary_complexes[key]["Asym ID2"     ].append( asym_ids2[k] )
								self.binary_complexes[key]["Uniprot positions1"].append( uni_pos1[i] )
								self.binary_complexes[key]["Uniprot positions2"].append( uni_pos2[k] )
								self.binary_complexes[key]["PDB positions1"].append( pdb_pos1[i] )
								self.binary_complexes[key]["PDB positions2"].append( pdb_pos2[k] )
								
								if self.prot2_frags_limit != None:							
									# Keep tab on the no. of fragments included per chain.
									counter[asym_ids2[k]] += 1

		# binary_complexes = []
		total_binary_complexes = 0
		for key in self.binary_complexes.keys():
			total_binary_complexes += len( self.binary_complexes[key]["PDB ID"] )
			with open( f"{self.binary_complexes_dir}{key}.json", "w" ) as w:
				json.dump( self.binary_complexes[key], w )
		
		self.logger["counts"]["total_binary_complexes"] = total_binary_complexes
		with open( self.binary_complexes_file, "w" ) as w:
			w.writelines( ",".join( list( self.binary_complexes.keys() ) ) )



	def parallelize_PDB_download( self ):
		"""
		Filter out unwanted PDBs and download relevant ones before proceeding downstream.
			See self.download_pdb_and_sift() for more.
		Log the following:
			Defectors --> if for any PDB, the Exception could not be resolved.
			If PDB ID does not exist.
			If a deprecated PDB ID  has been replaced with a superseeding one.
			PDB contains DNA/RNA.
			If SIFTS mapping could not be obtained.
			If its a monomeric PDB.
			If no valid (no DNA/RNA) chains remain in PDB.
		
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		self.pdb_info_df = pd.DataFrame()
		with open( self.merged_PDBs, "r" ) as f:
			all_pdbs = f.readlines()[0].split( "," )

		print( "Downloading PDB complexes...\n" )
		if os.path.exists( self.downloaded_pdbs_file ):
			print( "Using pre-selected PDBs..." )
			with open( self.downloaded_pdbs_file, "r" ) as w:
				unique_pdbs_to_download = w.readlines()[0].split( "," )
			self.logger["counts"]["all_PDBs"] = len( all_pdbs )
			self.logger["counts"]["unique_PDBs"] = len( unique_pdbs_to_download )
		else:
			unique_pdbs_to_download = pd.unique( all_pdbs )
			self.logger["counts"]["all_PDBs"] = len( all_pdbs )
			self.logger["counts"]["unique_PDBs"] = len( unique_pdbs_to_download )
		
		print( f"\nTotal unique PDBs = {len( unique_pdbs_to_download )}\n")
		
		all_uni_ids_in_pdb = set()
		# Just to shuffle large PDBs hoarded together.
		np.random.shuffle( unique_pdbs_to_download )
		for start in np.arange( 0, len( unique_pdbs_to_download), self.batch_size ):
			end = start + self.batch_size

			batch = unique_pdbs_to_download[start:end]
			with Pool( self.cores ) as p: 
				for result in tqdm.tqdm( p.imap_unordered( 
													self.download_pdb_and_sift, batch ), 
										total = len( batch ) ):
		
					if result[0] == 0:
						self.logger["defectors"][0].append( result[1] )
						self.logger["defectors"][1] += 1
						continue

					logs_dict = result[0]
					for key in logs_dict.keys():
						if key not in ["all_chains_in_pdb", "all_uni_ids_in_pdb"]:
							self.logger[key][0].extend( logs_dict[key] )
							self.logger[key][1] += len( logs_dict[key] )

					self.logger["counts"]["all_chains_in_pdb"] += logs_dict["all_chains_in_pdb"]
					all_uni_ids_in_pdb.update( logs_dict["all_uni_ids_in_pdb"] )
					
					if result[1] != None:
						result[1] = str( result[1] )
						if len( self.pdb_info_df ) == 0:
							self.pdb_info_df = result[2]
						else:
							self.pdb_info_df = pd.concat( [self.pdb_info_df, result[2]] )
		
		self.logger["counts"]["all_uni_ids_in_pdb"] = len( all_uni_ids_in_pdb )

		result_unique_pdbs = pd.unique( self.pdb_info_df["PDB ID"].values )
		if not os.path.exists( self.downloaded_pdbs_file ):
			with open( f"{self.downloaded_pdbs_file}", "w" ) as w:
				w.writelines( ",".join( result_unique_pdbs ) )

		self.pdb_info_df = self.pdb_info_df.reset_index( drop = True )
		self.pdb_info_df.to_hdf( self.pdb_info_file, key = "data", mode = "w" )

		# Some defectors may occur occasionally due to failure to 
		# 	fetch info. even upon multiple trials.
		# 	Rerun the script, if defectors persist, investigate further.
		print( "Total no. of defectors found: %s \n"%( self.logger["defectors"][1] ) )
		print( "Total PDBs downloaded: %s \n"%( len( result_unique_pdbs ) ) )



	def parallelize_uniprot_seq_download( self ):
		"""
		Parallelize downloading sequences for all unique Uniprot IDs.
		Log the following:
			Invalid Uni IDs.
			Obsolete Uni IDs.
			If the no. of Uni and PDB fragments do not match.
			If for any fragment, the Uni and PDB lengths do not match.
			Uni ID for which seq exceeds max_len.
		
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		print( "Downloading Uniprot sequences...\n" )
		uniprot_ids = pd.unique( ",".join( self.pdb_info_df["Uniprot ID"].values ).split( "," ) )

		with Pool( self.cores ) as p:
			for result in tqdm.tqdm( 
						p.imap_unordered( 
								self.fetch_uniprot_seq_and_name, uniprot_ids ),
						total = len( uniprot_ids )
						 ):

				uni_id, uni_seq, name, logs_dict = result
				for key in logs_dict.keys():
					self.logger[key][0].extend( logs_dict[key] )
					self.logger[key][1] += len( logs_dict[key] )

				if len( uni_seq ) != 0 and name != None:
					if not self.create_dataset:
						self.uniprot_seq[uni_id] = [uni_seq, name]
					else:
						self.uniprot_seq[uni_id] = uni_seq

		print( "Total Uniprot seq obtained: %s \n"%len( self.uniprot_seq ) )
		# print( "Invalid Uniprot IDs: %s"%self.logger["invalid_uni_id"][1] )
		print( "Obsolete Uniprot IDs:  %s"%self.logger["obsolete_uni_id"][1] )
		print( self.logger["obsolete_uni_id"][0] )
		
		with open( self.uniprot_seq_file, "w" ) as w:
			json.dump( self.uniprot_seq, w )



	def parallelize_dataset_creation( self ):
		"""
		Obtain all relevant info. for all selected entries in parallel.
		Log the following:
			Chains that were not mapped.
			Chains that were mapped to incorrect Uni ID.
			Chains with no disordered regions.
			PDBs with no IDRs.
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		print( "Initiating Dataset creation module...\n" )
		counter = 0
		selected_uni_ids = set()
		groups = self.pdb_info_df.groupby( "PDB ID" )
		groups = list( groups )
		
		for start in np.arange( 0, len(  groups ), self.batch_size ):
			end = start + self.batch_size
			batch = groups[start:end]

			with Pool( self.cores ) as p:
				for result in tqdm.tqdm( 
								p.imap_unordered( 
												self.where_the_magic_happens, batch ), 
								total = len( batch ) ):

					logs_dict, entry_dict = result

					for key in logs_dict.keys():
						if key not in ["obsolete_uni_id"]:
							self.logger[key][0].extend( logs_dict[key] )
							self.logger[key][1] += len( logs_dict[key] )
									
					# Do not consider PDBs for which no IDR chain was selected.
					pdb = list( entry_dict.keys() )[0]
					if len( entry_dict[pdb]["prot1"]["Asym ID"] ) == 0:
						self.logger["no_IDRs_in_PDB"][0].append( pdb )
						self.logger["no_IDRs_in_PDB"][1] += 1

					else:
						# with open( f"{self.idr_complexes_dir}{pdb}.json", "w" ) as w:
						# 	json.dump( entry_dict, w )
						self.idr_pdbs_selected[pdb] = None
						selected_uni_ids.update( entry_dict[pdb]["prot1"]["Uniprot ID"] )
						selected_uni_ids.update( entry_dict[pdb]["prot2"]["Uniprot ID"] )

				# if counter % 5000 == 0:
				print( "Selected complexes: ", len( self.idr_pdbs_selected ) ) 

		print( "Total Uniprot seq downloaded = ", len( self.uniprot_seq) )
		print( "Total selected Uniprot IDs = ", len( selected_uni_ids ) )
		with open( self.idr_pdbs_file, "w" ) as w:
			json.dump( self.idr_pdbs_selected, w )
		
		# Remove Uni ID that haven't been selected.
		self.uniprot_seq = {key: value for key, value in self.uniprot_seq.items() if key in selected_uni_ids}

		with open( self.sel_uniprot_seq_file, "w" ) as w:
			json.dump( self.uniprot_seq, w )



	def plot_expt_methods_count( self ):
		"""
		Create a pie plot for the experimental methods used for determining 
			structures of entries included in StrIDR.

		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""

		all_methods = {}
		for pdb in self.database_dict.keys():
			with open( f"{self.PDB_api_path}{pdb}_entry.json", "r" ) as f:
				entry_dict = json.load( f )
			method = entry_dict["exptl"][0]["method"]
			if method not in all_methods.keys():
				all_methods[method] = 0

			all_methods[method] += 1

		with open( "./All_methods.json", "w" ) as w:
			json.dump( all_methods, w )

		labels1, labels2 = [], []
		values1, values2 = [], []

		others = 0
		for key in all_methods.keys():
		    if all_methods[key] >3000:
		        labels1.append( f"{key} ({all_methods[key]})" )
		        values1.append( all_methods[key] )
		    else:
		        others += all_methods[key]
		        labels2.append( f"{key} ({all_methods[key]})" )
		        values2.append( all_methods[key] )

		labels1.append( f"Others ({others})" )
		values1.append( others )

		colors = ["green", "blue", "yellow", "red"]
		plt.pie( values1, labels = labels1, colors = colors )
		plt.savefig( "./Expt_sources.png", dpi = 300 )
		plt.close()

		colors = ["lightgreen", "maroon", "orange", "darkblue", "magenta", "cyan", "grey"]
		plt.pie( values2, labels = labels2, colors = colors )
		plt.savefig( "./Other_Expt_sources.png", dpi = 300 )
		plt.close()



	def plot_disorder_res_dist( self, counts, file_name ):
		fig, axis = plt.subplots( 1, 1, figsize = ( 10, 8 ) )
		axis.hist( counts, bins = 150 )
		axis.xaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis.yaxis.set_tick_params( labelsize = 14, length = 8, width = 2 )
		axis.set_ylabel( "Counts", fontsize = 16 )
		axis.set_xlabel( "No. of disordered residues", fontsize = 16 )
		label = " ".join( file_name.split( "_" ) )
		axis.set_title( f"{label}", fontsize = 16 )

		plt.savefig( f"./{file_name}.png", dpi = 300 )
		plt.close()



	def get_disorder_count_for_entries( self ):
		disorder_per_pdb = []
		disorder_per_chain = []

		for pdb in self.database_dict.keys():
			per_pdb = []
			for uni_id in self.database_dict[pdb].keys():
				num_disorder_res = self.database_dict[pdb][uni_id]["Count_disorder_residues"]

				if num_disorder_res[0] <= 1000:
					disorder_per_chain.extend( num_disorder_res )
				per_pdb.extend( num_disorder_res )
			if sum( per_pdb ) <= 1000:
				disorder_per_pdb.append( sum( per_pdb ) )
		
		## Uncomment below lines to print no of chains with >400 diosrdered residues.
		# x1 = np.array( disorder_per_chain )
		# x2 = np.array( disorder_per_pdb )
		# y1 = np.count_nonzero( np.where( x1 > 350, 1, 0 ) )
		# y2 = np.count_nonzero( np.where( x2 > 350, 1, 0 ) )
		# print( y1, "  ", y2 )
		# with open( "./Plot_logs.txt", "w" ) as w:
		# 	w.writelines( "Trimming the entries for creating plots.\nCutoff = 350\n" )
		# 	w.writelines( f"Entries with >350 disordered residues per chain: {y1}\n" )
		# 	w.writelines( f"Entries with >350 disordered residues per pdb: {y2}\n" )


		self.plot_disorder_res_dist( disorder_per_pdb, "Disordered_residues_per_PDB" )
		self.plot_disorder_res_dist( disorder_per_chain, "Disordered_residues_per_chain" )

		self.logger["counts"][">5_disorder_residue_pdbs"] = sum( [1 for count in disorder_per_pdb if count > 5] )
		self.logger["counts"][">5_disorder_residue_chains"] = sum( [1 for count in disorder_per_chain if count > 5] )



	def parallelize_database_creation( self ):
		"""
		Obtain all relevant info. for all selected entries in parallel.
		Log the following:
			Chains that were not mapped.
			Chains that were mapped to incorrect Uni ID.
			Chains with no disordered regions.
			PDBs with no IDRs.
		
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""
		print( "Initiating Database Creation module...\n" )
		# complexes_ = complexes

		if not os.path.exists( self.processed_entry_dir ):
			subprocess.call( ["mkdir", f"{self.processed_entry_dir}"] )
		
		counter = 0
		groups = self.pdb_info_df.groupby( "PDB ID" )
		groups = list( groups )

		count_disorder_res = []
		count_pdb_with_disorder = 0
		uni_ids_with_IDR_in_struct = set()

		for start in np.arange( 0, len(  groups ), self.batch_size ):
			end = start + self.batch_size
			batch = groups[start:end]

			with Pool( self.cores ) as p:
				for result in tqdm.tqdm( 
							p.imap_unordered( 
											self.the_forge, batch ), 
							total = len( batch ) ):
					# excess_uni_ids, not_mapped, incorrect_uni_map, entry_dict = result
					logs_dict, entry_dict = result
					# entry_dict = entry_dict[0]
					pdb = list( entry_dict.keys() )[0]

					counter += 1
					for key in logs_dict.keys():
						if key not in ["all_chains_post_download", "all_entities_post_download", 
										"uni_ids_with_IDR_in_struct", "pdbs_IDR_in_seq"]:
							
							self.logger[key][0].extend( logs_dict[key] )
							self.logger[key][1] += len( logs_dict[key] )
					self.logger["counts"]["pdbs_IDR_in_seq"] += logs_dict["pdbs_IDR_in_seq"]
					self.logger["counts"]["all_chains_post_download"] += logs_dict["all_chains_post_download"]
					self.logger["counts"]["all_entities_post_download"] += logs_dict["all_entities_post_download"]
					uni_ids_with_IDR_in_struct.update( logs_dict["uni_ids_with_IDR_in_struct"] )

					if entry_dict[pdb] == {}:
						self.logger["no_IDRs_in_PDB"][0].append( pdb )
						self.logger["no_IDRs_in_PDB"][1] += len( [pdb] )

					else:
						self.database_dict[pdb] = {}
						for uni_id in entry_dict[pdb].keys():
							self.database_dict[pdb][uni_id] = {}

							for key in entry_dict[pdb][uni_id].keys():
								self.database_dict[pdb][uni_id][key] = entry_dict[pdb][uni_id][key]

				print( "PDBs selected so far = ", len( self.database_dict.keys() ) )

		self.logger["counts"]["uni_ids_with_IDR_in_struct"] = len( uni_ids_with_IDR_in_struct )

		with open( f"{self.output_file_name}.json", "w" ) as w:
			json.dump( self.database_dict, w)

		print( "Total PDBs selected = ", len( self.database_dict.keys() ) )



	def save_logs( self ):
		"""
		Save the dataset/database.
		
		Input:
		----------
		Does not take any input arguments.

		Returns:
		----------
		None
		"""

		proc = subprocess.Popen( "hostname", shell = True, stdout = subprocess.PIPE )
		system = proc.communicate()[0]

		proc = subprocess.Popen( "date", shell = True, stdout = subprocess.PIPE )
		sys_date = proc.communicate()[0]

		self.tom = time.time()
		# Don't update total time for re-runs.
		if "total" not in self.logger["time_taken"].keys():
			self.logger["time_taken"]["total"] = self.tom-self.tim
		print( f"\nTotal time taken = {self.logger['time_taken']['total']/3600} hours" )

		with open( self.logger_file, "w" ) as w:
			json.dump( self.logger, w )

		with open( f"Logs_{self.output_file_name}.txt", "w" ) as w:
			w.writelines( "<---------------Logs Disobind--------------->\n" )
			w.writelines( f"Created on: System = {system} \t Date = {sys_date}\n" )
			w.writelines( f"CPU Cores used = {self.cores}\n" )
			w.writelines( f"Min protein length = {self.min_len}\n" )
			
			if self.create_dataset:
				w.writelines( f"Max no. of chains allowed in PDB = {self.max_num_chains}\n" )
				w.writelines( f"Max protein length = {self.max_len}\n" )
				w.writelines( f"Fragment if exceeding max length = {self.fragment_exceeding_maxlen} \n" )
				w.writelines( f"Min %disordered residues per IDR chain = {self.min_disorder_percent}\n" )
				w.writelines( "\tUnique Uniprot IDs = {}\n".format( len( self.uniprot_seq.keys() ) ) )
				# w.writelines( f"\tNon overlapping Uniprot ID pairs = {len( self.uniprot_pairs )}\n" )

			w.writelines( "\n\n-----------------------------------------------------------\n" )
			w.writelines( "\t\t Counts ------------\n" )
			w.writelines( "All PDBs collected = {}\n".format( self.logger["counts"]["all_PDBs"] ) )
			w.writelines( "Unique PDBs collected = {}\n".format( self.logger["counts"]["unique_PDBs"] ) )
			w.writelines( "Total PDBs downloaded = {}\n".format( self.logger["counts"]["PDB_download"] ) )
			w.writelines( "Total Uniprot seq downloaded = {}\n".format( len( self.uniprot_seq ) ) )
			w.writelines( "Total PDBs post removing unwanted Uni IDs = {}\n".format( self.logger["counts"]["PDBs_after_filtering_by_uniprot_ids"] ) )
			w.writelines( "PDBs with IDRs included = {}\n".format( self.logger["counts"]["IDR_PDBs_obtained"] ) )
			
			if self.create_dataset:
				w.writelines( "Total binary complexes = {}\n".format( self.logger["counts"]["total_binary_complexes"] ) )
				w.writelines( "Uniprot IDs post identifying IDRs = {}\n".format( self.logger["counts"]["Uniprot_IDs_remaining"] ) )
				w.writelines( "\tTotal Uniprot ID pairs = {}\n".format( self.logger["counts"]["uniprot_ID_pairs"] ) )
			else:
				# w.writelines( "PDBs with IDR in seq = {}\n".format( self.logger["counts"]["pdbs_IDR_in_seq"] ) )
				w.writelines( "PDBs with >5 disordered residues in struct = {}\n".format( 
															self.logger["counts"][">5_disorder_residue_pdbs"] ) )
				w.writelines( "Chains with >5 disordered residues in struct = {}\n".format( 
															self.logger["counts"][">5_disorder_residue_chains"] ) )
				w.writelines( "Entities with IDR in seq = {}\n".format( self.logger["entities_IDR_in_seq"][1] ) )
				w.writelines( "Entities with IDR in struct = {}\n".format( self.logger["entities_IDR_in_struct"][1] ) )
				w.writelines( "Chains with IDR in seq = {}\n".format( self.logger["chains_IDR_in_seq"][1] ) )
				w.writelines( "Chains with IDR in struct = {}\n".format( self.logger["chains_IDR_in_struct"][1] ) )
				w.writelines( "Uniprot IDs with IDR in struct = {}\n".format( self.logger["counts"]["uni_ids_with_IDR_in_struct"] ) )
				w.writelines( "Total chains in all downloaded PDBs = {}\n".format( self.logger["counts"]["all_chains_post_download"] ) )
				w.writelines( "Total entities in all downloaded PDBs = {}\n".format( self.logger["counts"]["all_entities_post_download"] ) )
				w.writelines( "Total chains in all PDBs before filtering = {}\n".format( self.logger["counts"]["all_chains_in_pdb"] ) )
				w.writelines( "Total unique Uniprot ID in all PDBs before filtering = {}\n".format( self.logger["counts"]["all_uni_ids_in_pdb"] ) )
			w.writelines( "\n-----------------------------------------------------------\n\n" )

			w.writelines( "\t\t Time taken ------------\n" )
			w.writelines( "Time taken for downloading PDBs = {} minutes\n".format( self.logger["time_taken"]["PDB_download"]/60 ) )
			w.writelines( "Time taken for downloading Uni seq = {} minutes\n".format( self.logger["time_taken"]["Uniprot_seq_download"]/60 ) )
			w.writelines( "Time taken for extracting relevant info from PDB = {} hours\n".format( self.logger["time_taken"]["IDR_PDBs_obtained"]/3600 ) )
			
			w.writelines( "Total = {} hours\n\n".format( self.logger["time_taken"]["total"]/3600 ) )
			w.writelines( "\n-----------------------------------------------------------\n\n" )
			
			for key in self.logger.keys():
				if key not in ["time_taken", "counts"]:
					w.writelines( "-->" + key + "\n" )
					w.writelines( "Count = " + str( self.logger[key][1] ) + "\n" )
					for i in self.logger[key][0]:
						w.writelines( str( i ) + "," )
					w.writelines( "\n========================================================\n" )


# <=======================================================> #
if __name__ == "__main__":
	parser = argparse.ArgumentParser( description = "Disobind mode: Download PDBs of IDRs in complexes and parse to obtain binary complexes of IDRs with a partner protein. \n StrIDR mode: Get the PDB structures, sequences, and associated info for the database of IDRs with structures. ")
	
	
	parser.add_argument( 
						"--create_dataset", "-d", dest = "d", 
						help = "If this flag is present, PDBs are parsed to create binary complexes of IDRs with a partner protein for the Disobind dataset. \
								If not, the PDBs are parsed to create the StrIDR database. ", 
						action="store_true", required = False, default = False )
	
	parser.add_argument( 
						"--max_cores", "-c", dest = "c", 
						help = "No. of cores to be used.", 
						type = int, required = False, default = 10 )
	
	cores = parser.parse_args().c

	create_dataset = parser.parse_args().d

	parse_pdbs_for_idrs( cores, create_dataset ).forward()

	print( "May the Force be with you..." )

# <=======================================================> #

