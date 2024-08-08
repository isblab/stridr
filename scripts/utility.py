###################### Contains helper functions #################
######### ------>"May the Force serve u well..." <------##########
##################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import pandas as pd
import os
import requests
import subprocess

from dataset.from_APIs_with_love import *


#########################################################
def convert_to_str( plist, add = "null", ):
	"""
	Convert an input list elements to char.

	Input:
	----------
	plist --> input list to be type casted to str.
	add --> convert "null" or nan to "null" by default.

	Returns:
	----------
	new_list --> input list type casted to str.
	"""
	new_list = []

	for i in range( len( plist ) ):
		if pd.isna( plist[i] ):
			new_list.append( add )
		elif  plist[i] == "null":
			new_list.append( add )
		else:
			try:
				new_list.append( str( int( float( plist[i] ) ) ) )
			except:
				nums = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
				tmp = ""
				for ch in plist[i]:
					if ch in nums:
						tmp += ch
				new_list.append( tmp )
				
	return new_list

# ## Test Cases
# x1 = [1, 2 , 3, 4, 5]
# ans1 = ["1", "2", "3", "4", "5"]
# x2 = [1.0, 2.0, 3.0, 4, 5]
# ans2 = ["1", "2", "3", "4", "5"]
# x3 = ["null", 2.0, 3.0, "null", 5]
# ans3 = ["0", "2", "3", "0", "5"]
# for xx, a in zip( [x1, x2, x3], [ans1, ans2, ans3] ):
# 	y = convert_to_str( xx )
# 	print(y)
	# if all([a[i] == y[i] for i in range( len( a ) )] ):
	# 	print( "Test Passed..." )
	# else:
	# 	raise Exception( "Test Failed..." )
# exit()


#########################################################
def get_intersection( pos_x, pos_y ):
	"""
	Obtain the intersecting positions in the input lists.

	Input:
	----------
	pos_x --> PDB/Uniprot seq pos.
	pos_y --> Fuzzy region pos.

	Returns:
	----------
	intersect --> sorted list of intersecting elements.

	e.g. pos_x = [1,2,3,4]; pos_y = [3,4,5,6]
			intersection --> [3,4]
	"""
	if len(pos_x) > len( pos_y ):
		intersect = set( pos_x ).intersection( pos_y )
	
	else:
		intersect = set( pos_y ).intersection( pos_x )

	return sorted( intersect )


#########################################################
def get_overlap( seq1, seq2 ):
	"""
	Obtain the overlapping residue positions in query and target.

	Input:
	----------
	seq1 --> uniprot position for query.
			could be str, float, int.
	seq2 --> uniprot position for target.
			must be a list of lists of type int.

	Returns:
	----------
	overlap --> list of overlapping residues in query and target.
	"""
	overlap = []
	seq1 = convert_to_str( seq1 )
	seq1 = list( map( int, map( float, seq1 ) ) )
	[overlap.append( get_intersection( seq1, x ) ) for x in seq2]
	return overlap


#########################################################
def check_for_overlap( uni_pos1, uni_pos2, ignore_boundary = True ):
	"""
	A cheaper way to check for overlap - won't return overlapping residues.
	Assumes that the list of residue positions is continous i.e. has no missing residues.
	
	Input:
	----------
	start1 --> first residue for entry1 (int).
	start2 --> first residue for entry2 (int).
	end1 --> last residue for entry1 (int).
	end2 --> last residue for entry2 (int).
	ignore_boundary --> if True, ignore overlap if only boundary residue is present.
		e.g. [1, 2, 3, 4] and [4, 5, 6, 7] --> no overlap if ignore_boundary is True.

	Returns:
	----------
	overlap --> bool
			True if there is overlap in Uniprot positions else False.

	e.g. entry1 --> [1, 2, 3, 4, 5, 6]; entry2 --> [5, 6, 7, 8, 9, 10]
		start1 = 1; end1 = 6; start2 = 5; end2 = 10
		= max( 0, ( min( end1, end2 ) - max( start1, start2 ) + 1 ) )
		= max( 0, ( min( 6, 10 ) - max( 1, 5 ) + 1 ) )
		= max( 0, ( 6 - 5 + 1 ) )
		= max( 0, 2 ) = 2
	+1 at the end ensures overlap for only 1 terminal overlapping residue.
	"""
	boundary = 0 if ignore_boundary else 1
	start1, end1 = uni_pos1[0], uni_pos1[-1]
	start2, end2 = uni_pos2[0], uni_pos2[-1]
	overlap = max( 0, ( min( end1, end2 ) - max( start1, start2 ) + boundary ) )
	# Regions are overlaping if overlap is positive.
	overlap = overlap > 0
	
	return overlap


#########################################################
def merge_residue_positions( res_pos1, res_pos2 ):
	"""
	Merge the given input lists/np.array containing protein residue positions.
	Assumes that the input has continous residue positions.

	Input:
	----------
	res_pos1 --> list/np.array of residue positions.
	res_pos2 --> list/np.array of residue positions.

	Returns:
	----------
	merged_pos --> np.array containing the merged sorted list.
	"""
	merged_pos = np.sort( np.unique( np.concatenate( ( res_pos1, res_pos2 ) ) ) )

	return merged_pos



#########################################################
def merged_seq_exceeds_maxlen( res_pos1, res_pos2, max_len ):
	"""
	Check if the input lists/np.array would exceed max_len upon merging.
	
	Input:
	----------
	res_pos1 --> list/np.array of residue positions.
	res_pos2 --> list/np.array of residue positions.
	max_len --> maximum length of the merged sequence allowed.

	Returns:
	----------
	True if the merged array exceeds max_len else False.
	"""
	merged_pos = merge_residue_positions( res_pos1, res_pos2 )

	if len( merged_pos ) > max_len:
		return True
	else:
		return False



#########################################################
def merge_overlapping_tuples( disorder_regions ):
	"""
	Merge overalpping regions.

	Input:
	----------
	disorder_regions --> list of tuples containing start end elements of a list.

	Returns:
	----------
	merged_regions --> list of tuples with overlapping tuples merged.
	"""
	merged_regions = []
	disorder_regions = sorted( disorder_regions )

	for i, x in enumerate( disorder_regions ):
		x1, x2 = x
		for y in disorder_regions[i+1:]:
			y1, y2 = y

			# x --> (1,7) and y --> (8,17) == (1,17)
			if abs( y2 - x1 ) == 1 or abs( x2 - y1 ) == 1:
				x1, x2 = min( x[0], y[0] ), max( x[1], y[1] )
				disorder_regions.remove( y )

			# x --> (1,7) and y --> (10,17) or x --> (20,30) and y --> (10,17)
			# Forget about the other non-overlapping ones.
			elif x1 > y2 or x2 < y1:
				continue
			
			# x --> (1,12) and y --> (10,17) == (1,17)
			elif x1 <= y1 and x2 <= y2:
				x1, x2 = x1, y2
				disorder_regions.remove( y )

			# x --> (10,20) and y --> (12,19) == (10,20)
			elif x1 <= y1 and x2 >= y2:
				x1, x2 = x1, x2
				disorder_regions.remove( y )
			
			# x --> (12,20) and y --> (10,17) == (10,20)
			elif x1 >= y1 and x2 >= y2:
				x1, x2 = y1, x2
				disorder_regions.remove( y )
			
			# x --> (11,17) and y --> (10,19) == (10,19)
			elif x1 >= y1 and x2 <= y2:
				x1, x2 = y1, y2
				disorder_regions.remove( y )
		merged_regions.append( ( x1, x2 ) )

	return merged_regions

## Test Cases
# # x1 --> Output: [(1,2), (5,10), (15,25)]
# x1 = [(1,2), (5, 10), (15, 25)] 
# # x2 --> Output: [(1,2), (5, 25)]
# x2 = [(1,2), (5, 10), (8, 25)]
# # x3 --> Output: [(1,2), (5, 34)]
# x3 = [(1,2), (5, 10), (8, 25), (26, 34)]
# # x4 --> Output: [(1,2), (5, 37)]
# x4 = [(1,2), (5, 10), (8, 25), (26, 34), (10, 20), (32, 37)]
# # x5 --> Output: [(1,37)]
# x5 = [(1,2), (5, 10), (8, 25), (26, 34), (1, 20), (32, 37)]

# print( merge_overlap( x5 ) )
# exit()


#########################################################
def consolidate_regions( disorder_regions, min_len ):
	"""
	Given a comma separated string of disorder regions, 
		merge adjacent overlapping regions.
	Include only regions larger than min_len.

	Input:
	----------
	disorder_regions --> comma separated string.
		e.g. "start1-end1,start2-end2"
	min_len --> (int) min length of the disordered region.

	Returns:
	----------
	merged_disorder_regions --> list of tuples containing 
		start and end residue positions for a disordered region. 
	"""

	merged_disorder_regions = []
	disorder_regions = [( int( x.split( "-" )[0] ), int( x.split( "-" )[1] ) ) for x in disorder_regions.split( "," )]
	disorder_regions = merge_overlapping_tuples( disorder_regions )
	merged_disorder_regions = [x for x in disorder_regions if ( x[1] - x[0] + 1 ) >= min_len]

	return merged_disorder_regions


#########################################################
def load_disorder_dbs( disprot_path = None, ideal_path = None, mobidb_path = None ):
	"""
	Load the csv files corresponding to the disordered protein databases
		DisProt, IDEAL, MobiDB.

	Input:
	----------
	disprot_path --> path for the DisProt csv file.
	ideal_path --> path for the IDEAL csv file.
	mobidb_path --> path for the MobiDB csv file.

	Returns:
	----------
	disprot --> csv file for the DisProt database.
	ideal --> csv file for the IDEAL database.
	mobidb --> csv file for the mobiDB database.
	"""
	disprot = pd.read_csv( disprot_path )
	ideal = pd.read_csv( ideal_path )
	mobidb = pd.read_csv( mobidb_path )

	return disprot, ideal, mobidb


def find_disorder_regions( disprot, ideal, mobidb, uni_ids, min_len = 1, return_ids = False ):
	"""
	Obtain the overlapping disordered region for the given residue positions.
	Overlap the uni positions with all disordered regions on disprot or IDEAL.

	Input:
	----------
	disprot --> csv file for the DisProt database.
	ideal --> csv file for the IDEAL database.
	mobidb --> csv file for the mobiDB database.
	uni_ids --> list of Uniprot IDs.
	min_len --> min length of disorder region to be considered.
	return_ids --> return the respective database IDs (Default = False).

	Returns:
	----------
	disordered_regions --> list of lists containing disordered regions positions.
	db_ids --> database identifier for the entry.
	"""
	overlap = []
	disorder_regions = []

	disprot_regions, disprot_ids = [], []
	ideal_regions, ideal_ids = [], []
	mobidb_regions, mobidb_ids = [], []
	
	for id_ in uni_ids:
		disprot_subset = disprot[disprot["Uniprot ID"].str.contains( id_ )]
		disprot_regions.extend( disprot_subset["Disorder regions"].tolist() )
		if return_ids:
			disprot_ids.extend( disprot_subset["Disprot ID"].tolist() )

		ideal_subset = ideal[ideal["Uniprot ID"].str.contains( id_ )]
		ideal_regions.extend( ideal_subset["Disorder regions"].tolist() )
		if return_ids:
			ideal_ids.extend( ideal_subset["IDP ID"].tolist() )

		mobidb_regions.extend( mobidb[mobidb["Uniprot ID"].str.contains( id_ )]["Disorder regions"].tolist() )
		if return_ids:
			mobidb_ids += uni_ids if len( mobidb_regions ) != 0 else []


	disorder_regions = ",".join( disprot_regions + ideal_regions + mobidb_regions )
	if return_ids:
		db_ids = f"{','.join( disprot_ids )}--{','.join( ideal_ids )}--{','.join( mobidb_ids )}"

	if len( disorder_regions ) != 0:
		disorder_regions = consolidate_regions( disorder_regions, min_len )
		disordered_pos = [np.arange( x[0], x[1] + 1 ) for x in disorder_regions]
	else:
		disordered_pos = []

	if return_ids:
		return disordered_pos, db_ids
	else:
		return disordered_pos



#########################################################
def change_basis( mapped_uni_pos, mapped_pdb_pos, target_pos, add = "null", forward = True ):
	"""
	Convert PDB positions to Uniprot positions and vice-versa.
	For a protein, residue positions in Uniprot and those in PDB often differ.
	Given mapped PDB and Uniprot positions, one 
			can just convert the PDB positions to Uniprot positions.

	Input:
	----------
	mapped_uni_pos --> Uniprot positions mapped to PDB using SIFTS (or PDBSWS).
	mapped_pdb_pos --> PDB positions mapped to Uniprot using SIFTS (or PDBSWS).
	target_pos --> must be a list of PDB positions if forward is true else a list of uni positions.
	add --> by default convert "null"/ nan in PDB input/ target to "null".
	forward --> If true converts the target PDB pos to Uniprot pos else converts the target Uniprot pos to PDB pos.

	Returns:
	----------
	modified_target --> target_pos in Uni/PDB basis as specified.
	"""
	modified_target, considered = [], []
	# Convert all elements to char.
	mapped_uni_pos = convert_to_str( mapped_uni_pos, add = "null" )
	mapped_pdb_pos = convert_to_str( mapped_pdb_pos, add = "null" )
	target_pos = convert_to_str( target_pos, add = add )

	for i in range( len( target_pos ) ):
		if target_pos[i] == "null":
			modified_target.append( "null" )
		
		# Uniprot --> PDB.
		elif not forward:
			idx = mapped_uni_pos.index( target_pos[i] )
			modified_target.append( mapped_pdb_pos[idx] )
		
		# PDB --> Uniprot
		elif forward:
			idx = mapped_pdb_pos.index( target_pos[i] )
			modified_target.append( mapped_uni_pos[idx] )

	return modified_target


# # ## Test Cases
# pdb = ["null", "null", "null", 11, 12, 13, "null", "null", 16, 17, 18, 19, 20]
# uni = [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
# indices = np.arange( 0, len( pdb ), 1 )

# ans1 = ["null", "null", "null", "38", "39", "40", "null", "null", "43", "44", "45", "46", "47"]
# ans2 = ["null", "11", "12", "13", "null", "null", "16"]

# new_basis = change_basis( uni, pdb, pdb, add = "null", forward = True )
# print( new_basis )
# if all( [ans1[i] == new_basis[i] for i in range( len( ans1 ) )] ):
# 	print( "Test1 Passed...\n" )
# else:
# 	raise Exception( "Test1 Failed...\n" )

# old_basis = change_basis( uni, pdb, new_basis[2:9], add = "null", forward = False )
# print( old_basis )
# if all( [ans2[i] == old_basis[i] for i in range( len( ans2 ) )] ):
# 	print( "Test2 Passed..." )
# else:
# 	raise Exception( "Test2 Failed..." )

# exit()


def change_basis2( mapping, target_pos, indices, add = "null" ):
	"""
	Convert PDB positions to Uniprot positions and vice-versa.
		This version maps according to an index list.
		Needs "null" to be removed priorly.

	Input:
	----------
	mapping --> Uniprot/PDB positions obtained using SIFTS (or PDBSWS).
				Should be, mapped_uni_pos:: PDB --> Uniprot conversion.
				Should be, mapped_pdb_pos:: Uniprot --> PDB conversion.
	target_pos --> list of residue position that need to be converted to Uniprot or PDB basis.
	indices --> an index list for the mapped positions.
	add --> by default convert "null"/ nan in PDB input/ target to "null".

	Returns:
	----------
	modified_target --> target_pos in Uni/PDB basis as specified.
	"""
	modified_target = []
	# Convert all elements to char.
	mapping = convert_to_str( mapping, add = "null" )
	target_pos = convert_to_str( target_pos, add = add )

	if "null" in target_pos:
		raise Exception( "Invalid position in target_pos...\n" )

	for i, idx in enumerate( indices ):
		modified_target.append( mapping[idx] )
		
	return modified_target

# # ## Test Cases
# pdb = ["8", "9", "10", 11, 12, 13, "14", "15", 16, 17, 18, 19, 20]
# uni = [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
# indices = np.arange( 0, len( pdb ), 1 )

# ans1 = ["35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47"]
# ans2 = ["10", "11", "12", "13", "14", "15", "16"]

# new_basis = change_basis2( uni, pdb, indices, add = "null" )
# print( new_basis )
# if all( [ans1[i] == new_basis[i] for i in range( len( ans1 ) )] ):
# 	print( "Test1 Passed...\n" )
# else:
# 	raise Exception( "Test1 Failed...\n" )

# old_basis = change_basis2( pdb, new_basis[2:9], indices[2:9], add = "null" )
# print( old_basis )
# if all( [ans2[i] == old_basis[i] for i in range( len( ans2 ) )] ):
# 	print( "Test2 Passed..." )
# else:
# 	raise Exception( "Test2 Failed..." )

# exit()


#########################################################
def add_residue_positions( query_pos, target_pos, add_null = True ):
	"""
	Add missing residues as "nulls" in target using query_pos as reference.
	e.g.
	query_pos = [11, 12, 13, 14, 15, 16, 17, 18]
	target_pos = [13, 14, 15, 16, 17] 
	Output target_pos = ["null", "null", 13, 14, 15, 16, 17, "null"]

	Input:
	----------
	query_pos --> reference list of positions.
	target_pos --> list containing missing residues.
	add_null --> if True, adds a "null" else adds the missing residue position.	

	Returns:
	----------
	positions --> target_pos with "nulls" added for missing residues as per query_pos.
	"""
	positions = []

	for i, pos in enumerate( query_pos ):
		add = "null" if add_null else pos
		if pos not in target_pos:
			positions.append( add )
		else:
			positions.append( pos )
	if len( positions ) != len( query_pos ):
		raise Exception( "Excessive residues added to the target..." )
	return positions

# ## Test Cases
# uni = np.arange( 11, 20, 1 )
# pdb1 = [14, 15, 16, 17, 18]
# ans1 = ["null", "null", "null", 14, 15, 16, 17, 18, "null"]
# pdb2 = [11, 12, 13, 14]
# ans2 = [11, 12, 13, 14, "null", "null", "null", "null", "null"]
# pdb3 = [18, 19]
# ans3 = ["null", "null", "null", "null", "null", "null", "null", 18, 19]
# for x, a in zip([pdb1, pdb2, pdb3], [ans1, ans2, ans3]):
# 	y1 = add_residue_positions( uni, x, True )
# 	if all( [a[i] == y1[i] for i in range( len( a ) )] ):
# 		print( "Test Passed..." )
# 	else:
# 		raise Exception( "Test Failed..." )
# 	# y2 = add_residue_positions( uni, x, False )
# exit()



#########################################################
def remove_nulls( positions, index_list = [] ):
	"""
	Remove all missing residues appearing as "null"/nan from the input.
	e.g. [1, 2, 3, 4, "null", "null", "null", 8, 9, 10]
		[1, 2, 3, 4, 8, 9, 10]

	Input:
	----------
	positions --> list containing residue positions.
	index_list --> list of indices for the positions list.

	Returns:
	----------
	new_pos --> input residue positions without "nulls".
	new_ind --> index list with missing residue indices removed.
	"""
	new_pos, new_ind = [], []
	positions = convert_to_str( positions, add = "null" )
	# positions = [int( pos ) for pos in positions if pos != "null"]
	for idx, pos in enumerate( positions ):
		if pos != "null":
			new_pos.append( pos )
			if len( index_list ) != 0:
				new_ind.append( idx )

	return new_pos, new_ind

# ## Test Cases
# x1 = ["null", "null", 11, 12, 13, 14, 15, 16.0, 17, 18.0, 20.0, "null", "null", "null", "null", "25.0", "26.0", "27.0", "28.0"]
# ans1 = [11, 12, 13, 14, 15, 16, 17, 18, 20, 25, 26, 27, 28]
# output = remove_nulls( x1 )
# print( output )
# if all( [y1 == y2 for y1, y2 in zip( output, ans1 )] ):
# 	print( "Test Passed..." )
# else:
# 	raise Exception( "Test Failed..." )
# exit()


def remove_nulls2( positions, index_list = [] ):
	"""
	Remove all missing residues appearing as "null"/nan from the input.
		This function returns a list of lists for each fragment obtained 
			post removing missing residues.
	e.g. [1, 2, 3, 4, "null", "null", "null", 8, 9, 10]
		[[1, 2, 3, 4], [8, 9, 10]]

	Input:
	----------
	positions --> list containing residue positions.
	index_list --> list of indices for the positions list.

	Returns:
	----------
	new_pos --> input residue positions without "nulls".
	new_ind --> index list with missing residue indices removed.
		Returned only if an input index list is provided.
	"""
	new_pos, new_ind, tmp_pos, tmp_ind = [], [], [], []
	positions = convert_to_str( positions, add = "null" )
	
	#TODO not really using the values in index_list anywhere 
	for idx, pos in enumerate( positions ):
		if pos != "null":
			tmp_pos.append( pos )
			if len( index_list ) != 0:
				tmp_ind.append( idx )
		else:
			if tmp_pos != []:
				new_pos.append( tmp_pos )
				if len( index_list ) != 0:
					new_ind.append( tmp_ind )
			tmp_pos, tmp_ind = [], []

	if tmp_pos != []:
		new_pos.append( tmp_pos )
		if len( index_list ) != 0:
			new_ind.append( tmp_ind )	

	return new_pos, new_ind

# ## Test Cases
# x1 = ["null", "null", 11, 12, 13, 14, 15, 16.0, 17, 18.0, 20.0, "null", "null", "null", "null", "25.0", "26.0", "27.0", "28.0"]
# ans1 = [["11", "12", "13", "14", "15", "16", "17", "18", "20"], ["25", "26", "27", "28"]]
# output, _ = remove_nulls2( x1 )
# print( output )
# if all( [y1 == y2 for y1, y2 in zip( output, ans1 )] ):
# 	print( "Test Passed..." )
# else:
# 	raise Exception( "Test Failed..." )
# exit()


#########################################################
def ranges( positions ):
	"""
	Get tuples of continous residue positions.
	e.g. [1, 2, 3, 4, 7, 8, 9, 10]
		[(1, 4), (7, 10)]

	Input:
	----------
	positions --> list containng residue positions.

	Returns:
	----------
	Returns list of tuples of start and end residues positions 
		for a continous set of residues.
	"""
	positions = list( map( int, positions ) )
	positions = sorted( set( positions ) )
	# Get start, end positions for each continous confident patch.
	gaps = [[x, y] for x, y in zip( positions, positions[1:] ) if x+1 < y]
	edges = iter( positions[:1] + sum( gaps, [] ) + positions[-1:] )
	return list( zip( edges, edges ) )


# ## Test Cases
# x1 = [11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 40, 41, 42, 43, 67, 68, 69, 70]
# ans1 = [(11, 18), (22, 25), (40, 43), (67, 70)]
# output = ranges( x1 )
# print( output )
# if all( [y1 == y2 for y1, y2 in zip( output, ans1 )] ):
# 	print( "Test Passed..." )
# else:
# 	raise Exception( "Test Failed..." )
# exit()



#########################################################
#########################################################
# MMSeqs2
def mmseqs_cluster( db_file, out_file, algo, min_seq_id, cluster_mode ):
	"""
	Use MMSeqs2 for performing seq based clustering.

	Input:
	----------
	db_file --> FASTA file containing all sequences to be clustered.
	out_file --> name for the output files generated.
	algo --> either of easy-cluster or easy-linclust.
	min_seq_id --> minimum seq identity threshold used for clustering.
			Sensitivity of clustering isadjusted accordingly.
	cluster_mode --> Algorithm used for clustering.
			0 --> greedy set (gives the least clusterd).
			1 --> connected component (gives most clusters).
			2 --> greedy incremental (same as CD-HIT).
	
	Returns:
	----------
	None
	"""
	with open( "./mmseqs2_redir.txt", 'w' ) as w:
		subprocess.call( 
						["mmseqs", f"{algo}",
						db_file, out_file, "tmp",
						"--min-seq-id", f"{min_seq_id}",
						"--cluster-mode", f"{cluster_mode}",
						], stdout = w )


def read_mmseqs_tsv_output( tsv_file ):
	"""
	Read the clustering output tsv file generated by MMSeqs2.

	Input:
	----------
	tsv_file --> PATH for the tsv file.

	Returns:
	----------
	df --> pandas dataframe for the input tsv file.
	"""
	df = pd.read_csv( tsv_file, index_col = False, header = None, sep = "\t" )

	return df


# MMalign
def mmalign( pdb1, pdb2, redir_file = "./mmalign_tmp.txt" ):
	"""
	Use MMalign to perform structural alignment fo pdb1 and pdb2.

	Input:
	----------
	pdb1 --> PATH for PDB1 file to be aligned (will consider this as reference).
	pdb2 --> PATH for PDB2 file to be aligned.
	redir_file --> file to store the MMalign output.


	Returns:
	----------
	None
	"""
	import os
	abs_path = os.getcwd()
	print( abs_path )
	mmalign_script = os.path.join( abs_path.split( "Database" )[0], "dataset/MMalign" )
	with open( redir_file, 'w' ) as w:
		subprocess.call( 
						[f"{mmalign_script}", f"{pdb1}", f"{pdb2}"
						], stdout = w )


# USalign
def usalign( usalign_script, pdb1, pdb2, chain1, chain2, mol, mm, ter, redir_file = "./usalign_tmp.txt" ):
	"""
	Use USalign to perform structural alignment fo pdb1 and pdb2.

	Input:
	----------
	script_path --> absolute path for the USalign executable.
	pdb1 --> PATH for PDB1 file to be aligned (will consider this as reference).
	pdb2 --> PATH for PDB2 file to be aligned.
	chain1 --> comma separated Auth asym IDs for PDB1 to be aligned.
	chain2 --> comma separated Auth asym IDs for PDB2 to be aligned.
	mol --> molecule type [auto, prot, RNA.
	mm --> multimeric laignment option.
		0: (default) alignment of two monomeric structures.
		1: alignment of two multi-chain oligomeric structures.
		2: alignment of individual chains to an oligomeric structure.
		Look at USalign -h option for more details.
	ter --> #chains to align.
		0: align all chains from all models.
		1: align all chains of the first model.
		2: (default) only align the first chain.
		Look at USalign -h option for more details.
	redir_file --> file to store the MMalign output.

	Returns:
	----------
	None
	"""
	# USalign -chain1 C,D,E,F 5jdo.pdb -chain2 A,B,C,D 3wtg.pdb -ter 0
	if len( chain1.split( "," ) ) == 1 and len( chain2.split( "," ) ) == 1:
		mm = 0
	
	elif len( chain1.split( "," ) ) > 1 and len( chain2.split( "," ) ) > 1:
		mm = 1
	
	else:
		mm = 2

	with open( redir_file, 'w' ) as w:
		subprocess.call( 
						[f"{usalign_script}", 
						"-chain1", f"{chain1}", f"{pdb1}", 
						"-chain2", f"{chain2}", f"{pdb2}",
						"-mol", f"{mol}",
						"-mm", f"{mm}",
						"-ter", f"{ter}"], 
						stdout = w,
						stderr = subprocess.STDOUT )


def get_aligned_TM_score( redir_file = "./align_tmp.txt" ):
	"""
	Read the MMalign/USalign output stored in a txt file.
		Line15: TM-score= 0.XXXXX (if normalized by length of Chain_1, i.e., LN=XX, d0=X.XX)
		Line16: TM-score= 0.XXXXX (if normalized by length of Chain_2, i.e., LN=XXX, d0=X.XX)

	Input:
	----------
	redir_file --> file containing the MMalign/UMalign output.

	Returns:
	----------
	tm1 --> TM-score normalized by length of Chain_1.
	tm1 --> TM-score normalized by length of Chain_2.
	"""
	with open( redir_file, "r" ) as f:
		output = f.readlines()
	tm1, tm2 = None, None
	for i in range( len( output ) ):
		if output[i].startswith( "TM-score=" ):
			line15 = output[i].split( " " )
			line16 = output[i+1].split( " " )
			tm1 = float( line15[1] )
			tm2 = float( line16[1] )
			break

	return tm1, tm2


