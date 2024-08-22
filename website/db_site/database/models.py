from django.db import models
from django.urls import reverse

class ChainID( models.Model ):
	asym_id = models.TextField()
	auth_asym_id = models.TextField()


class DisorderResidues( models.Model ):
	disorder_in_seq = models.TextField()
	disorder_in_struct = models.TextField()


class CrossRefs( models.Model ):
	disprot = models.TextField()
	ideal = models.TextField()
	mobidb = models.TextField()


class Uniprot( models.Model ):
	uni_id = models.CharField( max_length = 20 )
	name = models.CharField( max_length = 500 )


class PDB( models.Model ):
	pdb_id = models.CharField( max_length = 6, primary_key = True )
	uniprot = models.ManyToManyField( Uniprot, through = "UniprotPDB" )

	def get_absolute_url( self ):
		return reverse( "database:pdb_detail", kwargs = {"pdb_id": self.pdb_id} )


class UniprotPDB( models.Model ):
	pdb = models.ForeignKey( PDB, on_delete = models.CASCADE )
	uniprot = models.ForeignKey( Uniprot, on_delete = models.CASCADE )
	chain_ids = models.ManyToManyField( ChainID )
	disorder_residues = models.ManyToManyField( DisorderResidues )
	cross_refs = models.ManyToManyField( CrossRefs )

