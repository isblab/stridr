from django import forms

class SearchForm( forms.Form ):
	search_input = forms.CharField( 
					label = "Search", 
					max_length = 15,
					widget=forms.TextInput( attrs={"placeholder": "Enter PDB ID/Uniprot ID/Disprot ID/IDEAL ID/MobiDB ID",
													"id": "searchBar",
													"style": "width: 80%; margin: 0 auto; text-align: center; font-family: verdana"} ) )

	# def search( self ):
	# 	search_input = self.cleaned_data.get( "search_input" )
	# 	print( search_input, "------------------------" )

	# 	results = PDB.objects.filter(
	# 		Q( pdb_id__icontains = search_input) | Q( uniprot__uni_id__icontains = search_input )
	# 		).distinct()

	# 	return results
