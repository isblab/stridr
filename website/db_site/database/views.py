from django.shortcuts import render
from django.http import FileResponse, StreamingHttpResponse, HttpResponse
from wsgiref.util import FileWrapper
from django.views import View
from django.views.generic import ListView, DetailView, TemplateView
from django.core.paginator import Paginator
from database.models import PDB, Uniprot, UniprotPDB
from database.forms import SearchForm
from django.db.models import Q

import os, gzip, shutil

# Create your views here.
# def home( request ):
# 	return render(request, "home.html" )

def home( request ):
	form = SearchForm()
	return render(request, "home.html", {"form": form})


class PDBListView( ListView ):
	model = PDB
	template_name = "browse.html"
	context_object_name = "entry_list"
	paginate_by = 20

	def get_context_data( self, **kwargs ):
	    context = super().get_context_data( **kwargs )
	    context["form"] = SearchForm( self.request.GET )
	    return context

	def post(self, request, *args, **kwargs):
		form = SearchForm(request.POST)
		if form.is_valid():
			search_input = form.cleaned_data.get( "search_input" )
			return redirect( "database:search_view", search_input = search_input )



class PDBDetailView( DetailView ):
	model = PDB
	# model = UniprotPDB
	template_name = "pdb_detail.html"
	context_object_name = "pdb"

	def get_object( self, queryset = None ):
		pdb_id = self.kwargs["pdb_id"]
		return PDB.objects.get( pdb_id = pdb_id )



class SearchListView( ListView ):
	model = PDB
	context_object_name = "results"
	template_name = "search.html"
	paginate_by = 20
	
	def get( self, request ):
		form = SearchForm( request.GET )
		results = None

		if form.is_valid():
			search_input = form.cleaned_data.get( "search_input" )

			results = PDB.objects.filter(
				Q( pdb_id__icontains = search_input) | 
				Q( uniprot__uni_id__icontains = search_input ) |
				Q( uniprot__uniprotpdb__cross_refs__disprot__icontains = search_input ) |
				Q( uniprot__uniprotpdb__cross_refs__ideal__icontains = search_input ) |
				Q( uniprot__uniprotpdb__cross_refs__mobidb__icontains = search_input )
				).distinct()

		return render( request, "search.html", {"results": results} )


def download( request ):
	return render(request, "download.html")


class DownloadView( TemplateView ):
	template_name = "download.html"
	
	def get( self, request, filename ):
		base_path = "./database/static"
		file_path = os.path.join( base_path, filename )

		if os.path.exists( file_path ):
			with open( file_path, "rb" ) as f:
				file = f.read()
				if filename not in ["StrIDR_database.json", "Uniprot_seq.json"]:
					file = gzip.compress( file )
					filename = f"{os.path.basename( filename )}.gz"

				response = HttpResponse( file, content_type = "'application/octet-stream'" )
				response["Content-Disposition"] = f'attachment; filename={filename}'

				return response
		else:
			return HttpResponse( "File Not Found" )


