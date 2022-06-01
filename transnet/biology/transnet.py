"""Connecting transomics network
"""

__all__ = [
    "Transnet",
]

from transnet.api.kegg import *
from transnet.api.uniprot import *
from transnet.api.brenda import *
from transnet.api.string import *
from transnet.biology.elements import *
from transnet.biology.layers import *

class Transnet():
    """
    Transnet
    """

    def __init__(
        self,
        ncbi_organism: str = None,
        kegg_organism: str = None,
        organism_full: str = None,
        pathways: Pathways = None,
        transcriptome: Transcriptome = None,
        proteome: Proteome = None,
#        fluxome: Fluxome = None,
#        metabolome: Metabolome = None,
#        genome: Genome = None,
        ):
        self.ncbi_organism = ncbi_organism
        self.kegg_organism = kegg_organism
        self.organism_full = organism_full
        self.pathways = pathways
        self.transcriptome = transcriptome
        self.proteome = proteome
#        self.fluxome = fluxome
#        self.metabolome = metabolome
#        self.genome = genome

        self.initialize_apis()

    def __repr__(self):
        return f"<Transnet for {self.organism_full}>"

    def initialize_apis(self):
        """
        Initializes the APIs
        """
        self.brenda = BRENDA_api()
        print('APIs initialized')
    
    def generate(self):
        """
        Generates the transnet for the organism
        """
        print('- Generating transnet')
        print('---------------------')

        self.generate_pathways()
        
        self.generate_transcriptome()
        
        self.generate_proteome()
        
        print('- DONE!')
        print('-------------------')
        #self.fluxome = self.generate_fluxome()
        #self.metabolome = self.generate_metabolome()
        #self.genome = self.generate_genome()

    def generate_pathways(self):
        """
        Generates the pathways for the organism
        """
        print('- Generating pathways')
        self.pathways = Pathways(kegg_organism=self.kegg_organism)
        self.pathways.populate()
        self.pathways.fill_pathways()
        print('- Pathways generated')

    def generate_transcriptome(self):
        """
        Generates the transcriptome for the organism
        """
        print('- Generating transcriptome')
        self.transcriptome = Transcriptome(kegg_organism=self.kegg_organism)
        self.transcriptome.populate()
        self.transcriptome.fill_genes()
        print('- Transcriptome generated')

    def generate_proteome(self):
        """
        Generates the proteome for the organism
        """
        print('- Generating proteome')
        self.proteome = Proteome(ncbi_organism=self.ncbi_organism, brenda_api=self.brenda)
        self.proteome.populate()
        self.proteome.fill_proteins()
        print('- Proteome generated')

#    def generate_fluxome(self):
#        """
#        Generates the fluxome for the organism
#        """
#        self.fluxome = Fluxome()
#
#    def generate_metabolome(self):
#        """
#        Generates the metabolome for the organism
#        """
#        self.metabolome = Metabolome()
#
#    def generate_genome(self):
#        """
#        Generates the genome for the organism
#        """
#        self.genome = Genome()