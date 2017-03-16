import csv, requests, sys, json, re

# Usage: python other_variants_checker.py filename.csv (where filename.csv is a comma separated value file with headings: 'VARIANT_ID', 'CDNACHANGE', 
# 'AACHANGE', 'GENE, 'REFSEQ', 'ENSEMBL_TRANSCRIPT', 'ENSEMBL_GENE_ID', 'VARIANTSTATUS', 'ZYGOSITY', 'EXONINTRON', 'PROTEINEFFECT', 'MOI', 
# 'VARIANTCATEGORY', 'CDNA_START', 'DIST_FROM_CDNA', 'SEQUENCE'

file = sys.argv[1] # filename.csv
variants = []

class Variant:

# Define variables for data from input filename.csv

    def __init__(self, variant_id, cdnachange, aachange, gene, refseq, ensembl_transcript, ensembl_gene_id, variant_status, zygosity, exon_intron, protein_effect, moi, category, cdna_start, dist_from_cdna_start, cdna_end, dist_from_cdna_end, sequence):
        self.variant_id = variant_id
        self.cdnachange = cdnachange
        self.aachange = aachange
        self.gene = gene
        self.refseq = refseq
        self.ensembl_transcript = ensembl_transcript
        self.ensembl_gene_id = ensembl_gene_id
        self.status = variant_status
        self.zygosity = zygosity
        self.exon_intron = exon_intron
        self.protein_effect = protein_effect
        self.mode = moi
        self.category = category
        self.cdna_start = cdna_start
        self.dist_from_cdna_start = dist_from_cdna_start
        self.cdna_end = cdna_end
        self.dist_from_cdna_end = dist_from_cdna_end
        self.sequence = sequence


    def ensembl_capture(self,transcript,cdna_pos):
        '''
        Generate genomic co-ordinate from ensembl REST API using cDNA and transcript
        '''
        cdna_input = str(cdna_pos) + ".." + str(cdna_pos) + "?"
        server = "http://grch37.rest.ensembl.org"
        ext = "/map/cds/"+ transcript + "/" + cdna_input

        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        self.cdna_pos_coord = decoded['mappings'][0]['start']
        self.chr = decoded['mappings'][0]['seq_region_name']
        self.strand = decoded['mappings'][0]['strand']

    def cdna_start_coord_capture(self, transcript, cdna_start):
        ''' 
        Use the cdna_start position to retrieve the genomic co-ordinate
        '''
        self.ensembl_capture(transcript,cdna_start)
        self.cdna_start_coord = self.cdna_pos_coord

    def cdna_end_coord_capture(self, transcript, cdna_end):
        '''
        Use cdna_end co-ord (where it exists, skip otherwise) to retrieve the genomic co-ordinate
        '''
        if cdna_end == 'N/A':
            self.cdna_end_coord = 'N/A'
        else:
            self.ensembl_capture(transcript,cdna_end)
            self.cdna_end_coord = self.cdna_pos_coord


    def working_coords_capture(self,dist_from_cdna_start, dist_from_cdna_end):
         '''
         Generate the genomic co-ordinate for the starting cDNA position, start postion for ExAC database (one base before the cDNA position for 
         all variants that are not missense) and end position of variation
         '''
         if self.strand == 1:
            self.dna_start_pos = int(self.cdna_start_coord)+int(dist_from_cdna_start)
            print("DNA start pos = " + str(self.dna_start_pos))
            if self.cdna_end_coord == 'N/A':
                self.dna_end_pos = 'N/A'
            else: 
                self.dna_end_pos = int(self.cdna_end_coord)+int(dist_from_cdna_end)

         elif self.strand == -1:
            self.dna_start_pos = int(self.cdna_start_coord)-int(dist_from_cdna_start)        
            if self.cdna_end_coord == 'N/A':
                self.dna_end_pos = 'N/A'
            else: 
                self.dna_end_pos = int(self.cdna_end_coord)-int(dist_from_cdna_end)
         else:
             print("No strand information found")


         if self.dna_end_pos != 'N/A':
              if int(self.dna_start_pos) > int(self.dna_end_pos):
                  self.exac_start_pos = int(self.dna_end_pos)-1
              else:
                  self.exac_start_pos = int(self.dna_start_pos)-1
         else:    
              self.exac_start_pos = int(self.dna_start_pos)-1  
            
         print("Strand: " + str(self.strand),"\nDNA start pos: " + str(self.dna_start_pos) + "\nExAC start pos: " + str(self.exac_start_pos) + "\nDNA end pos: " + str(self.dna_end_pos))


    def base_capture(self,chr,base_pos):
        '''
        Capture the base at a specified genomic position      
        '''
        base_input = str(chr) + ":" + str(base_pos) + ".." + str(base_pos) + ":1?"

        server = "http://grch37.rest.ensembl.org"
        ext = "/sequence/region/human/" + base_input

        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()
 
        self.base = r.text

    def dna_start_base(self):
        '''
        Capture the base at the dna_start_pos position i.e. the first base of the variation      
        '''
        self.base_capture(self.chr,self.dna_start_pos)
        self.dna_start_base = self.base

    def dna_end_base(self):
        '''
        Capture the base at the dna_end_pos position i.e. the last base of the variation      
        '''
        if self.dna_end_pos == 'N/A':
            self.dna_end_base = 'N/A'
        else:
            self.base_capture(self.chr,self.dna_end_pos)
            self.dna_end_base = self.base

    def exac_start_base(self):
        '''
        Capture the base at the exac_start_pos position i.e. the first base of the variation that ExAC uses in its' nomenclature      
        '''
        self.base_capture(self.chr,self.exac_start_pos)
        self.exac_start_base = self.base

    def sequence_order(self):
        '''
        First genomic co-ordinate in the sequence_capture request must be a smaller number than the second genomic co-ordinate or it will raise an error therfroe the cdna_start and cdna_end must be ordered correctly
        '''
        if self.dna_end_pos != 'N/A':

            if self.dna_start_pos > self.dna_end_pos:
                self.start = self.dna_end_pos
                self.end =  self.dna_start_pos
            elif self.dna_start_pos < self.dna_end_pos:
                self.start = self.dna_start_pos
                self.end = self.dna_end_pos
            elif self.dna_start_pos == self.dna_end_pos:
                self.start = self.dna_start_pos
                self.end = self.dna_start_pos
            else:
                print("Strand not found.\nUnable to execute sequence_order function.\nProgram will abort.")
      
        else:
            # dna_end_pos = 'N/A' therefore the variant is a single nucleotide variant
            pass

    def sequence_capture(self):
        '''
        Capture sequences from Ensembl between two bases (the start and end positions of the cDNA, inclusive)
        '''
        self.sequence_order()
        server = "http://grch37.rest.ensembl.org"
        ext = "/sequence/region/human/" + str(self.chr) + ":" + str(self.start) + ".." + str(self.end) + ":1?"
        r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        self.seq_capture = r.text

#    def rev_comp(self):
#        '''
#        Create the complmentary DNA sequence
#        '''
#
#        self.sequence_capture()

#        comp_base_dict = {'A':'T','C':'G','G':'C','T':'A'}
#        print(self.seq_capture)
#        seq_list = list(self.seq_capture)
#        print(seq_list)
#        seq_list = [comp_base_dict[base] for base in seq_list]
#        self.seq_capture_rev = ''.join(seq_list)
#        print(self.seq_capture_rev)


    def wt_allele_construct(self,category):
        '''
        Retrieve and/or construct the WT allele sequence
        '''
        if category == "SN deletion":
            # Base at ExAC start position + base at cDNA start          
            self.wt_allele_seq = self.exac_start_base + self.dna_start_base

        elif category == "SN duplication":
            # Base at ExAC start position
            self.wt_allele_seq = self.exac_start_base
            
        elif category == ">1bp deletion":
            self.sequence_capture()
            # Base at ExAC start position + DNA sequence between cDNA start and cDNA end inclusive            
#            self.rev_comp()
            self.wt_allele_seq = self.exac_start_base + self.seq_capture

        elif category == ">1bp duplication":
            # Base at ExAC start position
            self.wt_allele_seq = self.exac_start_base
            
        elif category == "delins":
            self.sequence_capture()
#            self.rev_comp()
            # Base at ExAC start position + DNA sequence between cDNA start and cDNA end inclusive
            self.wt_allele_seq = self.exac_start_base + self.seq_capture

        elif category == "insertion":
            # Base at ExAC start position
            self.wt_allele_seq = self.exac_start_base

        else:
            print("Variant category could not be identified.\nProgram will abort")
            exit()


    def alt_allele_construct(self,category,sequence):
        '''
        Retrieve and/or construct the alt allele sequence
        '''
        if category == "SN deletion":
            # Base at ExAC start position
            self.alt_allele_seq = self.exac_start_base
       
        elif category == "SN duplication":
            # Base at ExAC start position + base at first genomic position of cDNA
            self.alt_allele_seq = self.exac_start_base + self.dna_start_base

        elif category == ">1bp deletion":
            # Base at ExAC start position
            self.alt_allele_seq = self.exac_start_base

        elif category == ">1bp duplication":
            # Base at ExAC start position + DNA sequence between cDNA start and cDNA end inclusive
            self.sequence_capture()
#            self.rev_comp()
            self.alt_allele_seq = self.exac_start_base + self.seq_capture

        elif category == "delins":
            # Base at ExAC start position + DNA sequence from input file          
            self.alt_allele_seq = self.exac_start_base + self.sequence

        elif category == "insertion":
            # Base at ExAC start position + DNA sequence from input file
            self.alt_allele_seq = self.exac_start_base + self.sequence

        else:
            print("Variant category could not be identified.\nProgram will abort")
            exit()


    def exac_capture(self):
        '''
        Generate ExAC API list, submit and receive output
        '''
        
        server = "http://exac.hms.harvard.edu"
        ext = "/rest/variant/variant/" + str(self.chr) + "-" + str(self.exac_start_pos) + "-" + str(self.wt_allele_seq) + "-" + str(self.alt_allele_seq)
        print(ext)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        
        if not r.ok:
             r.raise_for_status()
             sys.exit()

        data = r.json()
        
        try:
            self.hom_NFE = data["pop_homs"]["European (Non-Finnish)"]
            self.hom_EA = data["pop_homs"]["East Asian"]
            self.hom_Other = data["pop_homs"]["Other"]
            self.hom_AFR = data["pop_homs"]["African"]
            self.hom_LAT = data["pop_homs"]["Latino"]
            self.hom_SA = data["pop_homs"]["South Asian"]
            self.hom_FE = data["pop_homs"]["European (Finnish)"]
            self.tot_hom = data["hom_count"]

            self.pop_NFE = data["pop_acs"]["European (Non-Finnish)"]
            self.pop_EA = data["pop_acs"]["East Asian"]
            self.pop_Other = data["pop_acs"]["Other"]
            self.pop_AFR = data["pop_acs"]["African"]
            self.pop_LAT = data["pop_acs"]["Latino"]
            self.pop_SA = data["pop_acs"]["South Asian"]
            self.pop_FE = data["pop_acs"]["European (Finnish)"]
            self.allele_count = data["allele_count"]

            self.pop_total_NFE = data["pop_ans"]["European (Non-Finnish)"]
            self.pop_total_EA = data["pop_ans"]["East Asian"]
            self.pop_total_Other = data["pop_ans"]["Other"]
            self.pop_total_AFR = data["pop_ans"]["African"]
            self.pop_total_LAT = data["pop_ans"]["Latino"]
            self.pop_total_SA = data["pop_ans"]["South Asian"]
            self.pop_total_FE = data["pop_ans"]["European (Finnish)"]
            self.tot_allele_num = data["allele_num"]

            self.freq_NFE = float(data["pop_acs"]["European (Non-Finnish)"])/data["pop_ans"]["European (Non-Finnish)"]
            self.freq_EA = float(data["pop_acs"]["East Asian"])/data["pop_ans"]["East Asian"]
            self.freq_Other = float(data["pop_acs"]["Other"])/data["pop_ans"]["Other"]
            self.freq_AFR = float(data["pop_acs"]["African"])/data["pop_ans"]["African"]
            self.freq_LAT = float(data["pop_acs"]["Latino"])/data["pop_ans"]["Latino"]
            self.freq_SA = float(data["pop_acs"]["South Asian"])/data["pop_ans"]["South Asian"]
            self.freq_FE = float(data["pop_acs"]["European (Finnish)"])/data["pop_ans"]["European (Finnish)"]

            self.total_allele_freq = data["allele_freq"]

#            self.genomic_start_calc = self.genomic_start

        except:
            self.hom_NFE = "0"
            self.hom_EA = "0"
            self.hom_Other = "0"
            self.hom_AFR = "0"
            self.hom_LAT = "0"
            self.hom_SA = "0"
            self.hom_FE = "0"
            self.tot_hom = "0"

            self.pop_NFE = "0"
            self.pop_EA = "0"
            self.pop_Other = "0"
            self.pop_AFR = "0"
            self.pop_LAT = "0"
            self.pop_SA = "0"
            self.pop_FE = "0"
            self.allele_count = "0"

            self.pop_total_NFE = "0"
            self.pop_total_EA = "0"
            self.pop_total_Other = "0"
            self.pop_total_AFR = "0"
            self.pop_total_LAT = "0"
            self.pop_total_SA = "0"
            self.pop_total_FE = "0"
            self.tot_allele_num = "0"

            self.freq_NFE = "0"
            self.freq_EA = "0"
            self.freq_Other = "0"
            self.freq_AFR = "0"
            self.freq_LAT = "0"
            self.freq_SA = "0"
            self.freq_FE = "0"

            self.total_allele_freq = "0"


# Open custom report file provided
with open(file) as infile:
    csvreader = csv.reader(infile, delimiter = ',')
    # Skip header row
    next(infile)
    # Add data to Variant class
    for row in csvreader:
        print("Converting row " + row[0])
        row[0] = Variant(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17]) 
        row[0].cdna_start_coord_capture(row[5],row[13]) # uses cDNA start postion to retrieve genomic co-ordinate for start of sequence
        row[0].cdna_end_coord_capture(row[5],row[15]) # uses cDNA end postion to retrieve genomic co-ordinate for end of sequence
        row[0].working_coords_capture(row[14],row[16]) # use intronic distances form start and end cdna positions to calculate genomic start and end and exac start positions
        row[0].dna_start_base()
        row[0].dna_end_base()
        row[0].exac_start_base()
        row[0].wt_allele_construct(row[12])
        row[0].alt_allele_construct(row[12],row[17])
        row[0].exac_capture()
        variants.append(row[0])

# Open outfile for writing successfully retrieved variants from ExAC
with open(file+'_output.csv', 'w') as outfile:
    csvwriter = csv.writer(outfile, delimiter=',')
    # Add a header row
    csvwriter.writerow(["VARIANT_ID","CDNA_CHANGE", "AACHANGE","GENE", "REFSEQ", "ENSEMBL_TRANSCRIPT", "ENSEMBL_GENE_ID", "VARIANT_STATUS", "ZYGOSITY", "EXON_INTRON", "PROTEIN_EFFECT", "MOI", "VARIANT_CATEGORY", "CDNA_START","DIST_FROM_CDNA_START","CDNA_END","DIST_FROM_CDNA_END","SEQUENCE","CHROMOSOME","STRAND", "CDNA_START_POS", "DNA_START_BASE", "EXAC_START_POS", "EXAC_START_BASE", "CDNA_END_POS", "DNA_END_BASE", "WT_ALLELE_SEQ", "ALT_ALLELE_SEQ", "total_allele_count", "tot_allele_num","tot_hom", "tot_allele_freq", "hom_AFR", "hom_NFE", "hom_EA", "hom_FE", "hom_LAT", "hom_Other", "hom_SA", "pop_AFR", "pop_LAT", "pop_Other", "pop_SA","pop_total_AFR", "pop_total_NFE", "pop_total_EA", "pop_total_FE", "pop_total_LAT", "pop_total_Other","pop_total_SA","freq_AFR", "freq_NFE", "freq_EA", "freq_FE", "freq_LAT", "freq_Other", "freq_SA"])
    # Loop through created variant objects
    for variant in variants:
        row = [variant.variant_id, variant.cdnachange, variant.aachange, variant.gene, variant.refseq, variant.ensembl_transcript, variant.ensembl_gene_id, variant.status, variant.zygosity, variant.exon_intron, variant.protein_effect, variant.mode, variant.category, variant.cdna_start, variant.dist_from_cdna_start, variant.cdna_end, variant.dist_from_cdna_end, variant.sequence, variant.chr, variant.strand, variant.dna_start_pos, variant.dna_start_base, variant.exac_start_pos, variant.exac_start_base, variant.dna_end_pos, variant.dna_end_base, variant.wt_allele_seq, variant.alt_allele_seq, variant.allele_count, variant.tot_allele_num, variant.tot_hom, variant.total_allele_freq, variant.hom_AFR, variant.hom_NFE, variant.hom_EA, variant.hom_FE, variant.hom_SA, variant.pop_AFR, variant.pop_NFE, variant.pop_EA, variant.pop_FE, variant.pop_LAT, variant.pop_Other, variant.pop_SA, variant.pop_total_AFR, variant.pop_total_NFE, variant.pop_total_EA, variant.pop_total_FE, variant.pop_total_LAT, variant.pop_total_Other, variant.pop_total_SA, variant.freq_AFR, variant.freq_NFE, variant.freq_EA, variant.freq_FE, variant.freq_LAT, variant.freq_Other, variant.freq_SA]
        # Write row to outfile
        csvwriter.writerow(row)
