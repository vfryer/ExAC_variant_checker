import csv, requests, sys, json, re

# Usage: python exac_variant_finder.py filename.csv (where filename.csv is a comma separated value file with headings:
# 'Primary Key, Ensembl Transcript Number, cDNA position (from STARLiMS, without c. or nucleotide change), CDS WT allele, CDS alt allele, WT allele from STARLiMS, alt allele from STARLiMS, amin acid change, zygosity, variant status, protein effect, exon/intron, exonintron, DMUDBEXP, mode of inheritance, gene, STARLiMS RefSeq, cDNA change, cDNA position)

file = sys.argv[1]
variants = []

class Variant:

# Define variables for data from input filename.csv

    def __init__(self,gene,transcript,cDNA_change,cDNA_pos,cDNA_add,wt_allele,alt_allele,aa_change,refseq_id,ensembl_gene_id,zygosity,variant_status,mode,protein_effect, exon_intron):
        self.gene = gene
        self.transcript = transcript
        self.cDNA_change = cDNA_change
        self.cDNA_pos = cDNA_pos
        self.cDNA_add = cDNA_add
        self.wt_allele = wt_allele
        self.alt_allele= alt_allele
        self.aa_change = aa_change
        self.refseq_id = refseq_id
        self.ensembl_gene_id = ensembl_gene_id
        self.zygosity = zygosity
        self.variant_status = variant_status
        self.mode = mode
        self.protein_effect = protein_effect
        self.exon_intron = exon_intron

# Generate genomic co-ordinates from ensembl REST API using cDNA and transcript and add to list of variants

    def ensembl_capture(self,transcript,cDNA_pos):
        cDNA_input = str(cDNA_pos) + ".." + str(cDNA_pos) + "?"
        server = "http://grch37.rest.ensembl.org"
        ext = "/map/cds/"+ transcript + "/" + cDNA_input

        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        self.genomic_start = decoded['mappings'][0]['start']
        self.chromosome = decoded['mappings'][0]['seq_region_name']
        self.strand = decoded['mappings'][0]['strand']


# Generate ExAC API list, submit and receive output
    def exac_capture(self, wt_allele, alt_allele):
        server = "http://exac.hms.harvard.edu"
        ext = "/rest/variant/variant/" + str(self.chromosome) + "-" + str(self.genomic_start) + "-" + wt_allele + "-" + alt_allele
        
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

            self.genomic_start_calc = self.genomic_start

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

            self.genomic_start_calc = self.genomic_start

# Generate ExAC API list for a second time with calculated genomic start, submit and receive output
    def exac_capture_calc(self, genomic_start_calc, wt_allele, alt_allele):
        server = "http://exac.hms.harvard.edu"
        ext = "/rest/variant/variant/" + str(self.chromosome) + "-" + str(genomic_start_calc) + "-" + wt_allele + "-" + alt_allele
        
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

            self.genomic_start_calc = str(genomic_start_calc)

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

            self.genomic_start_calc = str(genomic_start_calc)

# Open custom report file provided
with open(file) as infile:
    csvreader = csv.reader(infile, delimiter = ',')
    # Skip header row
    next(infile)
    # Add data to Variant class
    for row in csvreader:
        print "Converting row " + row[0]
        row[0] = Variant(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15])
        row[0].ensembl_capture(row[2],row[4])
        if row[4] != row[5]: 
            if row[0].strand == 1:
                #print "forward strand"
                #print "strand +1 orig genomic start = " + str(row[0].genomic_start)
                genomic_start_calc = int(row[0].genomic_start) + int(row[5])
                #print genomic_start_calc
            elif row[0].strand == -1:
                #print "reverse strand"
                #print "strand -1 orig genomic start = " + str(row[0].genomic_start)
                genomic_start_calc = int(row[0].genomic_start) - int(row[5])
                #print genomic_start_calc
            # Add amplicon object name to amplicons array
            row[0].exac_capture_calc(genomic_start_calc,row[6],row[7])
            variants.append(row[0])
        elif row[4] == row[5]:
            row[0].ensembl_capture(row[2],row[4])
            row[0].exac_capture(row[6],row[7])
            variants.append(row[0])
        else:
            print "Out of loop"
            row[0].ensembl_capture(row[2],row[4])
            row[0].exac_capture(row[6],row[7])
            variants.append(row[0])

# Open outfile for writing successfully retrieved variants from ExAC
with open('output.csv', 'wb') as outfile:
    csvwriter = csv.writer(outfile, delimiter=',')
    # Add required header row
    csvwriter.writerow(["gene", "transcript", "cDNA_change", "cDNA_pos", "cDNA_add", "aa_change", "chromosome", "genomic_start", "genomic_start_calc", "strand", "wt_allele", "alt_allele", "zygosity", "variant_status", "mode", "refseq_id", "ensembl_gene_id", "exon_intron", "protein_effect", "total_allele_count", "tot_allele_num", "tot_hom", "tot_allele_freq", "hom_AFR", "hom_NFE", "hom_EA", "hom_FE", "hom_LAT", "hom_Other", "hom_SA", "pop_AFR", "pop_NFE", "pop_EA", "pop_FE", "pop_LAT", "pop_Other", "pop_SA","pop_total_AFR", "pop_total_NFE", "pop_total_EA", "pop_total_FE", "pop_total_LAT", "pop_total_Other","pop_total_SA", "freq_AFR", "freq_NFE", "freq_EA", "freq_FE", "freq_LAT", "freq_Other", "freq_SA"]) 
    # Cycle through created variant objects
    for variant in variants:
        row = [variant.gene, variant.transcript, variant.cDNA_change, variant.cDNA_pos, variant.cDNA_add, variant.aa_change, variant.chromosome, variant.genomic_start, variant.genomic_start_calc, variant.strand, variant.wt_allele, variant.alt_allele, variant.zygosity, variant.variant_status, variant.mode, variant.refseq_id, variant. ensembl_gene_id, variant.exon_intron, variant.protein_effect, variant.allele_count, variant.tot_allele_num, variant.tot_hom, variant.total_allele_freq, variant.hom_AFR, variant.hom_NFE, variant.hom_EA, variant.hom_FE, variant.hom_LAT, variant.hom_Other, variant.hom_SA, variant.pop_AFR, variant.pop_NFE, variant.pop_EA, variant.pop_FE, variant.pop_LAT, variant.pop_Other, variant.pop_SA, variant.pop_total_AFR, variant.pop_total_NFE, variant.pop_total_EA, variant.pop_total_FE, variant.pop_total_LAT, variant.pop_total_Other, variant.pop_total_SA, variant.freq_AFR, variant.freq_NFE, variant.freq_EA, variant.freq_FE, variant.freq_LAT, variant.freq_Other, variant.freq_SA]
	# Write row to outfile
        csvwriter.writerow(row)
