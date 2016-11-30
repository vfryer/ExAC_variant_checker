import csv, requests, sys, json

# Usage: python exac_variant_finder.py filename.csv (where filename.csv is a comma separated value file with headings:
# 'Primary Key, Ensembl Transcript Number, cDNA position (from STARLiMS, without c. or nucleotide change), CDS WT allele, CDS alt allele, WT allele from STARLiMS, alt allele from STARLiMS, amin acid change, zygosity, variant status, protein effect, exon/intron, exonintron, DMUDBEXP, mode of inheritance, gene, STARLiMS RefSeq, cDNA change, cDNA position)

file = sys.argv[1]
variants = []

class Variant:

# Define variables for data from input filename.csv

    def __init__(self,gene,transcript,cDNA,wt_allele,alt_allele,zygosity,variant_status,mode,protein_effect):
        self.gene = gene
        self.transcript = transcript
        self.cDNA = cDNA
	self.wt_allele = wt_allele
	self.alt_allele = alt_allele
        self.zygosity = zygosity
        self.variant_status = variant_status
        self.mode = mode
        self.protein_effect = protein_effect

# Generate genomic co-ordinates from ensembl REST API using cDNA and transcript and add to list of variants

    def ensembl_capture(self,transcript,cDNA):
        cDNA_input = cDNA + ".." + cDNA + "?"
        server = "http://grch37.rest.ensembl.org"
        ext = "/map/cds/"+ transcript + "/" + cDNA_input

        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        self.genomic_start = decoded['mappings'][0]['start']
        self.genomic_end = decoded['mappings'][0]['end']
        self.chromosome = decoded['mappings'][0]['seq_region_name']
        self.strand = decoded['mappings'][0]['strand']

# Generate ExAC API list, submit and receive output
    def exac_capture(self):
        server = "http://exac.hms.harvard.edu"
        ext = "/rest/variant/variant/" + str(self.chromosome) + "-" + str(self.genomic_start) + "-" + str(self.wt_allele) + "-" + str(self.alt_allele)
        
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
            self.total_alleles = data["allele_num"]

            self.freq_NFE = float(data["pop_acs"]["European (Non-Finnish)"])/data["pop_ans"]["European (Non-Finnish)"]
            self.freq_EA = float(data["pop_acs"]["East Asian"])/data["pop_ans"]["East Asian"]
            self.freq_Other = float(data["pop_acs"]["Other"])/data["pop_ans"]["Other"]
            self.freq_AFR = float(data["pop_acs"]["African"])/data["pop_ans"]["African"]
            self.freq_LAT = float(data["pop_acs"]["Latino"])/data["pop_ans"]["Latino"]
            self.freq_SA = float(data["pop_acs"]["South Asian"])/data["pop_ans"]["South Asian"]
            self.freq_FE = float(data["pop_acs"]["European (Finnish)"])/data["pop_ans"]["European (Finnish)"]
            self.total_allele_freq = data["allele_freq"]

        except:
            self.hom_NFE = "0"
            self.hom_EA = "0"
            self.hom_Other = "0"
            self.hom_AFR = "0"
            self.hom_LAT = "0"
            self.hom_SA = "0"
            self.hom_FE = "0"

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
            self.total_alleles = "0"

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
        print "Converting row " + row[0]
        row[0] = Variant(row[15], row[1], row[2], row[3], row[4], row[8], row[9], row[14], row[10])
        row[0].ensembl_capture(row[1],row[2])
        row[0].exac_capture()
        # Add amplicon object name to amplicons array
        variants.append(row[0])

 
# Open outfile for writing successfully retrieved variants from ExAC
with open('output.csv', 'wb') as outfile:
    csvwriter = csv.writer(outfile, delimiter=',')
    # Add required header row
    csvwriter.writerow(["gene", "transcript", "cDNA", "chromosome", "genomic_start", "genomic_end", "strand", "wt_allele", "alt_allele", "zygosity", "variant_status", "mode", "protein_effect", "allele_count", "tot_allele_freq", "total_alleles", "hom_NFE", "hom_EA", "hom_Other", "hom_AFR", "hom_LAT", "hom_SA", "hom_FE", "pop_NFE", "pop_EA", "pop_Other", "pop_AFR", "pop_LAT", "pop_SA", "pop_FE" ,"pop_total_NFE", "pop_total_EA", "pop_total_Other", "pop_total_AFR", "pop_total_LAT", "pop_total_SA", "pop_total_FE", "freq_NFE", "freq_EA", "freq_Other", "freq_AFR", "freq_LAT", "freq_SA", "freq_FE"]) 
    # Cycle through created variant objects
    for variant in variants:
        row = [variant.gene, variant.transcript, variant.cDNA, variant.chromosome, variant.genomic_start, variant.genomic_end, variant.strand, variant.wt_allele, variant.alt_allele, variant.zygosity, variant.variant_status, variant.mode, variant.protein_effect, variant.allele_count, variant.total_allele_freq, variant.total_alleles, variant.hom_NFE, variant.hom_EA, variant.hom_Other, variant.hom_AFR, variant.hom_LAT, variant.hom_SA, variant.hom_FE, variant.pop_NFE, variant.pop_EA, variant.pop_Other, variant.pop_AFR, variant.pop_LAT, variant.pop_SA, variant.pop_FE, variant.pop_total_NFE, variant.pop_total_EA, variant.pop_total_Other, variant.pop_total_AFR, variant.pop_total_LAT, variant.pop_total_SA, variant.pop_total_FE, variant.freq_NFE, variant.freq_EA, variant.freq_Other, variant.freq_AFR, variant.freq_LAT, variant.freq_SA, variant.freq_FE]
	# Write row to outfile
        csvwriter.writerow(row)
