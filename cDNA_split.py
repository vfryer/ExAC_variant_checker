# usage python cDNA_split.py input_filename output_filename

import re, csv, sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

# for initial test purposes only
# cDNA = "c.1468T>C"
# print("cDNA is:" + cDNA)

# open a .csv and write to a new .csv appending new data to the end of each line
with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
    reader = csv.reader(input_file)
    # skip headers
    headers = next(reader, None)
    writer = csv.writer(output_file)

    # write headers from input file into output file
    if headers:
        writer.writerow(headers)

    for row in reader:

        if '>' and 'c.0' in row[1]:
            writer.writerow(row+['Out of exome range'])
        elif '>' in row[1]:
            writer.writerow(row+['SNV'])
        elif '_+' in row[1]:
            writer.writerow(row+['Variation nomenclature incorrect'])
        elif '?' in row[1]:
            writer.writerow(row+['Unknown start/end position of variation'])
        elif 'c.(' in row[1]:
            writer.writerow(row+['Unknown start/end position of variation'])
        elif 'delins' in row[1]:
            writer.writerow(row+['delins'])
        elif 'het' in row[1]:
            cDNA = row[1]
            start = 'c.' 
            end = 'het'
            pos = ((cDNA.split(start))[1].split(end)[0])
            if 'del' in row[1]:
                row[1] = start + pos + 'del'
                writer.writerow(row+['SN deletion'])
            elif 'dup' in row[1]:
                row[1] = start + pos + 'dup'
                writer.writerow(row+['SN duplication'])
 
        elif '_' in row[1]:
            cDNA = row[1]
            start = 'c.'
            middle = '_'

            if 'del' in row[1]:
                end = 'del'
                first_pos = ((cDNA.split(start))[1].split(middle)[0])
                second_pos = ((cDNA.split(middle))[1].split(end)[0])
                if first_pos == second_pos:
                    row[1] = start + second_pos + end
                    writer.writerow(row+['SN deletion'])
                else:
                    writer.writerow(row+['>1bp deletion'])
            elif 'dup' in row[1]:         
                end = 'dup'
                first_pos = ((cDNA.split(start))[1].split(middle)[0])
                second_pos = ((cDNA.split(middle))[1].split(end)[0])
                if first_pos == second_pos:
                    row[1] = start + second_pos + end
                    writer.writerow(row+['SN duplication'])
                else:
                    writer.writerow(row+['>1bp duplication'])
            elif 'ins' in row[1]:         
                writer.writerow(row+['insertion'])
            else:
                writer.writerow(row+[0])
        elif 'del' in row[1]:
            writer.writerow(row+['SN deletion'])
        elif 'dup' in row[1]:
            writer.writerow(row+['SN duplication'])
        elif '' in row[1]:
            writer.writerow(row+['No variation to assess'])
        else:
            writer.writerow(row+['Variation nomenclature not understood'])

#cDNA_pos = map(int,re.findall(r'\d+',cDNA))
#print("cDNA position : " + str(cDNA_pos))

#wt_allele = cDNA.split('>',1)[0][-1]
#print("WT Allele: " + wt_allele)

#alt_allele = cDNA.split('>',1)[1]
#'input_file.txt'print("Alt allele: " + alt_allele)

