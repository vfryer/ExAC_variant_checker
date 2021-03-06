# usage python cDNA_split.py input_filename
# Parses file with data from STARLiMS "Variant_ID","cDNA_change","AA_change","Gene","RefSeq","Ensembl_transcript","Ensembl_gene_ID","Variant_status","Zygosity","Exon_Intron","Protein_effect","MOI"
# Outputs only variants that can be searched in ExAC (does not include deepintornic variants, variants >50bp from cDNA of those otherwise incorrectly named or not possible to parse)

import re, csv, sys

input_filename = sys.argv[1]

# open a .csv and write to a new .csv appending new data to the end of each line
with open(input_filename, 'r') as input_file, open(input_filename+'_parser_output.csv', 'w') as output_file:
    reader = csv.reader(input_file)
    # skip headers
    headers = next(reader, None)
    writer = csv.writer(output_file)

    # write headers from input file into output file
    if headers:
        writer.writerow(headers+['VARIANTCATEGORY']+['CDNA_START']+['DIST_FROM_CDNA_START']+['CDNA_END']+['DIST_FROM_CDNA_END']+['SEQUENCE'])

    for row in reader:

        cDNA = row[1]
        start = 'c.'
        print('Processing row ' + row[0])


        if '>' and 'c.0' in row[1]:
            pass

        elif 'c.-' in row[1]:
            pass

        elif '>' in row[1]:
            pos = row[1][2:-3]
            alt_allele = row[1][-1]

            if '-' in pos:
                # if the SNV is upstream of the nearest exon
                first_pos = pos.split('-')[0]
                inton_pos = pos.split('-')[1]

                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) > 200000:
                    pass
                else:
                    writer.writerow(row+['SNV'] + [first_pos] + ['-'+intron_pos] + ['N/A'] + ['0'] + [alt_allele])            
            elif '+' in pos:
                # if the SNV is downstream of the nearest exon
                first_pos = pos.split('+')[0]
                intron_pos = pos.split('+')[1]
                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) > 200000:
                    pass
                else:
                    writer.writerow(row+['SNV'] + [first_pos] + [intron_pos] + ['N/A'] + ['0'] + [alt_allele])            
            elif int(pos) >200000:
                pass
            else:
                writer.writerow(row+['SNV'] + [pos] + ['0'] + ['N/A'] + ['0'] + [alt_allele])            

        elif '_+' in row[1]:
            pass

        elif '?' in row[1]:
            pass

        elif 'c.(' in row[1]:
            pass

        elif 'delins' in row[1]:
            middle = '_'
            pos = ((cDNA.split(start))[1].split(middle)[0])
            second_pos = ((cDNA.split(middle))[1].split('delins')[0])
            sequence = cDNA.split('delins')[-1]

            if '-' in pos:
                first_pos = pos.split('-')[0]
                intron_pos = pos.split('-')[-1]
                if '-' in second_pos:
                    second_pos_split = second_pos.split('-')[0]
                    second_pos_intron = second_pos.split('-')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+['-'+intron_pos]+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                elif '+' in second_pos:
                    second_pos_split = second_pos.split('+')[0]
                    second_pos_intron = second_pos.split('+')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+['-'+intron_pos]+[second_pos_split]+[second_pos_intron]+[sequence])

                else:
                    # second position is exonic
                    if int(intron_pos) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+['-'+intron_pos]+[second_pos]+['0']+[sequence])

            elif '+' in pos:
                first_pos = pos.split('+')[0]
                intron_pos = pos.split('+')[-1]
                if '-' in second_pos:
                    second_pos_split = second_pos.split('-')[0]
                    second_pos_intron = second_pos.split('-')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+[intron_pos]+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                elif '+' in second_pos:
                    second_pos_split = second_pos.split('+')[0]
                    second_pos_intron = second_pos.split('+')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+[intron_pos]+[second_pos_split]+[second_pos_intron]+[sequence])
                else:
                    # second position is exonic
                    if int(intron_pos) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[first_pos]+['-'+intron_pos]+[second_pos]+['0']+[sequence])
            else:
                # the first postion is exonic
                if '-' in second_pos:
                    second_pos_split = second_pos.split('-')[0]
                    second_pos_intron = second_pos.split('-')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[pos]+['0']+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                elif '+' in second_pos:
                    second_pos_split = second_pos.split('+')[0]
                    second_pos_intron = second_pos.split('+')[-1]
                    if int(intron_pos) >=51:
                        pass
                    elif int(second_pos_intron) >=51:
                        pass
                    else:
                        writer.writerow(row+['delins']+[pos]+['0']+[second_pos_split]+[second_pos_intron]+[sequence])

                else:
                    # second position is exonic
                    writer.writerow(row+['delins']+[pos]+['0']+[second_pos]+['0']+[sequence])
                 
                
        elif 'het' in row[1]:
            end = 'het'
            pos = ((cDNA.split(start))[1].split(end)[0])
            if 'del' in row[1]:
                row[1] = start + pos + 'del'
                writer.writerow(row+['SN deletion']+[pos]+['0']+['N/A']+['0']+['N/A'])
            elif 'dup' in row[1]:
                row[1] = start + pos + 'dup'
                writer.writerow(row+['SN duplication']+[pos]+['0']+['N/A']+['0']+['N/A'])
 
        elif '_' in row[1]:
            middle = '_'

            if 'del' in row[1]:
                end = 'del'
                pos = ((cDNA.split(start))[1].split(middle)[0])
                second_pos = ((cDNA.split(middle))[1].split(end)[0])
                if pos == second_pos:
                    # if the first and second position are the same
                    if '-' in pos:
                        first_pos = pos.split('-')[0]
                        intron_pos = pos.split('-')[-1]
                        row[1] = start + second_pos + end
                        if int(intron_pos) >=51:
                            pass
                        elif int(first_pos) >=200000:
                            pass
                        else:
                            writer.writerow(row+['SN deletion']+[first_pos]+['-'+intron_pos]+['N/A']+['0']+['N/A'])
                    elif '+' in pos:
                        first_pos = pos.split('+')[0]
                        intron_pos = pos.split('+')[-1]
                        row[1] = start + second_pos + end
                        if int(intron_pos) >=51:
                            pass   
                        elif int(first_pos) >=200000:
                            pass   
                        else:
                            writer.writerow(row+['SN deletion']+[first_pos]+[intron_pos]+['N/A']+['0']+['N/A'])
                    else:
                        row[1] = start + second_pos + end
                        if int(pos) >=200000:
                            pass   
                        else:
                            writer.writerow(row+['SN deletion']+[pos]+['0']+['N/A']+['0']+['N/A'])
                else:
                    # if the first and second position are not the same
                        # if the second position is intronic
                        if '-' in second_pos:
                            second_pos_split = second_pos.split('-')[0]
                            second_pos_intron = second_pos.split('-')[-1]

                            # if the first position is also intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+['-'+intron_pos]+[second_pos_split]+['-'+second_pos_intron]+['N/A'])

                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+[intron_pos]+[second_pos_split]+['-'+second_pos_intron]+['N/A'])
                           
                            else:
                                # first position is exonic
                                first_pos = pos.split('+')[0]                       
                                if int(first_pos) >=200000:
                                    pass
                                else:
                                    if '-' in second_pos:
                                        writer.writerow(row+['>1bp deletion']+[pos]+['0']+[second_pos_split]+['-'+second_pos_intron]+['N/A'])
                                    else:
                                        writer.writerow(row+['>1bp deletion']+[pos]+['0']+[second_pos_split]+['-'+second_pos_intron]+['N/A'])

                        elif '+' in second_pos:
                            second_pos_split = second_pos.split('+')[0]
                            second_pos_intron = second_pos.split('+')[-1]
                       
                            # if the first position is also intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+['-'+intron_pos]+[second_pos_split]+[second_pos_intron]+['N/A'])

                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+[intron_pos]+[second_pos_split]+[second_pos_intron]+['N/A'])

                            else:
                                # first position is exonic
                                first_pos = pos.split('+')[0]                       
                                if int(first_pos) >=200000:
                                    pass
                                    # writer.writerow(row+['cDNA position too large'])
                                else:
                                    if '-' in second_pos:
                                        writer.writerow(row+['>1bp deletion']+[pos]+['0']+['-'+second_pos_split]+[second_pos_intron]+['N/A'])
                                    else:
                                        writer.writerow(row+['>1bp deletion']+[pos]+['0']+[second_pos_split]+[second_pos_intron]+['N/A'])

                        else:
                        # the second position is exonic (i.e. no '+' or '-')
                            # if the first position is intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+['-'+intron_pos]+[second_pos]+['0']+['N/A'])
                                
                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[first_pos]+[intron_pos]+[second_pos]+['0']+['N/A'])
                                 
                            else:
                            # if the first position is exonic
                                first_pos = pos.split('+')[0]                       
                                if int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp deletion']+[pos]+['0']+[second_pos]+['0']+['N/A'])


            elif 'dup' in row[1]:                
                end = 'dup'
                pos = ((cDNA.split(start))[1].split(middle)[0])
                second_pos = ((cDNA.split(middle))[1].split(end)[0])
                if pos == second_pos:
                    # if the first and second position are the same
                    if '-' in pos:
                        first_pos = pos.split('-')[0]
                        intron_pos = pos.split('-')[-1]
                        row[1] = start + second_pos + end
                        if int(intron_pos) >=51:
                            pass
                        elif int(first_pos) >=200000:
                            pass
                        else:
                            writer.writerow(row+['SN duplication']+[first_pos]+['-'+intron_pos]+['N/A']+['0']+['N/A'])
                    elif '+' in pos:
                        first_pos = pos.split('+')[0]
                        intron_pos = pos.split('+')[-1]
                        row[1] = start + second_pos + end
                        if int(intron_pos) >=51:
                            pass
                        elif int(first_pos) >=200000:
                            pass
                        else:
                            writer.writerow(row+['SN duplication']+[first_pos]+[intron_pos]+['N/A']+['0']+['N/A'])
                    else:
                        if int(pos) >=200000:
                            pass
                        else:
                            row[1] = start + second_pos + end
                            writer.writerow(row+['SN duplication']+[pos]+['0']+['N/A']+['0']+['N/A'])

                else:
                    # if the first and second position are not the same
                        # if the second position is intronic
                        if '-' in second_pos:
                            second_pos_split = second_pos.split('-')[0]
                            second_pos_intron = second_pos.split('-')[-1]
                            # if the first position is also intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+['-'+intron_pos]+[second_pos_split]+['-'+second_pos_intron]+['N/A'])

                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+[intron_pos]+[second_pos_split]+[second_pos_intron]+['N/A'])

                            else:
                                # first position is not intronic
                                first_pos = pos.split('+')[0]
                                if int(first_pos) >=200000:
                                    pass
                                else:
                                    if '-' in second_pos:
                                        writer.writerow(row+['>1bp duplication']+[pos]+['0']+[second_pos_split]+['-'+second_pos_intron]+['N/A'])
                                    else:
                                        writer.writerow(row+['>1bp duplication']+[pos]+['0']+[second_pos_split]+[second_pos_intron]+['N/A'])


                        elif '+' in second_pos:
                            second_pos_split = second_pos.split('+')[0]
                            second_pos_intron = second_pos.split('+')[-1]
                       
                            # if the first position is also intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+['-'+intron_pos]+[second_pos_split]+[second_pos_intron]+['N/A'])

                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+[intron_pos]+[second_pos_split]+[second_pos_intron]+['N/A'])

                            else:
                                # first position is not intronic
                                first_pos = pos.split('+')[0]
                                if int(first_pos) >=200000:
                                    pass
                                else:
                                    if '-' in second_pos:
                                        writer.writerow(row+['>1bp duplication']+[pos]+['0']+['-'+second_pos_split]+[second_pos_intron]+['N/A'])
                                    else:
                                        writer.writerow(row+['>1bp duplication']+[pos]+['0']+[second_pos_split]+[second_pos_intron]+['N/A'])
                        else:
                        # the second position is exonic (i.e. no '+' or '-')
                            # if the first position is intronic
                            if '-' in pos:
                                first_pos = pos.split('-')[0]
                                intron_pos = pos.split('-')[-1]
                                if int(intron_pos) >=51:
                                     pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+['-'+intron_pos]+[second_pos]+['0']+['N/A'])
                                
                            elif '+' in pos:
                                first_pos = pos.split('+')[0]
                                intron_pos = pos.split('+')[-1]
                                if int(intron_pos) >=51:
                                    pass
                                elif int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[first_pos]+[intron_pos]+[second_pos]+['0']+['N/A'])
                                 
                            else:
                            # if the first position is exonic
                                first_pos = pos.split('+')[0]                       
                                if int(first_pos) >=200000:
                                    pass
                                else:
                                    writer.writerow(row+['>1bp duplication']+[pos]+['0']+[second_pos]+['0']+['N/A'])

 
            elif 'ins' in row[1]:
                end = 'ins'
                pos = ((cDNA.split(start))[1].split(middle)[0])
                second_pos = ((cDNA.split(middle))[1].split(end)[0])
                sequence = cDNA.split('ins',1)[1]

                if '-' in pos:
                    first_pos = pos.split('-')[0]
                    intron_pos = pos.split('-')[-1]
                    if '-' in second_pos:
                        second_pos_split = second_pos.split('-')[0]
                        second_pos_intron = second_pos.split('-')[-1]
                        if int(intron_pos) >=51:
                            pass
                        elif int(second_pos_split) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+['-'+intron_pos]+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                    elif '+' in second_pos:
                        second_pos_split = second_pos.split('+')[0]
                        second_pos_intron = second_pos.split('+')[-1]
                        if int(intron_pos) >=51:
                            pass
                        elif int(second_pos_split) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+['-'+intron_pos]+[second_pos_split]+[second_pos_intron]+[sequence])

                    else:
                        # second position is exonic
                        if int(intron_pos) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+['-'+intron_pos]+[second_pos]+['0']+[sequence])

                elif '+' in pos:
                    first_pos = pos.split('+')[0]
                    intron_pos = pos.split('+')[-1]

                    if int(first_pos) > 200000:
                        pass
  
                    elif '-' in second_pos:
                        second_pos_split = second_pos.split('-')[0]
                        second_pos_intron = second_pos.split('-')[-1]
                        if int(intron_pos) >=51:
                            pass
                        elif int(second_pos) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+[intron_pos]+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                    elif '+' in second_pos:
                        second_pos_split = second_pos.split('+')[0]
                        second_pos_intron = second_pos.split('+')[-1]
                        if int(intron_pos) >=51:
                            pass
                        elif int(second_pos_split) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+[intron_pos]+[second_pos_split]+[second_pos_intron]+[sequence])
                    else:
                        # second position is exonic
                        if int(intron_pos) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[first_pos]+[intron_pos]+[second_pos]+['0']+[sequence])
                else:
                    # the first postion is exonic

                    if int(pos) >200000:
                        pass

                    elif '-' in second_pos:
                        second_pos_split = second_pos.split('-')[0]
                        second_pos_intron = second_pos.split('-')[-1]
                        if int(intron_pos) >=51:
                            pass
                        elif int(second_pos_split) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[pos]+['0']+[second_pos_split]+['-'+second_pos_intron]+[sequence])

                    elif '+' in second_pos:
                        second_pos_split = second_pos.split('+')[0]
                        second_pos_intron = second_pos.split('+')[-1]
                        if int(intron_pos) >=51:                           
                            pass
                        elif int(second_pos_split) >200000:
                            pass
                        elif int(second_pos_intron) >=51:
                            pass
                        else:
                            writer.writerow(row+['insertion']+[pos]+['0']+[second_pos_split]+[second_pos_intron]+[sequence])

                    else:
                        # second position is exonic
                        if int(second_pos) >200000:
                            pass
                        writer.writerow(row+['insertion']+[pos]+['0']+[second_pos]+['0']+[sequence])


        elif 'del' in row[1]:
            end = 'del'
            pos = ((cDNA.split(start))[1].split(end)[0])
            if '-' in pos:
                first_pos = pos.split('-')[0]
                intron_pos = pos.split('-')[-1]
                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN deletion']+[first_pos]+['-'+intron_pos]+['N/A']+['0']+['N/A'])
            elif '+' in pos:
                first_pos = pos.split('+')[0]
                intron_pos = pos.split('+')[-1]
                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN deletion']+[first_pos]+[intron_pos]+['N/A']+['0']+['N/A'])
            else:
                if int(pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN deletion']+[pos]+['0']+['N/A']+['0']+['N/A'])

        elif 'dup' in row[1]:
            end = 'dup'
            pos = ((cDNA.split(start))[1].split(end)[0])
            if '-' in pos:
                first_pos = pos.split('-')[0]
                intron_pos = pos.split('-')[-1]
                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN duplication']+[first_pos]+['-'+intron_pos]+['N/A']+['0']+['N/A'])
            elif '+' in pos:
                first_pos = pos.split('+')[0]
                intron_pos = pos.split('+')[-1]
                if int(intron_pos) >=51:
                    pass
                elif int(first_pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN duplication']+[first_pos]+[intron_pos]+['N/A']+['0']+['N/A'])
            else:
                if int(pos) >=200000:
                    pass
                else:
                    writer.writerow(row+['SN duplication']+[pos]+['0']+['N/A']+['0']+['N/A'])


        elif '' in row[1]:
            pass

        else:
            pass
