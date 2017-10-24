#CPTAC_Colon_Proteogenomic_Events_membrane = '12_14_2015_CPTAC_Colon_Proteogenomic_Events_membrane.txt'
#CSSHTMLBuild(CPTAC_Colon_Proteogenomic_Events_membrane)


def CSSHTMLBuildFilter(inputfile,spect_task_id=''):
    css_table_header = '''<!DOCTYPE html>\n\
    <html>\n\
    <head>\n\

    <link href="http://cdn.bootcss.com/bootstrap/3.3.0/css/bootstrap.min.css" type="text/css" rel="stylesheet">\n\
    <script type="text/javascript" src="https://dl.dropboxusercontent.com/u/63775276/sorttable.js"></script>\n\


    <style type="text/css">\n\
    table.sortable {\n\
        font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;\n\
        width: 100%;\n\
        border-collapse: collapse;\n\
    }\n\
    \n\
    table.sortable td, #mutated_peptides th {\n\
        font-size: 1em;\n\
        border: 2px solid #000000;\n\
        padding: 2px 3px 2px 3px;\n\
    }\n\
    \n\
    table.sortable tr:nth-child(odd){background-color: #f2f2f2}\n\
    \n\
    table.sortable th {\n\
        font-size: 1em;\n\
        text-align: center;\n\
        padding-top: 2px;\n\
        padding-bottom: 2px;\n\
        background-color: #4CAF50;\n\
        color: #ffffff;\n\
        border: 2px solid #000000;\n\
    }\n\
    \n\
    table.sortable td.center {\n\
        text-align: center;\n\
    }\n\
    \n\
    table.sortable td.canceratlas {\n\
        color: #8cd68f;\n\
        background-color: #8cd68f;\n\
        text-align: center;\n\
    }\n\annotfile
    \n\
    </style>\n\


    </head>\n\
    <body>\n\

       <div class="form-group">\n\
            <div class="col-sm-12">\n\
                <span class="help-block">Sort table with Left Click at Header. @ Click event to see spectra. You can select with <b>Ctrl</b>, <b>Shift</b>, <b>Common</b> Or <b>Ctrl + A</b>.</span>\n\
            </div>\n\
        </div>\n\


        <div class="form-group">\n\
            <div class="col-sm-12">\n\
                <h1 class="text-center">Proteogenomics Annotation Data Table Towards Target Selection</h1>\n\
            </div>\n\
        </div>\n\
        <div class="form-group" style="border: 1px solid #ddd">\n'''

    table_column_name ='''<table id="table1"  class="sortable">\n\
    	
    <thead>\n\
      <tr>\n\
        <th><a href=" " title="Numeric identifier for the novel event. Click in a ID to display the events on spectra viewer." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Event ID<br>[@]</a></th>\n\
        <th><a href=" " title="Novel peptide event label. For example, mutation, insertion, deletion, novel splice ... " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Event Type</a></th>\n\
        <th><a href=" " title="Sequence of novel peptide formated as Prefix-AA.Left-Peptide-segment:Right-Peptide-segment.Postfix-AA, where the mutation event could be either the entire Left or Right Peptide-segment or it could be at the right end of Left-Peptide-segment or left end of Right-Peptide-segment." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Peptide</a></th>\n\
        <th><a href=" " title="Reference AA sequence is the part of the reference sequence differed from novel peptide sequence." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Ref_aa</a></th>\n\


        <th><a href=" " title="Mutated AA sequence is the part of the novel sequence differed from reference peptide sequence. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Mut_aa</a></th>\n\
        <th><a href=" " title="Location of the first AA of Ref_aa sequence in a protein. " style="background-color:#FFFFFF;color:#000000;text-decoration:none"> Ref_Location</a></th>\n\

        <th><a href=" " title="Location of the first AA of Mut_aa sequence in novel peptide. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Mut_loc_at_peptide</a></th>\n\
        <th><a href=" " title="Number of mapped gene to the genomic location of the mutation event. 1 is good and represent unambiguous mapping.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Total Number of<br>Overlapping Gene</a></th>\n\

        <th><a href=" " title="Name of the Overlapping gene. Each row holds one gene name and corresponding values based on that gene. Separate rows are created for each gene when multiple overlapping genes. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Overlapping<br>Gene</a></th>\n\
        <th><a href=" " title="A protein for the mapped gene that have a matched reference sequence for the novel event peptide. The matched reference sequences are determined from a regular expression based sequence match between novel event peptide sequence and reference protein sequence. For example, a novel event peptide Prefix-AA.Left-Peptide-segment:Right-Peptide-segment2.Postfix-AA make 3 regular expression- 1) Prefix-AA.[1-Left-Peptide-segment-1]?[2-Right-Peptide-segment].Postfix-AA, 2) Prefix-AA.[1-Left-Peptide-segment-1]?, and 3) ?[2-Right-Peptide-segment].Postfix-AA." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Overlapping Protein<br>(UniProtKB_AC)</a></th>\n\


        <th><a href=" " title="The keywords are retrieved from Uniprot database entry for the protein.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">UniProtKB <br> Key Words</a></th>\n\
        <th><a href=" " title="Protein position is determined, from different domain positions data that are available at Uniprot database. To be a membrane bound protein the protein should have a TRANSMEM domain entry. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Protein_position</a></th>\n\

        <th><a href=" " title=" A link to the database" style="background-color:#FFFFFF;color:#000000;text-decoration:none">Protein Atlas<br>(Cancer Atlas)</a></th>\n\

        <th><a href=" " title="Chromosome for the genetic variant. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Chromosome</a></th>\n\

        <th><a href=" " title="Genomic location of the novel peptide. Different segments are separated by ';'. In case a point mutation, the location of the mutation is the end coordinate of the segment before letter 'M'.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Genomic_Location</a></th>\n\

        <th><a href=" " title="This number indicates how many novel peptides are identified from this gene. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Novel peptides<br>in overlapping gene</a></th>\n\
        <th><a href=" " title="This number indicates how many known peptides are identified from this gene. This indicates overall expression of the gene.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Known peptides<br>in overlapping gene</a></th>\n\

        <th><a href=" " title="Total number of spectra that are matched to this novel event peptide.     " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Total<br>Spec Count</a></th>\n\
        <th><a href=" " title="This number indicate how many genomic location for the mapped RNA-seq reads. Should be 1 to avoid any ambiguity." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Number of possible<br>genomic locations</a></th>\n\

        <th><a href=" " title="This indicate total charge of Ref_aa." style="background-color:#FFFFFF;color:#000000;text-decoration:none">WT_Charge</a></th>\n\
        <th><a href=" " title="This indicate total charge of Mut_aa.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">MT_Charge</a></th>\n\
        <th><a href=" " title="This is the charge difference " style="background-color:#FFFFFF;color:#000000;text-decoration:none">MT-WT_Charge</a></th>\n\
        <th><a href=" " title="This indicates total surface area of Ref_aa " style="background-color:#FFFFFF;color:#000000;text-decoration:none">WT_Surface</a></th>\n\
        <th><a href=" " title="This indicates total surface area of Mut_aa " style="background-color:#FFFFFF;color:#000000;text-decoration:none">MT_Surface</a></th>\n\
        <th><a href=" " title="Indicates difference of surface" style="background-color:#FFFFFF;color:#000000;text-decoration:none">MT-WT_Surface</a></th>\n\
        <th><a href=" " title="This value indicates the hydrophobicity of the mutation  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Kyte-Doolittle<br>Hydrophobic(>0)</a></th>\n\
        <th><a href=" " title="This value indicates the hydrophilic property of the mutation. Note: See the pdf figure for a detailed plot. The plot include all the AA of mutated peptide. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Hopp-Woods<br>Hydrophilic(>0)</a></th>\n\
        <th><a href=" " title="Protein coding mutations that are mapped at the mutation location of the protein are collected from the Uniprot. Click the link to open dbSNP page. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">dbSNP Overlap</a></th>\n\
        <th><a href=" " title="This measure the Reproducibility of novel event over spectra samples. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Event Recurrence<br> Over Number of<br>Spectra Samples</a></th>\n\
        <th><a href=" " title="This measure the Reproducibility of novel event over RNA-seq  samples.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">    Event Recurrence<br>Over Number of<br>RNA_Samples</a></th>\n\
        <th><a href=" " title="Indicates the maximum number of read that are aligned to mutation in any RNA-seq samples. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Maximum Read<br>Depth</a></th>\n\
        <th><a href=" " title='This is determined, first from the Uniprot data. The Ref_Location is compared against the domain annotations. If  Ref_Location is in a Extracellular domain of a membrane bound protein then the event is labelled as Targetable. For some of the protein Uniprot don't have annotation for Extracellular domain, may be because they are not annotated in the plasma membrane yet. For them if the protein have a TRANSMEM domain we used protter to map the protein and peptide against membrane boundary. Then labelled them as Targetable-Peptide if the novel event is at the external part of the membrane. For such entry we need to do literature search to see whether they are expressed in the plasma membrane or the membrane is internal to the cell.' style="background-color:#FFFFFF;color:#000000;text-decoration:none">Mutation Event<br>Position</a></th>\n\
        <th><a href=" " title='Average of log2 ratio of gene and EGFR expression over a set of samples. If multiple values, then the last number is for normal samples and all the other values are for tumor samples. Please see the pdf file for the plots over all the samples.
For example, in
3.47523245003 
0.767078751683
the value 3.47523245003 indicates that on average gene is 2^3.475 times more expressed than EGFR in tumor samples while 0.767078751683 indicates that for normal samples.'   style="background-color:#FFFFFF;color:#000000;text-decoration:none">log2(Gene/EGFR)<br>Protein Expression</a></th>\n\
        <th><a href=" " title='This indicates average rank of the gene expression over a set of samples. The rank average are computed from reciprocal of average of reciprocal rank to minimize the outlier effect. Smaller the rank number higher is the expression value. For multiple rows, the last row is for normal samples and others are for tumor samples. 
For example, in 

1231.66849879,3774.92781345
2566.27261344,3184.10412998

the values 1231.66849879,3774.92781345 are from tumor samples and 2566.27261344,3184.10412998 are from normal samples. The value 1231.66849879 gives the average rank of gene in tumor samples and 3774.92781345 is the average rank of EGFR in tumor samples. Similarly 2566.27261344 and 3184.10412998 are rank of gene and EGFR respectively for normal samples.'  style="background-color:#FFFFFF;color:#000000;text-decoration:none">Gene Rank, EGFR Rank <br> Protein Expression</a></th>\n\
        <th><a href=" " title="This field indicates average if log2 ratio between mutation and reference allele from RNA-seq tumor samples. The detailed plots are available in the pdf file. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">In-Silico PCR</a></th>\n\ 

      </tr>\n\
     </thead>\n\
     <tbody>\n'''
     
     

    #<th>Reference_Peptide</th>\n\
    htmlfile=inputfile+'.html' 
        
    with open(inputfile,'r') as data, open(htmlfile,'w') as html_file:
        
        data_lines = data.readlines()[1:]
        html_file.write(css_table_header)

        html_file.write(table_column_name)

            
        for item in data_lines:
            
            item_list = item.split('\t')
            
            html_file.write('  <tr>\n')

            for i, item in enumerate(item_list):
                #skip reference peptide. It is in the text file
                if i==3:
                    continue
                if i>38:
                    break;

                #Event_ID
                if i!=15 and i!=13 and i!=36:
                    dt = item.strip('\n')
                    if dt=='':
                        dt='-'
                        
                    if i==0:
                        if spect_task_id!='':
                            dt='<a href="http://proteomics2.ucsd.edu/ProteoSAFe/result.jsp?task='+spect_task_id+'&view=group_by_event#%7B%22%23Num_lowerinput%22%3A%22'+dt+'%22%2C%22%23Num_upperinput%22%3A%22'+dt+'%22%7D" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==9:
                        dt='<a href="'+item_list[15]+'" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==10:
                        dt='<a href="http://www.uniprot.org/uniprot/'+dt+'" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==12: 
                        dt='<a href="'+item_list[13]+'"  target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==14:
                        dt='<a href="'+dt+'" target="_blank"><FONT COLOR="#8cd68f">Atlas</FONT></a>'
                    if i==28:
                        dt=dt.split(';')[0]
                        hval=float(dt)
                        if hval>0:
                            hlabel="Hydrophobic;"
                            dt=hlabel+dt
                            
                    if i==35 or i==37:
                        combinedt=''
                        for part in dt.split(';'):
                            combinedt+=part+'<br>'
                        dt=combinedt
                            

                    if i==29:
                        dt=dt.split(';')[0]
                        hval=float(dt)
                        if hval>0:
                            hlabel="Hydrophilic;"
                            dt=hlabel+dt

                    if i==38:
                        dt=dt.split(';')[0]
                            
                    if i==30:
                        dt=dt.split('#')
                        dbl=''
                        for d in dt:
                            if dbl=='':
                                
                                if 'dbSNP:rs' in d:
                                    snpid=d.split(':')[1].translate(None, '.)rs')
                                    d='<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='+snpid+'">'+d+'</a>'
                                
                                dbl=d
                            else:
                                
                                if 'dbSNP:rs' in d:
                                    snpid=d.split(':')[1].translate(None, '.)rs')
                                    d='<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='+snpid+'">'+d+'</a>'
                                dbl=dbl+'<br><br>'+d
                        dt=dbl
                    if i == 4 or i == 5: 
		            html_file.write('    <td style="word-break:break-all;"  class="center">%s</td>\n'%(dt))
                    else:
		            html_file.write('    <td class="center">%s</td>\n'%(dt))
                        
            html_file.write('  </tr>\n')
        html_file.write('</tbody>\n</table>\n\n</div>\n</body>\n</html>\n')









