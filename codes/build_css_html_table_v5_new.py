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
                <h1 class="text-center">Proteogenomics Event Validation (PEV) Results</h1>\n\
            </div>\n\
        </div>\n\
        <div class="form-group" style="border: 1px solid #ddd">\n'''

    table_column_name ='''<table id="table1"  class="sortable">\n\
    	
    <thead>\n\
      <tr>\n\
        <th><a href=" " title="Event identifier: Same ID as in -i input event file." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Event ID<br>[@]</a></th>\n\
        <th><a href=" " title="Novel peptide sequence from input:  Peptide sequence from -i input event file." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Novel Peptide</a></th>\n\
        <th><a href=" " title="Rescored Novel Peptide: Novel_Peptide sequence from MSGF+ search. Should be the same as Novel_Peptide or NA (indicates not found from MSGF+ search)." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Rescored Novel Peptide</a></th>\n\
        <th><a href=" " title="Alternative Reference Sequence: Reference sequence reported from MSGF+ search." style="background-color:#FFFFFF;color:#000000;text-decoration:none">Alternative Reference Sequence</a></th>\n\
        <th><a href=" " title="UniProt ID for Alternative Reference: Protein ID (UniProt) for reported reference sequence. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">UniProtID</a></th>\n\
        <th><a href=" " title="Variant Type: Description of variation type from -i input event file. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Variant Type</a></th>\n\

        <th><a href=" " title="PEV Score: Computed PEV score. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">PEV Score</a></th>\n\
        <th><a href=" " title="PEV Category: PEV+ when PEV_Score>0. This indicates proteogenomics (novel) peptide is better than reference peptides. PEV- when PEV_Score<0. This indicates reference peptide is better than novel peptide. PEVzero when PEV_Score==0. This indicates novel peptide and reference peptides are equally posible match. PEVexact-reference indicates that input novel peptide is matched exactly with a reference peptide sequence. PEVunscored indicates all the MSGF+ PSMs are reported with a spectrum level P-value greater than 0.05 (insignificant). 1 is good and represent unambiguous mapping.  " style="background-color:#FFFFFF;color:#000000;text-decoration:none">PEV Category</a></th>\n\
        <th><a href=" " title="All Higher Scored Alternative Reference: When one or more reference peptide is better than novel peptide this shows the better scored reference peptides. " style="background-color:#FFFFFF;color:#000000;text-decoration:none">Alternatives</a></th>\n\
        <th><a href=" " title="All Higher Scored Alternative Reference ProteinID: Protein ID (UniProt) for all the better scored reference peptides." style="background-color:#FFFFFF;color:#000000;text-decoration:none">AlternativesID</a></th>\n\


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
                #Event_ID
                dt = item.strip('\n')
                if dt=='':
                    dt='-'
                    
                html_file.write('    <td class="center">%s</td>\n'%(dt))
                        
            html_file.write('  </tr>\n')
        html_file.write('</tbody>\n</table>\n\n</div>\n</body>\n</html>\n')









