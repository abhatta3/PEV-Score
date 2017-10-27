#CPTAC_Colon_Proteogenomic_Events_membrane = '12_14_2015_CPTAC_Colon_Proteogenomic_Events_membrane.txt'
#CSSHTMLBuild(CPTAC_Colon_Proteogenomic_Events_membrane)
import heapq

def getb(name,val,charge,ms_dict,fixed_mod,s):
    tlist=[]
    for i in range(0,s):
            l=name[i]
            if l in ms_dict.keys():
                    oldval=val
                    val=float(ms_dict[l])+oldval
                    if fixed_mod.has_key(l):
                        val=val+float(fixed_mod[l])
                    m=(val+charge)/charge
                    tlist.append(m)
    return tlist

#place 1 to output vector if two picks differ by 1 aminoacid
def get_position_mass(pep,scan,pm,ms_dict,fixed_mod):
    pep=pep.translate(None, '1234567890:._+-*!@#$?')
    rpep=pep[::-1]
    s=len(pep)
    if fixed_mod.has_key('*'):
        additional_mass_at_n_terminal=float(fixed_mod['*'])
    else:
        additional_mass_at_n_terminal=0.0

    position_mass_b={}
    position_mass_y={}
    position_mass_t={}


    for i in range(0,s):
        position_mass_t[i]=0
        
    for charge in range(1,3):
        blist=getb(pep,additional_mass_at_n_terminal,charge,ms_dict,fixed_mod,s)
        ylist=getb(rpep,18.0,charge,ms_dict,fixed_mod,s)

        ylist=ylist[::-1]
        for i in range(0,s):
            position_mass_b[i]=0
            position_mass_y[i]=0
            if  i<len(blist) and i<len(ylist):
                if i==0:
                    bi=int(blist[i])
                    if scan.has_key(bi):
                        bival=scan[bi]['x']
                        inv=abs(float(bival)-blist[i])
                        if inv<pm:
                            position_mass_b[i]=1
                    elif scan.has_key(bi+1):
                        bival=scan[bi+1]['x']
                        inv=abs(float(bival)-blist[i])
                        if inv<pm:
                            position_mass_b[i]=1
                elif i==(len(ylist)-1):
                    yi=int(ylist[i-1])
                    if scan.has_key(yi):
                        yival=scan[yi]['x']
                        inv=abs(float(yival)-ylist[i-1])
                        if inv<pm:
                            position_mass_y[i]=1
                    elif scan.has_key(yi+1):
                        yival=scan[yi+1]['x']
                        inv=abs(float(yival)-ylist[i-1])
                        if inv<pm:
                            position_mass_y[i]=1                    
                else:
                    bi=int(blist[i])
                    if scan.has_key(bi):
                        bival=scan[bi]['x']
                        inv=abs(float(bival)-blist[i])
                        if inv<pm: # and position_mass_b[i-1]==1:
                            position_mass_b[i]=1
                    elif scan.has_key(bi+1):
                        bival=scan[bi+1]['x']
                        inv=abs(float(bival)-blist[i])
                        if inv<pm:# and position_mass_b[i-1]==1:
                            position_mass_b[i]=1
                    yi=int(ylist[i-1])
                    if scan.has_key(yi):
                        yival=scan[yi]['x']
                        inv=abs(float(yival)-ylist[i-1])
                        if inv<pm: # and position_mass_y[i-1]==1:
                            position_mass_y[i]=1
                    elif scan.has_key(yi+1):
                        yival=scan[yi+1]['x']
                        inv=abs(float(yival)-ylist[i-1])
                        if inv<pm: # and position_mass_y[i-1]==1:
                            position_mass_y[i]=1

        #update tags         
        for i in range(0,s-1):
            val=position_mass_b[i]
            if val > 0:
                inext=i+1
                if position_mass_b.has_key(inext):
                    nval=position_mass_b[inext]
                    if nval > 0:
                        position_mass_t[inext]=1 #only b

        for i in range(0,s-1):
            val=position_mass_y[i]
            if val > 0:
                inext=i+1
                if position_mass_y.has_key(inext):
                    nval=position_mass_y[inext]
                    if nval > 0:
                        if position_mass_t[i]==0:
                            position_mass_t[i]=-1 #only y
                        elif  position_mass_t[i]==1:
                            position_mass_t[i]=2 #both b and y
    #print pep, position_mass_t
    return  position_mass_t


    
def make_psm_pattern(mapscan,scanfile,ms_dict,fixed_mod):
    delta=0.5
    scan_i=-1
    eventsupport={}
    with open(scanfile,"r") as sf:
        matched=False
        for line in sf:
           if 'BEGIN IONS' in line:
               scan_i+=1
               scand={}
               slist=[]
               stemp={}
               scount=0
               if mapscan.has_key(scan_i):
                   matched=True
                   event=mapscan[scan_i]['ev']
                   pep=mapscan[scan_i]['ns']
                   pepref=mapscan[scan_i]['rs']
               else:
                   matched=False
           elif 'END IONS' not in line and matched==True and '=' not in line:
               line=line.strip().split()
               if len(line)>=2:
                   xval=float(line[0])
                   yval=float(line[1])
                   stemp[scount]={'y':yval,'x':xval}
                   scount+=1
                   slist.append(yval)
           elif 'END IONS' in line and matched==True:
               l50=heapq.nlargest(50, slist)[-1]
               for i in stemp:
                   yval=stemp[i]['y']
                   xval=stemp[i]['x']
                   if yval>=l50:
                       li=int(xval)
                       if scand.has_key(li):
                           oyval=scand[li]['y']
                           if oyval<yval:
                               scand[li]={'y':yval,'x':xval}
                       else:
                           scand[li]={'y':yval,'x':xval}
               pm=delta
               pd=get_position_mass(pep,scand,pm,ms_dict,fixed_mod)
               pdref=get_position_mass(pepref,scand,pm,ms_dict,fixed_mod)
               eventsupport[event]={'np':pd,'rp':pdref}
    return eventsupport



def dt_format(supd,dt):
    nstr=dt[:2]
    i=-1
    for d in dt[2:-2]:
        if d >='A' and d<='Z':
            i+=1
            if supd.has_key(i):
                sbit=supd[i]
                if sbit==-1:
                    nstr+='<font color="red" size="5">'+d+'</font>' #y ion
                elif sbit==1:
                    nstr+='<font color="blue" size="5">'+d+'</font>' #b ion
                elif sbit==2:
                    nstr+='<font color="green" size="5">'+d+'</font>' #both b abd y ion
                else:
                    nstr+=d
            else:
                nstr+=d
        else:
            nstr+=d
    nstr+=dt[-2:]
    return nstr



def CSSHTMLBuildFilter(inputfile,singlemgf,pevscan,ms_dict,fixed_mod):
    css_table_header = '''<!DOCTYPE html>\n\
    <html>\n\
    <head>\n\


    <script src="../sorttable.js"></script>\n\

    <style type="text/css">\n\

/* Sortable tables */
table.sortable thead {\n\
    background-color:#eee;\n\
    color:#666666;\n\
    font-weight: bold;\n\
    cursor: default;\n\
}\n\
 table.sortable th:not(.sorttable_sorted):not(.sorttable_sorted_reverse):not(.sorttable_nosort):after {
    content: '\\25B4\\25BE'
}\n\

#customers {\n\
    font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;\n\
    border-collapse: collapse;\n\
    width: 100%;\n\
}\n\

#customers td, #customers th {\n\
    border: 1px solid #ddd;\n\
    padding: 8px;\n\
}\n\

#customers tr:nth-child(even){background-color: #f2f2f2;}\n\
#customers tr:hover {background-color: #ddd;}\n\
#customers th {\n\
    padding-top: 12px;\n\
    padding-bottom: 12px;\n\
    text-align: left;\n\
    background-color: #4CAF50;\n\
    color: white;\n\
}\n\
\n\
    </style>\n\

    </head>\n\
    <body>\n\


  <h1  align="center">Proteogenomics Event Validation (PEV) Results</h1>\n\

  <h2  align="left">Character strings are colored to represent b/y ion fragment matches. [Blue: b ions] [Red: y ions] [Green: both b and y ions]</h2>\n\

<table class="sortable" id="customers">\n\
    	
    <thead>\n\
      <tr>\n\
        <th><a href=" " title="Event identifier: Same ID as in -i input event file." >EventID</a></th>\n\
        <th><a href=" " title="Novel peptide sequence from input:  Peptide sequence from -i input event file.">Novel_Peptide</a></th>\n\
        <th><a href=" " title="Rescored Novel Peptide: Novel_Peptide sequence from MSGF+ search. Should be the same as Novel_Peptide or NA (indicates not found from MSGF+ search).">Rescored Novel Peptide</a></th>\n\
        <th><a href=" " title="Alternative Reference Sequence: Reference sequence reported from MSGF+ search.">Alternative_Reference_Sequence</a></th>\n\
        <th><a href=" " title="UniProt ID for Alternative Reference: Protein ID (UniProt) for reported reference sequence.">UniProtID</a></th>\n\
        <th><a href=" " title="Variant Type: Description of variation type from -i input event file. ">Variant_Type</a></th>\n\
        <th><a href=" " title="PEV Score: Computed PEV score.">PEV_Score</a></th>\n\
        <th><a href=" " title="PEV Category: PEV+ when PEV_Score>0. This indicates proteogenomics (novel) peptide is better than reference peptides. PEV- when PEV_Score<0. This indicates reference peptide is better than novel peptide. PEVzero when PEV_Score==0. This indicates novel peptide and reference peptides are equally posible match. PEVexact-reference indicates that input novel peptide is matched exactly with a reference peptide sequence. PEVunscored indicates all the MSGF+ PSMs are reported with a spectrum level P-value greater than 0.05 (insignificant). 1 is good and represent unambiguous mapping.">PEV_Category</a></th>\n\
        <th><a href=" " title="All Higher Scored Alternative Reference: When one or more reference peptide is better than novel peptide this shows the better scored reference peptides." >Alternatives</a></th>\n\
        <th><a href=" " title="All Higher Scored Alternative Reference ProteinID: Protein ID (UniProt) for all the better scored reference peptides.">AlternativesID</a></th>\n\
      </tr>\n\
     </thead>\n\
     <tbody>\n'''
     
    #<th>Reference_Peptide</th>\n\
    htmlfile=inputfile+'.html' 
    pv_seq_bit=make_psm_pattern(pevscan,singlemgf,ms_dict,fixed_mod)
            
        
    with open(inputfile,'r') as data, open(htmlfile,'w') as html_file:
        data_lines = data.readlines()[1:]
        html_file.write(css_table_header)
        for item in data_lines:
            item_list = item.split('\t')
            html_file.write('  <tr>\n')
            for i, item in enumerate(item_list):
                if i==0:
                    eid=item
                dt = item.strip('\n')
                if dt=='':
                    dt='-'
                if i==2 or i==3:
                    if dt!='NA':
                        if pv_seq_bit.has_key(eid):
                            if i==2:
                                supd=pv_seq_bit[eid]['np']
                            elif i==3:
                                supd=pv_seq_bit[eid]['rp']

                            dt=dt_format(supd,dt)
                    
                html_file.write('    <td align="center">%s</td>\n'%(dt))
                        
            html_file.write('  </tr>\n')
        html_file.write('</tbody>\n</table>\n\n</div>\n</body>\n</html>\n')



#CSSHTMLBuildFilter("../results/result.txt")




