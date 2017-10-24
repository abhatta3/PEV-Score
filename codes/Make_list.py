#https://pythonhosted.org/pyteomics/data.html
import os

def make_list(mgf_name,inputevent,scanf,d_pep=8,d_index=1,d_fname=0,i_ev=0,i_seq=2,i_type=1):
    output=[]
    pep={}
    event_type={}
    novel={}
    for line in inputevent:
        if '#' in line:
            continue
        line=line.strip().split('\t')
        if len(line)<3:
            continue
        if '.' in line[i_seq][1] and '.' in line[i_seq][-2]:
            pepseq=line[i_seq][1:-1].translate(None, '+1234567890._:-*!@#$?')
        else:
            pepseq=line[i_seq].translate(None, '+1234567890._:-*!@#$?')
         
        novel[pepseq]=line[i_seq]
        if pep.has_key(pepseq):
            old=pep[pepseq]
            old.append(line[i_ev])
            pep[pepseq]=old
        else:
            pep[pepseq]=[line[i_ev]]
        event_type[pepseq]=line[i_type]
    output.append("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n") 

    #NPATH=os.path.dirname(os.path.abspath(scanf))
    if os.path.isdir(scanf):
        for scanf_filename in os.listdir(scanf):
            scanf_filename=scanf+'/'+scanf_filename
            with open(scanf_filename,"r") as sf:
                for line in sf:
                    if '#' in line:
                        continue
                    line=line.strip().split('\t')
                    if len(line)<8:
                        continue
                    if '.' in line[d_pep][1] and '.' in line[d_pep][-2]:
                        pepseq=line[d_pep][1:-1].translate(None, '+1234567890._:-*!@#$?')
                    else:
                        pepseq=line[d_pep].translate(None, '+1234567890._:-*!@#$?')
                    fname=line[d_fname]
                    if '=' in line[d_index]:
                        index=line[d_index].split('=')[1]
                    else:
                        index=line[d_index]
    
                    if pep.has_key(pepseq):
                        pev=pep[pepseq]
                        pis=pev[0]
                        if '/' in fname:
                            fname=fname.split("/")[-1]
                        if mgf_name.has_key(fname):
                            fname=mgf_name[fname] 
                        output.append("%s\t%s\t%s\t%s\t%s\n" %(pis,index,fname,line[d_pep],event_type[pepseq]))
            
    else:
        with open(scanf,"r") as sf:
            for line in sf:
                if '#' in line:
                    continue
                line=line.strip().split('\t')
                if len(line)<8:
                    continue
                if '.' in line[d_pep][1] and '.' in line[d_pep][-2]:
                    pepseq=line[d_pep][1:-1].translate(None, '+1234567890._:-*!@#$?')
                else:
                    pepseq=line[d_pep].translate(None, '+1234567890._:-*!@#$?')
                fname=line[d_fname]
                if '=' in line[d_index]:
                    index=line[d_index].split('=')[1]
                else:
                    index=line[d_index]

                if pep.has_key(pepseq):
                    pev=pep[pepseq]
                    pis=pev[0]
                    if '/' in fname:
                        fname=fname.split("/")[-1]
                    if mgf_name.has_key(fname):
                        fname=mgf_name[fname] 
                    output.append("%s\t%s\t%s\t%s\t%s\n" %(pis,index,fname,line[d_pep],event_type[pepseq]))
    return output


