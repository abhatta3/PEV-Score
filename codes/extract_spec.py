#https://pythonhosted.org/pyteomics/data.html
import os
import os.path

def extracter(inputf,spec_dir,spec_outfile,s_type):
    soutput=open(spec_outfile,"w")
    output=[]
    event_scan={}
    for line in inputf:
        if '#' in line:
            continue
        line=line.strip().split('\t')
        if len(line)>=5:
            index=int(line[1])
            fname=line[2]
            pev=line[0]
            peptide=line[3]
            etype=line[4]
            if event_scan.has_key(fname):
                sdic=event_scan[fname]
                sdic[index]={'pev':pev,'pseq':peptide,'etype':etype}
                event_scan[fname]=sdic
            else:
                temp={'pev':pev,'pseq':peptide,'etype':etype}
                sdic={index:temp}
                event_scan[fname]=sdic
    newscanno=0
    output.append("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n")            
    for fname in event_scan:
        scanlist=event_scan[fname]
        #PATH='mgffile/'+fname
        if os.path.isdir(spec_dir):
            PATH=spec_dir+'/'+fname
            print "AS DIRECTORY"
        else:
            #PATH=os.path.dirname(os.path.abspath(spec_dir))
            PATH=fname
            print "AS PATH"
        if s_type==2:
            from write_single_mgf import read_mzml
            newscanno=read_mzml(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile)
        elif s_type==3:
            from write_single_mgf import read_mzxml
            newscanno=read_mzxml(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile)
        else:
            from write_single_mgf import read_mgf            
            newscanno=read_mgf(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile)

    soutput.close()
    return output
                             
