import os
import os.path
from pyteomics import mzml
from pyteomics import mzxml
                
def read_mzml(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile):
    if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
        with mzml.read(PATH) as reader:
            for scanindex,spectrum in enumerate(reader):
                if scanlist.has_key(scanindex):
                    try:
                        pev=event_scan[fname][scanindex]['pev']
                        output.append("%s\t%d\t%s\t%s\t%s\n" %(pev,newscanno,spec_outfile,event_scan[fname][scanindex]['pseq'],event_scan[fname][scanindex]['etype']))
                        newscanno+=1
                        charge=int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        mz=spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        soutput.write("BEGIN IONS\nTITLE=controllerType=0 controllerNumber=1 scan=%d\nCHARGE=%d+\nPEPMASS=%s\n" %(newscanno,charge,mz))
                        for x,y in zip(spectrum['m/z array'],spectrum['intensity array']):
                            soutput.write("%s %s\n" %(x,y))
                        soutput.write("END IONS\n\n")
                    except:
                        print("Error reading mzML file")
    return newscanno

def read_mzxml(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile):
    if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
        with mzxml.read(PATH) as reader:
            for scanindex,spectrum in enumerate(reader):
                if scanlist.has_key(scanindex):
                    try:
                        pev=event_scan[fname][scanindex]['pev']
                        output.append("%s\t%d\t%s\t%s\t%s\n" %(pev,newscanno,spec_outfile,event_scan[fname][scanindex]['pseq'],event_scan[fname][scanindex]['etype']))
                        newscanno+=1
                        charge=int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        mz=spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        soutput.write("BEGIN IONS\nTITLE=controllerType=0 controllerNumber=1 scan=%d\nCHARGE=%d+\nPEPMASS=%s\n" %(newscanno,charge,mz))
                        for x,y in zip(spectrum['m/z array'],spectrum['intensity array']):
                            soutput.write("%s %s\n" %(x,y))
                        soutput.write("END IONS\n\n")
                    except:
                        print("Error reading mzML file")
                    
    return newscanno


def read_mgf(PATH,scanlist,event_scan,fname,output,soutput,newscanno,spec_outfile):
    flag=0
    if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
        with open(PATH,"r") as mgf:
            scanindex=0
            for line in mgf:
                if len(line)>0:
                    if 'BEGIN IONS' in line:
                        if scanlist.has_key(scanindex):
                            pev=event_scan[fname][scanindex]['pev']
                            output.append("%s\t%d\t%s\t%s\t%s\n" %(pev,newscanno,spec_outfile,event_scan[fname][scanindex]['pseq'],event_scan[fname][scanindex]['etype']))
                            newscanno+=1
                            soutput.write(line)
                            flag=1
                        scanindex+=1
                    elif flag==1 and 'END IONS' not in line:
                        soutput.write(line)
                    elif flag==1 and 'END IONS' in line:
                        soutput.write(line)
                        flag=0
                    
    return newscanno
