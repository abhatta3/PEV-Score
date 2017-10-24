import os, sys, getopt
import time
from utility import getannotdictionary
from novel_filtering import Filter
from Make_list import make_list
from extract_spec import extracter
from get_keywords import get_key
from modification_mutation_search import generate_alternatives
from read_mzid import mzid_tsv_convert
from Rescore import rescore

def print_help():
    print 'python pev.py -i <input event file> -s <input spectra files>  -d <input database search files> -p <parameter file> -o <outputfile>'

def get_prm(line,param_d):
    line=line.split(":")[1]
    pk_l=line.strip().split('=')
    kv=pk_l[0]
    kd=int(pk_l[1])
    if param_d.has_key(kv):
        param_d[kv]=kd
    
    return param_d

def write_mod_file(line,fl,fixed_mod):
    fl.write(line.split(":")[1])
    if 'c:' in line:
        line=line.split(":")[1].split(',')
        fixed_mod[line[1]]=line[0]
    
    return fixed_mod
    
def main(argv):
    input_event = ''
    input_spectra = ''
    input_psm = ''
    paramfile=''
    fixed_mod={}
    param_d={'s_type':1,'d_pep':8,'d_index':1,'d_fname':0,'i_ev':0,'i_seq':2,'i_type':1}
    
    try:
        opts, args = getopt.getopt(argv,"hi:s:d:p:o:")
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in "-i":
            input_event = arg
        elif opt in "-p":
            paramfile = arg
        elif opt in "-o":
            outfile = arg
        elif opt in "-d":
            input_psm = arg
        elif opt in "-s":
            input_spectra = arg
    if input_event == '' or outfile=='' or input_psm=='' or input_spectra=='' or paramfile=='':
        print("Error: Wrong parameter. File name missing") 
        print_help()
        sys.exit()
    cpath=os.path.dirname(os.path.realpath(sys.argv[0]))
    reffilename=cpath+"/data/uniprot-all.tab"
    rmssfilename=cpath+"/data/AA_mass"
    par_map_mgf={}
    #for cluster run when a param.xml describe file mapping data
    if os.path.exists('params/params.xml'):
        with open('params/params.xml',"r") as param: 
            for line in param:
                if '"upload_file_mapping">mgffile' in line:
                    line=line.strip().split(">")[1]
                    new_mgf=line.split("|")[0]
                    old_mgf=line.split("/")[-2].split("<")[0]
                    par_map_mgf[old_mgf]=new_mgf
    ms_dict=getannotdictionary(rmssfilename)
    #read parameter from file
    #d_pep=8,d_index=1,d_fname=0,i_ev=0,i_seq=2,i_type=1
    #,itraq_tag
    modfile=open(outfile+"/Modifications_msgf.txt","w")
    with open(paramfile,"r") as pmdata:
        for line in pmdata:
            if "#"==line[0]:
                continue
            elif "p:" in line:
                param_d=get_prm(line,param_d)
            elif 'c:' in line or 'v:' in line or 'm:' in line:
                    fixed_mod=write_mod_file(line,modfile,fixed_mod)
            elif  'n:' in line:
                msgfpar=line.strip().split(':')[1]
                

    modfile.close()
    print fixed_mod
    (novel_list,ref_list)=Filter(input_event,reffilename)
    print("Novel filtering done.")
    make_list_psm_list=make_list(par_map_mgf,novel_list, input_psm, param_d['d_pep'],param_d['d_index'],param_d['d_fname'],param_d['i_ev'],param_d['i_seq'],param_d['i_type'])
    
    print("PSM list done")    
    singlemgf=outfile+"/single.mgf"
    
    extracted_psm_list=extracter(make_list_psm_list,input_spectra,singlemgf,param_d['s_type'])

    print("single mgf done\n")

    start_time = time.time()
    
    key_out=get_key(extracted_psm_list,singlemgf,ms_dict,fixed_mod,4)

    print("Keyword Executes in (seconds) ---\t%s\n" %(time.time() - start_time))
    
    
    singlefasta=outfile+"/single.fa"
    
    alt_outf=outfile+"/Alternatives.txt"

    start_time = time.time()
    
    generate_alternatives(key_out,reffilename,singlefasta,alt_outf,cpath)

    print("Alternative Executes in (seconds) ---\t%s\n" %(time.time() - start_time))


    mzidout=outfile+"/msgfout.mzid"

    start_time = time.time()
    
    os.system("java -Xmx3500M -jar "+cpath+"/msgfplus/MSGFPlus.jar  -d  "+singlefasta+" -s  "+singlemgf+"  -o "+mzidout+" "+msgfpar+" -mod  "+outfile+"/Modifications_msgf.txt")

    print("MSGF+ search Executes in (seconds) ---\t%s\n" %(time.time() - start_time))

    
    #print("MSGF search done")
    tsvout=outfile+"/msgfout.tsv"
    mzid_tsv_convert(mzidout,tsvout)
    print("mzid to tsv done")
    outputfile=outfile+"/result.txt"

    start_time = time.time()
    
    rescore(tsvout,outputfile,extracted_psm_list,singlemgf,novel_list,key_out,fixed_mod,ref_list) #,pep_by_kw)

    print("Rescore Executes in (seconds) ---\t%s\n" %(time.time() - start_time))

            
if __name__ == "__main__":
    start_time = time.time()
    main(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
