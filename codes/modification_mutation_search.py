#USAGE: python modification_mutation_search.py outputcptac.txt alternate.fasta80.tab new_mod_alter_output80.txt
#PATH:/Users/anb013/Source_Code_GIT/Proteogenomics_Validation/Data
#most_common_mod={15.994915:'M',14.015650:'K',-17.026549:'Q',79.966331:'S,T,Y',43.005814:'*',42.010565:'*',0.984016:'N,Q'}
# Oxidation	+15.994915	M	OPTIONAL
# Lysine Methylation	+14.015650	K	OPTIONAL
# Pyroglutamate Formation	-17.026549	Q	OPTIONAL, N-TERMINAL
#Phosphorylation	+79.966331	STY	OPTIONAL
#N-terminal Carbamylation	+43.005814	*	OPTIONAL, N-TERMINAL
#N-terminal Acetylation	+42.010565	*	OPTIONAL, N-TERMINAL
#Deamidation	+0.984016	NQ	OPTIONAL
#(theoretical mass) / ((1/errorPPM)*1,000,000) = error(Da)event

def rg(string,kwd,kwd_count,flank,pos):
    aa_replace={'N':'GG', 'Q':'K,GA', 'R':'GV', 'C':'SS','W':'GE,AD,SV','Y':'M','P':'D,N','E':'Q,D,T','D':'N,S','I':'L,P','K':'Q','L':'I,P','F':'M','H':'D'}
    double_aa_replace={'GG':'N', 'GA':'Q', 'GV':'R','GE':'W', 'SS':'C','AD':'W','SV':'W'}
    start=pos
    end=flank-pos
    if kwd_count.has_key(string):
        fl=kwd_count[string]
        fl=fl+1
        kwd[string][fl]=(start,end)
        kwd_count[string]=fl
        return (kwd,kwd_count)
    else:
        kwd[string]={1:(start,end)}
        kwd_count[string]=1
    prev=''
    for i,w in enumerate(string):
        dk=prev+w
        if aa_replace.has_key(w):
            aalist=aa_replace[w].split(',')
            for aa in aalist:
                newstring=string[:i]+aa+string[i+1:]
                if kwd_count.has_key(newstring):
                    fl=kwd_count[newstring]
                    fl=fl+1
                    kwd[newstring][fl]=(start,end)
                    kwd_count[newstring]=fl
                else:
                    kwd[newstring]={1:(start,end)}
                    kwd_count[newstring]=1
        if double_aa_replace.has_key(dk):
            aa=double_aa_replace[dk]
            if i-1>0:
                newstring=string[:i-1]+aa+string[i+1:]
            else:
                newstring=aa+string[i+1:]
            if kwd_count.has_key(newstring):
                fl=kwd_count[newstring]
                fl=fl+1
                kwd[newstring][fl]=(start,end)
                kwd_count[newstring]=fl
            else:
                kwd[newstring]={1:(start,end)}
                kwd_count[newstring]=1
        prev=w
    return (kwd,kwd_count)

def create_pep_table(inputf,outf_fasta): #,delta):
    complete_KW={}
    KW_count={}
    for line in inputf:
        if '#' in line:
            continue
        line=line.strip().split('\t')
        kword=line[5]
        pos=int(line[6])
        flank=int(len(line[2]))
        (complete_KW,KW_count)=rg(kword,complete_KW,KW_count,flank,pos)
        outf_fasta.write(">%s|Novel\n%s\n" %(line[0],line[2]))
    #print KW_count
    return complete_KW #temp_dic_str

def print_kw(sequence,kw,string_index,ckw,outf,outf_fasta,ID,length,all_pep):#,all_pep_s):
    i=string_index
    if ckw.has_key(kw):
        for flank in ckw[kw]:
            (start,end)=ckw[kw][flank]
            start=i-start-5
            if start<0:
                start=0
            end=i+end+5
            if end>length:
                end=length
            s=sequence[start:end]
            if not all_pep.has_key(s):
                outf.write("%s|%s\t%s\n" %(kw,ID,s))
                outf_fasta.write(">%s|Modified\n%s\n" %(ID,s))
                all_pep[s]=1
                #all_pep_s[s[2:-2]]=1

def get_match_seq(sequence,ckw,outf,outf_fasta,ID,all_pep):#,all_pep_s):
    length=len(sequence)
    for i in range(0,length):
        kw=sequence[i:i+3]
        kw4=sequence[i:i+4]
        kw5=sequence[i:i+5]
        if len(kw)<3:
            break
        if len(kw)==3:
            print_kw(sequence,kw,i,ckw,outf,outf_fasta,ID,length,all_pep)#,all_pep_s)
        if len(kw4)==4:
            print_kw(sequence,kw4,i,ckw,outf,outf_fasta,ID,length,all_pep)#,all_pep_s)
        if len(kw5)==5:
            print_kw(sequence,kw5,i,ckw,outf,outf_fasta,ID,length,all_pep)#,all_pep_s)

def generate_alternatives(inputf,seqfile,singlefasta,alt_file,cpath):
    FASTA=seqfile
    outf=open(alt_file,"w")
    outf.write("#Keyword\tAlternative_Reference_Peptide\n")
    outf_fasta=open(singlefasta,"w")
    ckw=create_pep_table(inputf,outf_fasta)
    all_pep={}
    #all_pep_s={}        

    with open(FASTA,"r") as fastaf:
        fi=0
        for seq_record in fastaf:
            if fi==0:
                fi=1
                continue
            seq_record=seq_record.strip().split('\t')
            Header = seq_record[0]
            #UNIPROT_ID = Header
            Protein_sequence = seq_record[1]
            get_match_seq(Protein_sequence,ckw,outf,outf_fasta,Header,all_pep)#,all_pep_s)

    outf_fasta.close()
    outf.close()
