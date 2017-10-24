def Filter(inputevent,filename):
    novel={}
    with open(inputevent,"r") as inputfile:
        for line in inputfile:
            if '#' in line:
                continue
            line=line.strip().split('\t')
            if len(line)<3:
                continue
            if '.' in line[2][1] and '.' in line[2][-2]:
                pepseq=line[2].split('.')[1].translate(None, '1234567890._:+-*!@#$?')
            else:
                pepseq=line[2].translate(None, '1234567890._:+-*!@#$?')
                
            novel[line[0]]={'type':line[1], 'seq':pepseq,'n':True,'id':'','seqold':line[2]}
    

    with open(filename,"r") as sf:
        fi=0
        for line in sf:
            if fi==0:
                fi=1
                continue
            line=line.strip().split("\t")
            if len(line)<2:
                continue
            for event in novel:
                if novel[event]['n']==True:
                    if novel[event]['seq'] in line[1]:
                        novel[event]['n']=False
                        novel[event]['id']=line[0]
    output=[]
    output_ref=[]

    output.append("#Event\tEvent_Type\tPeptide\n") 
    output_ref.append("#UniProtID\tEvent\tEvent_Type\tPeptide\n") 
             
    for event in novel:
        if novel[event]['n']==True:
            output.append("%s\t%s\t%s\n" %(event,novel[event]['type'],novel[event]['seqold']))
        else:
            output_ref.append("%s\t%s\t%s\t%s\n" %(novel[event]['id'],event,novel[event]['type'],novel[event]['seqold']))
    return (output, output_ref)
    
    

    
    
