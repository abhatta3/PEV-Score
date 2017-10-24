def mzid_tsv_convert(infile,outfile):
    outf=open(outfile,"w")
    accn={}
    pep_details={}
    with open(infile,"r") as inf:
        s=0
        for line in inf:
            if '<DBSequence length=' in line:
                nline=line.strip().split()
                for word in nline:
                    if 'accession=' in word:
                        word=word.split('"')[1]
                        accname=word
                    elif 'id=' in word:
                        word=word.split('"')[1]
                        accid=word
                if not accn.has_key(accid):
                    accn[accid]=accname
            if '<PeptideEvidence ' in line:
                nline=line.strip().split()
                for word in nline:
                    if 'dBSequence_ref=' in word:
                        word=word.split('"')[1]
                        accnid=word
                    elif 'peptide_ref=' in word:
                        word=word.split('"')[1]
                        peptide=word
                    elif 'pre=' in word:
                        word=word.split('"')[1]
                        peptide_pre=word
                    elif 'post=' in word:
                        word=word.split('"')[1]
                        peptide_post=word
                if not pep_details.has_key(peptide):
                    if accn.has_key(accnid):
                        accnname=accn[accnid]
                    else:
                        accnname='-'
                    pep_details[peptide]={'name':accnname,'pre':peptide_pre,'post':peptide_post}
                        
                
                
            if "<SpectrumIdentificationResult spectrumID=" in line and s==0:
                s=1
                nline=line.strip().split('"')
                outf.write("%s\t" %(nline[1]))
                charge=''
                experimentalMassToCharge=''
                calculatedMassToCharge=''
                peptide=''
                msgfscore=''
                msgfeval=''
                pepname=''
            elif s==1:
                if "<SpectrumIdentificationItem " in line:
                    nline=line.strip().split()
                    for word in nline:
                        if 'chargeState=' in word and charge=='':
                            word=word.split('"')[1]
                            charge=word
                        if 'experimentalMassToCharge=' in word and experimentalMassToCharge=='':
                            word=word.split('"')[1]
                            experimentalMassToCharge=word
                
                        if 'calculatedMassToCharge=' in word:
                            word=word.split('"')[1]
                            calculatedMassToCharge+=word+'|'
    
                        if 'peptide_ref=' in word:
                            pepid=word.split('"')[1]
                            word=pepid.split('_')[1]
                            if pep_details.has_key(pepid):
                                name=pep_details[pepid]['name']
                                pre=pep_details[pepid]['pre']
                                post=pep_details[pepid]['post']
                            else:
                                name='-'
                                pre='-'
                                post='-'
                            peptide+=pre+'.'+word+'.'+post+'|'
                            pepname+=name+';'
    
                elif 'MS-GF:RawScore' in line:
                    nline=line.strip().split()
                    for word in nline:
                        if 'value=' in word:
                            word=word.split('"')[1]
                            msgfscore+=word+'|'
                elif 'MS-GF:SpecEValue' in line:
                    nline=line.strip().split()
                    for word in nline:
                        if 'value=' in word:
                            word=word.split('"')[1]
                            msgfeval+=word+'|'
                elif '</SpectrumIdentificationResult>' in line:
                    s=0
                    outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pepname,charge,experimentalMassToCharge,calculatedMassToCharge,peptide,msgfscore,msgfeval))
                
            