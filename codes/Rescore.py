# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:15:18 2017
@author: anindya
"""
import math

def get_list_strig(bm_seq):
    m_seq=''
    if len(bm_seq)<1:
        m_seq=';'
    else:
        for w in bm_seq:
            if len(w)>1:
                m_seq+=w+';'
    return m_seq

def fmod(seq,itq):
    if not itq.has_key('*'):
        val=''
    else:
        v=float(itq['*'])
        if v>0:
            val='+'+str(v)
        else:
            val=str(v)
 
    for l in seq:
        if itq.has_key(l):
            v=float(itq[l])
            if v>0:
                v_l='+'+str(v)
            else:
                v_l=str(v)
            
            val=val+v_l
        else:
            val=val+l      
    return val

    
def rescore(infile,outfile,scan_map_file,singlemgf,novel_list,key_out,fixed_mod,ref_list):
    support={}
    for line in key_out:
        if '#' in line:
            continue
        line=line.strip().split("\t")
        support[line[0]]=line[4]

    
    event_list={}
    event_pep={}
    for line in novel_list:
        if '#' in line:
            continue
        line=line.strip().split("\t")
        if len(line)<2:
            continue
        event_list[line[0]]=line[1]
        event_pep[line[0]]=line[2]

    scan_event={}
    for line in scan_map_file:
        if '#' in line:
            continue
        line=line.strip().split("\t")
        if len(line)<2:
            continue
        scan_event[line[1]]=line[0]

    output=open(outfile,"w")
    outputi=open(outfile+".index","w")

    reject={}
    event_dic={}    
    with open(infile,"r") as inf:
        for line in inf:
            rline=line
            line=line.strip().split("\t")

            index=line[0]
            indval=line[0].split('=')[1]
            evid=scan_event[indval]         
            nf=0
            mf=0
            ni=-1
            mi=-1
            uid='NA'
            event=evid
            for i,word in enumerate(line[1].split(';')):
                if evid+'|Novel' in word and nf==0:
                    nf=1
                    ni=i
                if 'Modified' in word and mf==0:
                    mf=1
                    mi=i
                    uid=word.split('|')[0]
                if mf==1 and nf==1:
                    break
            score=line[-1].split('|')
            try:
                if ni==-1:
                    novelscore=1.0
                else:
                    novelscore=float(score[ni])
                if mi==-1:
                    modscore=1.0
                else:
                    modscore=float(score[mi])
            except:
                reject[index]=rline
                continue
            if novelscore<0.01 or modscore<0.01:
                sc=float(-1.0*float(math.log(novelscore/modscore)))
            else:
                reject[index]=rline
                continue
                
            seq=line[5].split('|')
            if ni==-1:
                nov_seq='NA'
            else:
                nov_seq=seq[ni].translate(None, '[')
            if mi==-1:
                mod_seq='NA'
            else:
                mod_seq=seq[mi].translate(None, '[')
            mss=line[4].split('|')
            if ni==-1:
                nov_mss=1
            else:
                nov_mss=mss[ni]
    
            if mi==-1:
                mod_mss=1
            else:
                mod_mss=mss[mi]
    
            better_mod=''
            a_bmseqid=''
    
            if novelscore<modscore:
                new_nb=1
                new_mb=0
                new_eqb=0
            elif novelscore>modscore:
                new_nb=0
                new_mb=1
                new_eqb=0
                better_mod=mod_seq
                a_bmseqid=uid
            elif novelscore==modscore:
                new_nb=0
                new_mb=0
                new_eqb=1
            if event_dic.has_key(event):
                oldscore=event_dic[event]['pev']
                nb_count=event_dic[event]['nb']+new_nb # novel best
                mb_count=event_dic[event]['mb']+new_mb # mod best
                eqb_count=event_dic[event]['eqb']+new_eqb # nov=mod
                bmseq=event_dic[event]['bm_seq']        # mod sequences better than nov
                bmseqid=event_dic[event]['bm_seqid']
                if len(better_mod)>3:
                    if len(bmseq)>3:
                            bmseq+=better_mod+'|'
                            bmseqid+=a_bmseqid+'|'
                    else:
                        bmseq=better_mod+'|'
                        bmseqid=a_bmseqid+'|'
                if sc>oldscore:
                    event_dic[event]={'ind':index,'nseq':nov_seq,'mseq':mod_seq,'pev':sc,'neval':novelscore,'meval':modscore,'charge':line[2],'indexmass':line[3],'nmss':nov_mss,'mmss':mod_mss,'nb':nb_count,'mb':mb_count,'eqb':eqb_count,'bm_seq':bmseq,'uid':uid,'bm_seqid':bmseqid}
                else:
                    event_dic[event]['nb']=nb_count
                    event_dic[event]['mb']=mb_count
                    event_dic[event]['eqb']=eqb_count
                    event_dic[event]['bm_seq']=bmseq
                    event_dic[event]['bm_seqid']=bmseqid
            else:
                if len(better_mod)>1:
                    bmseq=better_mod+'|'
                    bmseqid=a_bmseqid+'|'
                else:
                    bmseq='NA'
                    bmseqid='NA'
                    
                event_dic[event]={'ind':index,'nseq':nov_seq,'mseq':mod_seq,'pev':sc,'neval':novelscore,'meval':modscore,'charge':line[2],'indexmass':line[3],'nmss':nov_mss,'mmss':mod_mss,'nb':new_nb,'mb':new_mb,'eqb':new_eqb,'bm_seq':bmseq,'uid':uid,'bm_seqid':bmseqid}

    output.write("#Event_ID\tNovel_Peptide\tRescored_Novel_Peptide\tAlternative_Reference_Sequence\tUniProt_ID_Alternative_Reference\tVariant_Type\tPEV_Score\tPEV_Category\tAll_Higher_Scored_Alternative_Reference\tAll_Higher_Scored_Alternative_Reference_ProteinID\n")
    outputi.write("#Event\tNovel_peptide_seq\tAlternative_reference_peptide_seq\tindex\n")

    for event in event_dic:
        nseq=event_dic[event]['nseq']
        mseq=event_dic[event]['mseq']
        pev=event_dic[event]['pev']
        bm_seq=event_dic[event]['bm_seq']
        bmseqid=event_dic[event]['bm_seqid']
        auid=event_dic[event]['uid']
        bm_seq=list(set(bm_seq.split('|')))
        bmseqid=list(set(bmseqid.split('|')))
        m_seq=get_list_strig(bm_seq)
        mid_seq=get_list_strig(bmseqid)
        if pev<0:
            cl='PEV-'
        elif pev==0:
            cl='PEVzero'
        else:
            cl='PEV+'
            mseq='NA'
            auid='NA'
            
        if len(m_seq)<5 or cl=='PEV+':
            m_seq=''
            mid_seq=''

        if cl=='PEV-' or cl=='PEVzero':
            ind=str(int(event_dic[event]['ind'].split('=')[1])+1)
            if nseq!='NA':
                i_nseq=fmod(nseq[2:-2],fixed_mod)
            else:
                i_nseq='NOTFOUND'
    
            if mseq!='NA':            
                i_mseq=fmod(mseq[2:-2],fixed_mod)
            else:
                i_mseq='NOTFOUND'
            outputi.write("%s\t%s\t%s\t%s\n" %(event,i_nseq,i_mseq,ind))
            

        output.write("%s\t%s\t%s\t%s\t%s\t%s\t%0.4f\t%s\t%s\t%s\n" %(event,event_pep[event],nseq,mseq,auid,event_list[event],pev,cl,m_seq,mid_seq))
        

    for rline in ref_list[1:]:
        rline=rline.strip().split("\t")
        output.write("%s\t%s\t%s\t%s\t%s\t%s\tNA\tPEVexact-reference\t\t\n" %(rline[1],rline[3],rline[3],rline[3],rline[0],rline[2]))

    for ev in event_pep:
        if not event_dic.has_key(ev):
            output.write("%s\t%s\tNA\tNA\tNA\t%s\tNA\tPEVunscored\t\t\n" %(ev,event_pep[ev],event_list[ev]))

    output.close()
    
    
        
    
        
    
