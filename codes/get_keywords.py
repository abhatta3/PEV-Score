import heapq
from utility import get_position_mass

def getk(bit,pep,l=4):
    ml=len(pep)
    ms=0
    kw=pep[-l:]
    mpos=ml-l
    for i in range(0,ml-l):
        ts=0
        for j in range(i,i+l):
            ts+=int(bit[j])
        if ts>=ms:
            kw=pep[i:i+l]
            ms=ts
            mpos=i
    r=(kw,mpos)
    return r

def get_key(inputfile,scanfile,ms_dict,fixed_mod,klength):
    #ignore equivalent mass residues
    #ignore  phosphorylation residue serine
    #replace_d={'GG':2,'GA':2,'GV':2,'GE':2,'AD':2,'SV':2,'SS':2,'W':1,'N':1,'Q':1,'K':1,'R':1}
    #pair_replace={'N':'GG','Q':'GA', 'K':'GA', 'R':'GV', 'W':'GE,AD,SV', 'C':'SS'}
    #aa_replace={'M': 'D,V', 'N':'L,M', 'Q':'K', 'R':'L','K':'A,S,N', 'E':'D,F', 'D':'T', 'C':'H'}
    #out_kw_f=open("keywords.txt","w") 
    outf = []
    delta=0.5
    mapscan={}
    eventdic={}
    peptide={}
    ev_ms={}
    for line in inputfile:
        if '#' in line:
            continue
        if len(line)<1:
            break
        line=line.strip().split("\t")
        mapscan[int(line[1])]=line[0]
        if '.' in line[3][1] and '.' in line[3][-2]:
            eventdic[line[0]]=line[3][1:-1].translate(None, '1234567890:._+-*!@#$?')
        else:
            eventdic[line[0]]=line[3].translate(None, '1234567890:._+-*!@#$?')
        peptide[eventdic[line[0]]]=line[3]
    scan_i=-1
    eventsupport={}
    with open(scanfile,"r") as sf:
        matched=False
        for line in sf:
           if 'BEGIN IONS' in line:
               pepmass=0.0
               pcharge=1.0
               scan_i+=1
               scand={}
               slist=[]
               stemp={}
               scount=0
               if mapscan.has_key(scan_i):
                   matched=True
                   event=mapscan[scan_i]
                   pep=eventdic[event]
               else:
                   matched=False
           elif 'END IONS' not in line and matched==True and '=' in line:
               if "PEPMASS" in line:
                   tempmass=line.strip().split('=')[1].split()
                   tempmassdata=tempmass[0]
                   pepmass=float(tempmassdata)
               elif "CHARGE" in line:
                   pcharge=float(line.strip().split('=')[1].split('+')[0])
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
               adjusted_mass=pepmass*pcharge-pcharge+1
               if eventsupport.has_key(event):
                   eventsupport[event].append(pd)
               else:
                   eventsupport[event]=[pd]
               if ev_ms.has_key(event):
                   if ev_ms[event]<adjusted_mass:
                       ev_ms[event]=adjusted_mass
               else:
                   ev_ms[event]=adjusted_mass
    
    outf.append("#Event_ID\tInput_Peptide(I)\tInput_Peptide(II)\tParent_Mass\tSupport_bits\tKey_Words\n")
    for event in eventsupport:
        pep=eventdic[event]
        last=len(pep)
        pbit={}
        for i in range(0,last):
            pbit[i]=0
        adjusted_mass=ev_ms[event]
        for support in eventsupport[event]:
            for j in support:
                pbit[j]+=support[j]
        (kw,peppos)=getk(pbit,pep,klength)
        #print pep,kw,peppos,pbit  
        sup_bit=''
        for i in range(0,last):
            sup_bit+=str(pbit[i])+','

        #ends here
        outf.append("%s\t%s\t%s\t%f\t%s\t%s\t%d\n" %(event,peptide[pep],pep,adjusted_mass,sup_bit,kw,peppos))
        #print("%s\t%s\t%s\t%f\t%s\t%s\t%d\n" %(event,peptide[pep],pep,adjusted_mass,sup_bit,kw,peppos))
        #out_kw_f.write("%s\t%s\t%d\t%s\t%d\n" %(pep,kw,peppos,sup_bit,len(pep)-len(sup_bit)))

    return outf

        
    
