def getannotdictionary(afile):
    with open(afile,"r") as af:
        d={}
        for line in af:
            if '#' in line:
                continue
            line=line.strip().split("\t")
            d[line[2]]=line[3]
    return d	

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
                        if inv<pm:# and position_mass_b[i-1]==1:
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
                        position_mass_t[inext]+=1

        for i in range(0,s-1):
            val=position_mass_y[i]
            if val > 0:
                inext=i+1
                if position_mass_y.has_key(inext):
                    nval=position_mass_y[inext]
                    if nval > 0:
                        position_mass_t[i]+=1
    #print pep, position_mass_t
            
    return  position_mass_t
