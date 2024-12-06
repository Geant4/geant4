# %%
# To run, change the links and definitions in this script and 
# type in the terminal containing the file 
# $python3 chromosAnalysis.py

# python version tested 3.10.12
# developed in visual studio code (jupyter notebook)
# Author contact, Konstantinos Chatzipapas, konstantinos.chatzipapas@cern.ch, 08/12/2023

# %%
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt

# %%
# Define filenames
iRootFile = "molecular-dna.root"
print("Loaded file: ",iRootFile )

# %%
# Initial parameters (need be inserted manually, if modifications are made in the geometry)
eVtoJ = 1.60218e-19
r3 = 10575e-9 * 3450e-9 * 10575e-9   # a * b * c
mass = 997 * 4 * 3.141592 * r3 / 3   # waterDensity * 4/3 * pi * r3 in kg

# %%
# Length of chromosomes in base pairs (bp) (need be inserted manually, if geometry file is modified)
lc =np.zeros((24,1))
lc[1]  = (335601410)
lc[2]  = (424230628)
lc[3]  = (330358245)
lc[4]  = (214431219)
lc[5]  = (225844963)
lc[6]  = (187075736)
lc[7]  = (104619404)
lc[8]  = (523147774)
lc[9]  = (286308941)
lc[10] = (308642731)
lc[11] = (281141907)
lc[12] = (286680368)
lc[13] = (297671931)
lc[14] = (242439583)
lc[15] = (517729277)
lc[16] = (363900456)
lc[17] = (407437606)
lc[18] = (352426730)
lc[19] = (165113373)
lc[20] = (137667917)
lc[21] = (132362463)
lc[22] = (159987865)
lc[23] = (99027882)

# Calculation of cumulative length and total length of the cell DNA
sum_lc=0
slc =np.zeros((24,1))
for j in range(len(lc)):
    sum_lc += lc[j]
    slc[j] = sum_lc
print("Total length of DNA in the cell(bp): ", sum_lc)
#print(slc) #Length of each chromosome

# %%
# Calculate the damage yield
number = 0
f = root.TFile(iRootFile)
gTree = f.Get("tuples/primary_source")
number += gTree.GetEntries()
    
nEntries = 0
fTree = f.Get("tuples/damage")
nEntries += fTree.GetEntries()
#print (nEntries)

c_DSBBPID = {}
for i in range(1,len(lc)):
    c_DSBBPID[i] = []


SSB  = 0
ttSSB = 0
ttDSB = 0

c_SSB =np.zeros(24)
c_DSB =np.zeros(24)
cc_SSB=np.zeros(24)
cc_DSB=np.zeros(24)

c_SSBp  = np.zeros(24)
c_DSBp  = np.zeros(24)
c_SSBpp = np.zeros(24)
c_DSBpp = np.zeros(24)

DSBd = np.zeros(24)
DSBi = np.zeros(24)
DSBm = np.zeros(24)
DSBh = np.zeros(24)

SSBd = np.zeros(24)
SSBi = np.zeros(24)
SSBm = np.zeros(24)

# Calculate total number of SSB, DSB, as well as initial break positions for the calculation of fragments
DD = 0
for e in fTree:
    if (e.TypeClassification=="DSB" or e.TypeClassification=="DSB+" or e.TypeClassification=="DSB++"):
        DD += 1
        ttDSB += e.DirectBreaks + e.IndirectBreaks
        for i in range(1,len(lc)):
            if (e.BasePair>=slc[i-1] and e.BasePair<slc[i]):
                cc_DSB[i] += 1    #e.DirectBreaks + e.IndirectBreaks
                c_DSBBPID[i].append((e.Event, e.BasePair, e.TypeClassification))


    if (e.TypeClassification=="SSB" or e.TypeClassification=="SSB+" or e.TypeClassification=="2SSB"):
        SSB  += 1
        ttSSB += e.DirectBreaks + e.IndirectBreaks
        for i in range(len(lc)):
            if (e.BasePair>=slc[i-1] and e.BasePair<slc[i]):
                cc_SSB[i] += 1    #e.DirectBreaks + e.IndirectBreaks

#print(DD,SSB)  #Used for testing

# Calculate damage complexity for each chromosome
for e in fTree:
    for i in range(1,len(lc)):
        if (e.BasePair>=slc[i-1] and e.BasePair<slc[i]):
            if e.TypeClassification=="DSB":
                c_DSB[i] += 1
            elif e.TypeClassification=="DSB+":
                c_DSBp[i] += 1
            elif e.TypeClassification=="DSB++":
                c_DSBpp[i] += 1
            elif e.TypeClassification=="SSB":
                c_SSB[i] += 1
            elif e.TypeClassification=="SSB+":
                c_SSBp[i] += 1
            elif e.TypeClassification=="2SSB":
                c_SSBpp[i] += 1
            
for e in fTree:
    for i in range(1,len(lc)):
        if (e.BasePair>=slc[i-1] and e.BasePair<slc[i]):
            if e.SourceClassification=="DSBd":
                DSBd[i] += 1
            elif e.SourceClassification=="DSBi":
                DSBi[i] += 1
            elif e.SourceClassification=="DSBm":
                DSBm[i] += 1
            elif e.SourceClassification=="DSBh":
                DSBh[i] += 1
            elif e.SourceClassification=="SSBd":
                SSBd[i] += 1
            elif e.SourceClassification=="SSBi":
                SSBi[i] += 1
            elif e.SourceClassification=="SSBm":
                SSBm[i] += 1

# Damage classification, this is not needed
repDSB = 0
irrDSB =0
tDSB = 0

gTree = f.Get("tuples/classification")
for e in gTree:
    repDSB += e.DSB
    irrDSB += e.DSBp + e.DSBpp
    tDSB += e.DSB + e.DSBp + 2*e.DSBpp


tSSB = 0
ffTree = f.Get("tuples/source")
for e in ffTree:
    tSSB += e.SSBd + e.SSBi + e.SSBm

totalDSB = repDSB + irrDSB

#print (ttDSB,tDSB,totalDSB,repDSB,irrDSB,"  ",SSB,ttSSB,tSSB,tSSB/ttDSB, tSSB/i)  #Used for testing

acc_edep=0
hTree= f.Get("tuples/chromosome_hits")
for e in hTree:
    acc_edep += (e.e_chromosome_kev + e.e_dna_kev )*1e3

dose = acc_edep * eVtoJ / mass
#print("Accumulated deposited energy in the cell(MeV): ", acc_edep)
#print("Absorbed Dose to the cell(Gy): ", dose)

# %%
#print(c_DSBBPID)
#print(c_DSBBPID[22])

# %%
# Printing damage classification for each chromosome
cSSBsum   =0
cSSBpsum  =0
cSSBppsum =0
cDSBsum   =0
cDSBpsum  =0
cDSBppsum =0
ccSSBsum  =0
ccDSBsum  =0


for i in range(1,len(lc)):
    print("\nchromosome",i," total DSB/Gy/GBp: ", cc_DSB[i]/(lc[i]/1e+9)/dose," --> DSB:", c_DSB[i]/(lc[i]/1e+9)/dose, "DSB+:", c_DSBp[i]/(lc[i]/1e+9)/dose, "DSB++:", c_DSBpp[i]/(lc[i]/1e+9)/dose)
    print(  "chromosome",i," total SSB/Gy/GBp: ", cc_SSB[i]/(lc[i]/1e+9)/dose," --> SSB:", c_SSB[i]/(lc[i]/1e+9)/dose, "SSB+:", c_SSBp[i]/(lc[i]/1e+9)/dose, "SSB++:", c_SSBpp[i]/(lc[i]/1e+9)/dose)
    print(  "chromosome",i," DSBd: ", DSBd[i]/(lc[i]/1e+9)/dose,"DSBi:", DSBi[i]/(lc[i]/1e+9)/dose, "DSBm:", DSBm[i]/(lc[i]/1e+9)/dose, "DSBh:", DSBh[i]/(lc[i]/1e+9)/dose)
    print(  "chromosome",i," SSBd: ", SSBd[i]/(lc[i]/1e+9)/dose,"SSBi:", SSBi[i]/(lc[i]/1e+9)/dose, "DSBm:", SSBm[i]/(lc[i]/1e+9)/dose)
    #print("chromosome",i," total DSB/Gy/GBp: ", cc_DSB[i]/(lc[i]/1e+9)/dose, "  \tchromosome",i," total SSB/Gy/GBp: ", cc_SSB[i]/(lc[i]/1e+9)/dose)
    cDSBsum += c_DSB[i]+c_DSBp[i]+c_DSBpp[i]
    cSSBsum += c_SSB[i]+c_SSBp[i]+c_SSBpp[i]
    cSSBpsum  += c_SSBp[i]
    cSSBppsum += c_SSBpp[i]
    cDSBpsum  += c_DSBp[i]
    cDSBppsum += c_DSBpp[i]
    
    ccDSBsum += cc_DSB[i]
    ccSSBsum += cc_SSB[i]
#print("total DSB/Gy/GBp:",ccDSBsum/(sum_lc/1e+9)/dose, "\ttotal SSB/Gy/GBp:",ccSSBsum/(sum_lc/1e+9)/dose)
print("\nIf the geometry is considered as a whole nucleus:")
#print("total DSB/Gy/GBp:",cDSBsum/(sum_lc/1e+9)/dose, "\ttotal SSB/Gy/GBp:",cSSBsum/(sum_lc/1e+9)/dose)
print("total DSB/Gy/GBp:",cDSBsum/(sum_lc/1e+9)/dose, " --> DSB:", (cDSBsum-cDSBpsum-cDSBppsum)/(sum_lc/1e+9)/dose, "DSB+:", cDSBpsum/(sum_lc/1e+9)/dose, "DSB++:", cDSBppsum/(sum_lc/1e+9)/dose)
print("total SSB/Gy/GBp:",cSSBsum/(sum_lc/1e+9)/dose, " --> SSB:", (cSSBsum-cSSBpsum-cSSBppsum)/(sum_lc/1e+9)/dose, "SSB+:", cSSBpsum/(sum_lc/1e+9)/dose, "SDSB++:", cSSBppsum/(sum_lc/1e+9)/dose)
print("SSB/DSB ratio: ", cSSBsum/cDSBsum)
print("Total dose (Gy): ", dose)
#print(cDSBsum,ccDSBsum,cSSBsum,ccSSBsum)

print("\nThank you for using molecularDNA and supporting the Geant4-DNA collaboration.!")
