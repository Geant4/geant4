#!/usr/bin/env python
# coding: utf-8

# In[1]:


# To run, change the links and definitions in this script and 
# type in the terminal containing the file 
# $python3 createSDD.py

# Author contact, Konstantinos Chatzipapas, chatzipa@cenbg.in2p3.fr, 14/03/2022


# In[2]:


import ROOT as root
import numpy as np


# In[3]:


# Define filenames
outputFile = "SDD_minFormatMolecularDNA.txt"
iMacFile   = "human_cell.mac"
iRootFile  = "molecular-dna.root"


# In[4]:


# Definitions
# They may be adapted in the code, if such information exists in the simulation definition file

# Define if it is a simulation of a single-track irradiation (0), a delivered dose (1) or a fluence (2)
doseFluence = 1
amountDoseFluence = 0

# Define dose rate
doseRate = float(0)

# Define the number of chromosomes and the length in MBp
nbChromo = 1
chromoLength = 6405.886128

# Define DNA density (BPs/um3)
dnaDensity = 1.2*1e07

# Define cell cycle phase (use G0 (1), G1 (2), S (3), G2 (4) and M (5), together with percentage, e.g [3,0.7] indicates a cell 70% of the way through S phase)
cellCyclePhase = "0, 0.0"

# Define DNA structure :  whole nucleus (0), a heterochromatin region (1), 
                        # euchromatin region (2), 
                        # a mixed (heterochromatin and euchromatin) region (3), 
                        # single DNA fiber (4), DNA wrapped around a single histone (5), 
                        # DNA plasmid (6) or a simple circular (7) or straight (8) DNA section.
dnaStructure = 4
nakedWet = 1            # naked (0), wet (1)

# Define if the real experiment was in-vitro (0) or in-vivo (1)
vitroVivo = 0

# Define the proliferation status, quiescent (0) or proliferating (1)
proliferation = 0

# Define the microenvironment
temperature = 27  # degrees
oxygen = 0.0      # molarity

# Damage definition
directIndirect = 1   # direct effects only (0) or including chemistry (1)
bpNm = 0             # BPs (0) or in nm (1)
bpThres = 10.0       # the distance in BPs or nm between backbone lesions that are considered DSBs
baseLesions = -1     # This value then determines the distance (in BP or nm) beyond the outer backbone damages where base damages are also stored in the same site (float).

# Default, if nothing exists in mac file
sourceTropus = "iso"


# In[5]:


# Initialize several parameters
title    = "No Title Included"
incidentParticles = " "
particle = " "
energy   = " "
energyUnits = " "
sourceShape = " "
sourceTropus = " "
sourceX = "0"
sourceY = "0"
sourceZ = "0"
sourceUnits = "nm"
targetType  = " "
targetShape = " "
targetRadius = 0
targetX  = 0
targetY  = 0
targetZ  = 0
r3 = 0
worldSize = 0
directDamageL = " "
directDamageU = " "
time = " "
timeUnit = " "



# In[6]:


with open(iMacFile, 'r') as infile:
    for line in infile:
        spl = line.split(" ")
        if spl[0] == "###":
            title = spl[1]
            #print(title)
        if spl[0] =="/run/beamOn":
            incidentParticles = int(spl[1])
            #print(incidentParticles)
        if spl[0] =="/gps/particle":
            spl2 = spl[1].split("\n")
            particle = spl2[0]
            #print(particle)
        if spl[0] =="/gps/energy":
            sourceType = 1
            energyDistribution = "M, 0"
            particleFraction = float(1)
            energy = float(spl[1])
            spl2 = spl[2].split("\n")
            energyUnits = spl2[0]
            #print(energy,energyUnits)
        if spl[0] =="/gps/pos/shape":
            spl2 = spl[1].split("\n")
            sourceShape = spl2[0]
            #print(sourceShape)
        if spl[0] =="/gps/ang/type":
            spl2 = spl[1].split("\n")
            sourceTropus = spl2[0]
            #print("type",sourceTropus)
        if spl[0] =="/gps/pos/centre":
            sourceX = spl[1]
            sourceY = spl[2]
            sourceZ = spl[3]
            spl2 = spl[4].split("\n")
            sourceUnits = spl2[0]
        if spl[0] =="/chromosome/add":
            targetType = spl[1]
            targetShape = spl[2]
            if targetShape == "sphere":
                #print (spl)
                if spl[4]=="nm\n":
                    targetRadius = float(spl[3])/1000
                    r3 = float(spl[3])*1e-09
                    #print ("r3=", r3)
                elif spl[4]=="um\n":
                    targetRadius = float(spl[3])
                    r3 = float(spl[3])*1e-06
                    #print ("r3=", r3)
            if targetShape == "ellipse":
                #print (spl)
                if spl[9]=="nm":
                    targetX = float(spl[3])/1000
                    targetY = float(spl[4])/1000
                    targetZ = float(spl[5])/1000
                    r3 = float(spl[3])*1e-09*float(spl[4])*1e-09*float(spl[5])*1e-09
                    #print ("r3=", r3)
                elif spl[9]=="um":
                    targetX = float(spl[3])
                    targetY = float(spl[4])
                    targetZ = float(spl[5])
                    r3 = float(spl[3])*1e-06*float(spl[4])*1e-06*float(spl[5])*1e-06
                    #print ("r3=", r3)
            #print(targetType,targetShape,targetRadius,targetX,targetY,targetZ)
        if spl[0] =="/world/worldSize":
            worldSize = float(spl[1])/2
            if spl[2] =="nm\n":
                worldSize = worldSize/1000
            #print(worldSize)
        if spl[0] =="/dnadamage/directDamageLower":
            directDamageL = spl[1]
            #print(directDamageL)
        if spl[0] =="/dnadamage/directDamageUpper":
            directDamageU = spl[1]
            #print(directDamageU)
        if spl[0] =="/scheduler/endTime": 
            time = spl[1]
            spl2 = spl[2].split("\n")
            timeUnit = spl2[0]
            #print(time, timeUnit)


# In[7]:


# Creating SDD file
with open(outputFile, 'w') as ofile:
# Starting with header part of SDD file    
    ofile.write("\nSDD version, SDDv1.0;\n")
    ofile.write("Software, MolecularDNA;\n")
    ofile.write("Author contact, Konstantinos Chatzipapas, chatzipa@cenbg.in2p3.fr, "
                "14/03/2022;\n")
    ofile.write("***Important information*********************************************\n"
                "To provide some extra information on the quantification of DSB, "
                "DSB+ and DSB++, the last column of data section, includes the values "
                "of 4, 5 that correspond to DSB+ and DSB++ respectively;\n"
                "*********************************************************************\n")
    #Adjust description
    ofile.write("Simulation Details, "+title+". DNA damages from direct and indirect effects;\n")
    ofile.write("Source, Monoenergetic "+sourceShape+" "+particle+" "+sourceTropus+"tropic,")
    ofile.write(" centered at "+sourceX+" "+sourceY+" "+sourceZ+" "+sourceUnits)
    ofile.write(" of a cell nucleus, containing a free DNA segment. Energy: ")
    ofile.write(str(energy)+" "+energyUnits+";\n") 
    ######
    ofile.write("Source type, "+str(sourceType)+";\n")   # Needs improvement
    ofile.write("Incident particles, "+str(incidentParticles)+";\n")
    ofile.write("Mean particle energy, "+str(energy)+" MeV;\n")
    ofile.write("Energy distribution, "+energyDistribution+";\n")   # Needs improvement
    ofile.write("Particle fraction, "+str(particleFraction)+";\n")   # Needs improvement
    
    f = root.TFile(iRootFile)
    
    eVtoJ = 1.60218e-19

    mass = 997 * 4 * 3.141592 * r3 / 3     # density * 4/3 * pi * r3

    acc_edep = 0
    dose = 0
    nbEntries = 0
    ffTree = f.Get("tuples/chromosome_hits")
    nbEntries += ffTree.GetEntries()
    for ev in ffTree:
        acc_edep += (ev.e_chromosome_kev + ev.e_dna_kev) *1e3  # eV


    # Calculate the absorbed dose
    amountDoseFluence = acc_edep * eVtoJ / mass    # Dose in Gy
    
    
    
    ofile.write("Dose or fluence, "+str(doseFluence)+", "+str(amountDoseFluence)+";\n")   # Needs improvement
    ofile.write("Dose rate, "+str(doseRate)+";\n")
    ofile.write("Irradiation target, Simple free DNA fragment in a "+targetShape+" with ")
    if targetShape == "sphere":
        ofile.write("radius "+str(targetRadius)+" um;\n")
    else: #targetShape == "ellipse":
        ofile.write("dimensions "+str(targetX)+" "+str(targetY)+" "+str(targetZ)+" um;\n")
    ofile.write("Volumes, 0,"+str(worldSize)+","+str(worldSize)+","+str(worldSize)+","+str(-worldSize)+","+str(-worldSize)+","+str(-worldSize)+",")
    if targetShape == "sphere":
        ofile.write(" 1,"+str(targetRadius)+","+str(targetRadius)+","+str(targetRadius)+";\n")
    elif targetShape == "ellipse":
        ofile.write(" 1,"+str(targetX)+","+str(targetY)+","+str(targetZ)+";\n")
    else:
        ofile.write(" 0,"+str(targetX)+","+str(targetY)+","+str(targetZ)+","+str(+targetX)+","+str(+targetY)+","+str(+targetZ)+";\n")
    ofile.write("Chromosome sizes, "+str(nbChromo)+", "+str(chromoLength)+";\n")
    ofile.write("DNA Density, "+str(dnaDensity)+";\n")
    ofile.write("Cell Cycle Phase, "+cellCyclePhase+";\n")
    ofile.write("DNA Structure, "+str(dnaStructure)+", "+str(nakedWet)+";\n")
    ofile.write("In vitro / in vivo, "+str(vitroVivo)+";\n")
    ofile.write("Proliferation status, "+str(proliferation)+";\n")
    ofile.write("Microenvironment, "+str(temperature)+", "+str(oxygen)+";\n")
    ofile.write("Damage definition, "+str(directIndirect)+", "+str(bpNm)+", "+str(bpThres)+", "+str(baseLesions)+", "+str(directDamageL)+";\n")
    ofile.write("Time, "+str(time)+" "+timeUnit+";\n")

    number = 0
    gTree = f.Get("tuples/primary_source")
    number += gTree.GetEntries()
    
    nEntries = 0
    fTree = f.Get("tuples/damage")
    nEntries += fTree.GetEntries()
    #print (nEntries)
    
    ofile.write("Damage and primary count, "+str(nEntries)+", "+str(number)+";\n")
    ofile.write("Data entries, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n")
    ofile.write("Data was generated in minimal output format and used as an example. ")
    ofile.write("Please modify description to match different conditions;\n")
    ofile.write("\n***EndOfHeader***;\n\n")

#############################################################################################    
# Writing the data part of SDD file    
    currentEvent = 0
    Primaryflag = True
    Eventflag = False

    for entry in fTree:
        if entry.Event!=currentEvent:
            Eventflag = True
        currentEvent = entry.Event
        #print(entry.SourceClassification)
        DSB = 0
        if (entry.TypeClassification == "DSB"):
            if (entry.SourceClassification == "DSBd"):
                DSB = 1
            elif (entry.SourceClassification == "DSBi"):
                DSB = 2
            elif (entry.SourceClassification == "DSBh" or entry.SourceClassification == "DSBm"):
                DSB = 3
                #print (entry.TypeClassification, entry.SourceClassification)
        elif (entry.TypeClassification == "DSB+"):
            DSB = 4
        elif (entry.TypeClassification == "DSB++"):
            DSB = 5
            #if (entry.SourceClassification == "DSBh" or entry.SourceClassification == "DSBm"):
                #print (entry.TypeClassification, entry.SourceClassification)
        else:
            DSB = 0
            
        #print(DSB)
        if Primaryflag:
            #print("True")
            ofile.write("2, "+str(entry.Event)+"; ")
            ofile.write(str(entry.Position_x_um)+", "+str(entry.Position_y_um)+", "+str(entry.Position_z_um)+"; ")
            ofile.write(str(entry.BaseDamage)+", "+str(entry.StrandDamage)+", "+str(DSB)+";\n")
            Primaryflag = False
            Eventflag = False
        else:
            if Eventflag:
                ofile.write("1, "+str(entry.Event)+"; ")
                ofile.write(str(entry.Position_x_um)+", "+str(entry.Position_y_um)+", "+str(entry.Position_z_um)+"; ")
                ofile.write(str(entry.BaseDamage)+", "+str(entry.StrandDamage)+", "+str(DSB)+";\n")
                Eventflag = False
            else:
                ofile.write("0, "+str(entry.Event)+"; ")
                ofile.write(str(entry.Position_x_um)+", "+str(entry.Position_y_um)+", "+str(entry.Position_z_um)+"; ")
                ofile.write(str(entry.BaseDamage)+", "+str(entry.StrandDamage)+", "+str(DSB)+";\n")


# In[8]:

print("Output file: ", outputFile)



