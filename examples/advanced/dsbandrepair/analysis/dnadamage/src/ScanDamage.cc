//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file ScanDamage.cc
/// \brief Implementation of the ScanDamage class

#include "ScanDamage.hh"

#include "TSystemDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
TRandom gRandomGen; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ScanDamage::ScanDamage(): fSkipScanningIndirectDamage(false)
{
    fThresholdEnergy = 17.5; //eV
    fEdepSumInNucleus = 0;
    fProbabilityForIndirectSB = 0.40; // 40%
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > ScanDamage::ExtractDamage(){
    ReadCellandVoxelDefFilePaths();
    fMergedTables.clear();
    fDamage.clear();
    RetrieveVoxelBp();
    FillVoxelData();
    ScanDamageFromPhys();
    if (!fSkipScanningIndirectDamage) ScanDamageFromChem();
    MergeDamageFromPhysChem();
    for (const auto& [pChrom,table] : fMergedTables) {
        unsigned int evt;
        unsigned int strand;
        ullint cpyNb;
        unsigned int isBase;
        double time;
        double edep;
        double origin;
        Damage::DamageType pType;
        Damage::DamageCause pOrigin;
        std::map<unsigned int,std::vector<Damage> > perChromoDamage;
        for (const auto& v : table) {
            evt = v.at(0);
            strand = v.at(1);
            cpyNb = v.at(2);
            isBase = v.at(3);
            time = v.at(4);
            edep = v.at(5);
            origin = v.at(6);
            if(isBase==0)
                pType=Damage::DamageType::fBackbone;
            else
                pType=Damage::DamageType::fBase;
            if(origin==0)
                pOrigin=Damage::DamageCause::fDirect;
            else
                pOrigin=Damage::DamageCause::fIndirect;
            if (perChromoDamage.find(evt) == perChromoDamage.end()) {
                std::vector<Damage> dm{Damage(pType,pChrom,evt,strand,cpyNb,Position(0,0,0),
                                            pOrigin,Damage::DamageChromatin::fUnspecified)};
                perChromoDamage.insert({evt,dm});
            } else perChromoDamage[evt].push_back(Damage(pType,pChrom,evt,strand,cpyNb,Position(0,0,0),
                                            pOrigin,Damage::DamageChromatin::fUnspecified));
        }
        fDamage.insert({pChrom,perChromoDamage});
    }
    return fDamage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::RetrieveVoxelBp()
{
    fBpPerVoxel.clear();
    
    for (const auto & entry : fVoxelDefFilesList) {
        std::ifstream file(entry);
        if(!file.good() )
        {
            std::cerr<<"**** Fatal Error *****"<<std::endl;
            std::cerr<<"ScanDamage::RetrieveVoxelBp: No file named "<<entry<<std::endl;
            std::cerr<<"*************** *****"<<std::endl;
            exit(EXIT_FAILURE);
        }
        std::string voxelName = "noName";
        std::string line;
        bool foundName = false;
        bool foundNumOfBp = false;
        while(std::getline(file, line)
            && !foundNumOfBp)
        {
            std::istringstream iss(line);
            std::string flag;
            iss >> flag;
            std::string charac;
            iss >> charac;
            // Look for the name of the voxel
            if(flag=="_Name")
            {
                voxelName = charac;
                foundName = true;
            }
            // Look for the flag "_Number"
            // And the characteristic "voxelBasePair"
            if(flag=="_Number" && charac=="voxelBasePair")
            {
                int numOfBp;
                iss >> numOfBp;
                if(!foundName)
                {
                    std::cerr<<"*** Fatal Error ***"<<std::endl;
                    std::cerr<<"ScanDamage::RetrieveVoxelBp: The number of bp was found before the name "
                    <<"of the voxel... This is an unexpected case."<<std::endl;
                    std::cerr<<"******"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                else
                {
                    fBpPerVoxel[voxelName] = numOfBp;
                    std::cout<<voxelName<<" has "<<numOfBp<<" bp"<<std::endl;
                    foundNumOfBp = true;
                }
            }
        }
        file.close();
    }
    
    if (fBpPerVoxel.size() == 0) {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"ScanDamage::RetrieveVoxelBp: No Bp found in voxel definition files. \n Or make"
        <<" sure that file imp.info exists in working directory!"<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::FillVoxelData()
{   
    std::ifstream file(fCellDefFilePath);
    if(!file.good() )
    {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"FillVoxelData: No file named "<<fCellDefFilePath<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
    ullint bpCount = 0;
    unsigned int voxelCount = 0;
    int chromo_previous = 0;
    // Read the file line by line
    std::string line;
    while(std::getline(file, line) )
    {
        std::istringstream iss(line);
        std::string flag;
        iss >> flag;
        // If the flag correspond to the placement of a voxel
        if(flag == "_pl")
        {
            std::string voxelName;
            iss >> voxelName;
            int chromo;
            iss >> chromo;
            int domain;
            iss >> domain;
            // If we change of chromosome then reset the number of bp.
            // Each chromosome starts at 0 bp.
            if(chromo != chromo_previous)
            {
                bpCount = 0;
                chromo_previous = chromo;
            }
            // Fill the data structure
            fVoxels.push_back( VoxelData(chromo, domain, bpCount) );
            int numBpinthisvoxel = int(fBpPerVoxel[voxelName]);
            bpCount += numBpinthisvoxel;
            if (fChromosomeBpMap.find(chromo) == fChromosomeBpMap.end()) {
                fChromosomeBpMap.insert({chromo,numBpinthisvoxel});
            } else {
                fChromosomeBpMap[chromo] += numBpinthisvoxel;
            }
            voxelCount++;
        }
    }
    file.close();
    

    if (fVoxels.size() == 0) {
        std::cerr<<"**** Fatal Error *****"<<std::endl;
        std::cerr<<"ScanDamage::FillVoxelData: NofVoxels info found in files "<<fCellDefFilePath<<std::endl;
        std::cerr<<"*************** *****"<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout<<"Num of voxels: "<< fVoxels.size()<<" placed in Cell Nucleus."<<std::endl;
    std::cout<<"=====Choromosome sizes====="<<std::endl;
    std::cout<<"Chromosome ID\t Number of Bp"<<std::endl;
    for (auto const& [chrom, nBp] : fChromosomeBpMap) {
         std::cout<<chrom<< "\t"<<nBp<<std::endl;
    }
    std::cout<<"==========================="<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ScanDamageFromPhys()
{
    std::cout<<"===== Start Scanning Damages From Phys =====\n";
    fEdepSumInNucleus = 0;
    fphysTables.clear();
    fphysSlectedTables.clear();
    fs::path currentP{"phys_output"};
    fs::file_status s = fs::file_status{};
    auto isExist = fs::status_known(s) ? fs::exists(s) : fs::exists(currentP);
    if (isExist) {
        bool isFoundRootFiles = false;
        for (const auto entry : fs::directory_iterator(currentP)) {
            if (entry.path().extension() == ".root") {
                std::cout <<"ScanDamageFromPhys(): Processing file: "<< entry.path().filename()<< std::endl;
                AnaPhysRootFile(entry.path());
                if (!isFoundRootFiles) isFoundRootFiles=true;
            }
        }
        if (!isFoundRootFiles) {
            std::cout<<"=====>> No root files found in folder \"phys_ouput\"!!! Skip Scanning Damages From Phys =====\n";
        }
        if (fphysTables.size() > 0) {
            SortPhysTableWithSelection();
        }
    } else {
        std::cout<<"=====>> Cannot find folder \"phys_ouput\"!!! Skip Scanning Damages From Phys =====\n";
    }

    std::cout<<"===== End Scanning Damages From Phys =====\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ScanDamageFromChem()
{
    std::cout<<"===== Start Scanning Damages From Chem =====\n";
    fs::path currentP{"chem_output"};
    fs::file_status s = fs::file_status{};
    auto isExist = fs::status_known(s) ? fs::exists(s) : fs::exists(currentP);
    if (isExist) {
        bool isFoundRootFiles = false;
        for (const auto entry : fs::directory_iterator(currentP)) {
            if (entry.path().extension() == ".root") {
                AnaChemRootFile(entry);
                if (!isFoundRootFiles) isFoundRootFiles=true;
            }
        }
        if (!isFoundRootFiles) {
            std::cout<<"=====>> No root files found in folder \"chem_ouput\"!!! Skip Scanning Damages From Chem =====\n";
            fSkipScanningIndirectDamage = true;
        }
        if (fchemTables.size() > 0) {
            SortChemTableWithSelection();
        }
    } else {
        std::cout<<"=====>> Cannot find folder \"chem_ouput\"!!! Skip Scanning Damages From Chem =====\n";
        fSkipScanningIndirectDamage = true;
    }
    
    std::cout<<"===== End Scanning Damages From Chem =====\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootFile(const std::string fileName)
{
    TFile* f = new TFile(fileName.c_str());
    if(f->IsZombie() ){
        // File is corrupted
        std::cerr<<"*********** Warning *************"<<std::endl;
        std::cerr<<"The file "<<fileName<<" seems to be corrupted..."<<std::endl;
        std::cerr<<"We will skip it."<<std::endl;
        std::cerr<<"**********************************"<<std::endl;
        return;
    }
    AnaPhysRootTree1(f);
    AnaPhysRootTree2(f);
    f->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaChemRootFile(fs::directory_entry entry)
{
    std::string fnameWithoutExtension = entry.path().stem().string();
    auto [eventNumber, voxelNumber] = GetEventNberAndVoxelNberFromChemRoot(fnameWithoutExtension);
    
    // Voxel data retrieval
    VoxelData& voxelData = fVoxels.at(voxelNumber);

    int chromo = voxelData.fChromosome;
    ullint firstBpNum = voxelData.fFirstBpCopyNum;
    // *******************
    // Analyse of the ntuple to detect SB and associate them with a bpNumCorrected
    // *******************

    // Load the file, the directory and the ntuple
    TFile f(entry.path().c_str());
    if(f.IsZombie() ){
        corruptedFiles++;
        // File is corrupted
        std::cerr<<"*********** Warning *************"<<std::endl;
        std::cerr<<"The file "<<entry.path().string()<<" seems to be corrupted..."<<std::endl;
        std::cerr<<"We will skip it."<<std::endl;
        std::cerr<<"Number of corrupted files: "<< corruptedFiles<<std::endl;
        std::cerr<<"**********************************"<<std::endl;
    } else {
        TDirectoryFile *d = dynamic_cast<TDirectoryFile*> (f.Get("ntuple") );
        TTree* chemTree = (TTree*) d->Get("ntuple_2");

        if( (int) chemTree->GetEntries() >0)
        {
            double strand;
            double copyNumber;
            double xp;
            double yp;
            double zp;
            double time;
            double base;
            chemTree->SetBranchAddress("strand", &strand);
            chemTree->SetBranchAddress("copyNumber", &copyNumber);
            chemTree->SetBranchAddress("xp", &xp);
            chemTree->SetBranchAddress("yp", &yp);
            chemTree->SetBranchAddress("zp", &zp);
            chemTree->SetBranchAddress("time", &time);
            chemTree->SetBranchAddress("base", &base);
            unsigned int entryNumber = (int) chemTree->GetEntries();
            for (unsigned int e=0;e<entryNumber;e++)
            {
                chemTree->GetEntry(e);
                ullint cpNumCorrected = firstBpNum+int(copyNumber);
                std::vector<ullint> newLine{ 
                    (ullint)eventNumber,(ullint)strand,cpNumCorrected,
                    (ullint)base,(ullint)(time*1000000)};
                auto itr = fchemTables.find(chromo);
                if ( itr == fchemTables.end()) {
                    Table tbforThisChro{newLine};
                    fchemTables.insert({chromo,tbforThisChro});
                } else (itr->second).push_back(newLine);
            }
        }
    }//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::tuple<unsigned int, unsigned int> ScanDamage::GetEventNberAndVoxelNberFromChemRoot(
    const std::string fileNameWithoutExtension)
{
    unsigned int evnN, volxelN;
    auto fristPos = fileNameWithoutExtension.find_first_of("_");
    auto secondPos = fileNameWithoutExtension.substr(fristPos+1).find_first_of("_");
    auto lastPos = fileNameWithoutExtension.find_last_of("_");
    evnN = std::stoul(fileNameWithoutExtension.substr(fristPos+1,secondPos));
    volxelN = std::stoul(fileNameWithoutExtension.substr(lastPos+1));
    return  {evnN, volxelN};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootTree1(TFile* f)
{
    TDirectoryFile* d = dynamic_cast<TDirectoryFile*> (f->Get("ntuple") );
    TTree* tPhys = dynamic_cast<TTree*> (d->Get("ntuple_1") );

    if( tPhys->GetEntries()  > 0)
        {
        double flagParticle;
        double flagParentID;
        double flagProcess;
        double x;
        double y;
        double z;
        double edep;
        double eventNumber;
        double volumeName;
        double copyNumber;
        double lastMetVoxelCopyNum;

        tPhys->SetBranchAddress("flagParticle", &flagParticle);
        tPhys->SetBranchAddress("flagParentID", &flagParentID);
        tPhys->SetBranchAddress("flagProcess", &flagProcess);
        tPhys->SetBranchAddress("x", &x);
        tPhys->SetBranchAddress("y", &y);
        tPhys->SetBranchAddress("z", &z);
        tPhys->SetBranchAddress("edep", &edep);
        tPhys->SetBranchAddress("eventNumber", &eventNumber);
        tPhys->SetBranchAddress("volumeName", &volumeName);
        tPhys->SetBranchAddress("copyNumber", &copyNumber);
        tPhys->SetBranchAddress("lastMetVoxelCopyNum", &lastMetVoxelCopyNum);

        unsigned int entryNumber =  tPhys->GetEntries() ;

        // Loop on all the "lines" (ie entry) of the ntuple
        for(unsigned int e=0; e<entryNumber; e++)
        {
            // Set all the variables to the values corresponding to the entry number
            tPhys->GetEntry(e);
            // Check if the process is an ionisation
            // Only ionisation should trigger the removal of a DNA molecule from the chemical step
            if(flagProcess == 13 // e-_DNAIonisation
                    || flagProcess == 113 // e-_DNAPTBIonisation
                    || flagProcess == 18 // proton_DNAIonisation
                    || flagProcess == 21 // hydrogen_DNAIonisation
                    || flagProcess == 24 // alpha_DNAIonisation
                    || flagProcess == 27 // alpha+_DNAIonisation
                    || flagProcess == 31 // helium_DNAIonisation
                    || flagProcess == 12 // e-_DNAExcitation
                    || flagProcess == 112 // e-_DNAPTBExcitation
                    || flagProcess == 15 // e-_DNAVibExcitation
                    || flagProcess == 17 // proton_DNAExcitation
                    || flagProcess == 20 // hydrogen_DNAExcitation
                    || flagProcess == 23 // alpha_DNAExcitation
                    || flagProcess == 26 // alpha+_DNAExcitation
                    || flagProcess == 30 // helium_DNAExcitation
                    ) {
                // Check the interaction happened in a dna molecule or its hydration shell
                if(volumeName == 1 // d1
                        || volumeName == 11 // p1
                        || volumeName == 2 // d2
                        || volumeName == 22 // p2
                        || volumeName == 7 // d1_w
                        || volumeName == 71 // p1_w
                        || volumeName == 8 // d2_w
                        || volumeName == 81 // p2_w
                        )
                {
                    // *************

                    // Retrieve the voxel copy number
                    double voxelCopyNumber  = lastMetVoxelCopyNum;
                    if (voxelCopyNumber >= 0 && voxelCopyNumber<fVoxels.size()){
                        // Chromosome, domain and firstNucleotideNum
                        const VoxelData& voxelData = fVoxels.at(size_t(voxelCopyNumber) );
                        int chromo = voxelData.fChromosome;
                        int domain = voxelData.fDomain;
                        ullint firstBpCN = voxelData.fFirstBpCopyNum;
                        ullint cpNumCorrected = firstBpCN+int(copyNumber);

                        // Get the event number
                        double eventNum = eventNumber;

                        // Determine the strand
                        double strand (-1);
                        if(volumeName==1
                                || volumeName==11
                                || volumeName==7
                                || volumeName==71
                                || volumeName==6 // ade
                                || volumeName==9 // ade
                                || volumeName==4 // gua
                                || volumeName==10) // gua
                            strand = 1;
                        else if(volumeName==2
                                || volumeName==22
                                || volumeName==8
                                || volumeName==81
                                || volumeName==5 // thy
                                || volumeName==12 // thy
                                || volumeName==3 // cyto
                                || volumeName==13) // cyto
                            strand = 2;
                        // Check if the chromo has already been registered
                        std::vector<ullint> newLine{ 
                            (ullint)eventNum,(ullint)strand,cpNumCorrected,
                            (ullint)volumeName,(ullint)flagProcess,(ullint)(edep*1000000)};
                            // *1000000 in edep is taken back after. This is done 
                            // because the number is "unsigned long long int" instead of "double".
                        auto itr = fphysTables.find(chromo);
                        if ( itr == fphysTables.end()) {
                            Table tbforThisChro{newLine};
                            fphysTables.insert({chromo,tbforThisChro});
                        } else (itr->second).push_back(newLine);
                    }
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::AnaPhysRootTree2(TFile* f)
{
    TDirectoryFile* d2 = dynamic_cast<TDirectoryFile*> (f->Get("ntuple") );
    TTree* tPhys2 = dynamic_cast<TTree*> (d2->Get("ntuple_3") );

    if( int(tPhys2->GetEntries() ) > 0)
    {
        double edep;
        double eventNumber;

        tPhys2->SetBranchAddress("edep", &edep);
        tPhys2->SetBranchAddress("eventNumber", &eventNumber);

        unsigned int entryNumber = int( tPhys2->GetEntries() );

        // Loop on all the "lines" (ie entry) of the ntuple
        for(unsigned int e=0; e<entryNumber; e++)
        {
            tPhys2->GetEntry(e);
            fEdepSumInNucleus += edep;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::SortPhysTableWithSelection()
{
    for(const auto [chrom, physTable] : fphysTables)
    {
        // Final table
        Table physTableWithSelection;

        std::map<ullint,std::map<ullint,std::map<ullint, ullint > > > energyMap;
        
        // Loop on all the lines of the table
        for(unsigned int line=0, eline=physTable.size(); line<eline; ++line)
        {
            ullint eventNum = physTable[line][0];
            ullint strand = physTable[line][1];
            ullint copyNumber = physTable[line][2];
            ullint volumeFlag = physTable[line][3];
            ullint processFlagr = physTable[line][4];
            ullint energy = physTable[line][5];

            // Cumulate the energy value
            energyMap[eventNum][strand][copyNumber] += energy;
        }

        int notDuplicatedLine = 0;

        // Loop on all the events
        std::map<ullint, std::map<ullint, std::map<ullint, ullint> > >::iterator iit = energyMap.begin();
        std::map<ullint, std::map<ullint, std::map<ullint, ullint> > >::iterator iite = energyMap.end();
        for(; iit!=iite;++iit)
        {
            ullint eventNum = iit->first;

            // Loop on all the strands
            std::map<ullint, std::map<ullint, ullint> >::iterator itt = iit->second.begin();
            std::map<ullint, std::map<ullint, ullint> >::iterator itte = iit->second.end();
            for(; itt!=itte;++itt)
            {
                ullint strand = itt->first;
                // Loop on all the copy numbers
                std::map<ullint, ullint>::iterator ittt = itt->second.begin();
                std::map<ullint, ullint>::iterator ittte = itt->second.end();
                for(; ittt!=ittte;++ittt)
                {
                    ullint copyNumber = ittt->first;
                    double currentE = double(ittt->second) / 1000000; // eV
                    // Energy condition(s) are set here
                    bool fill = false;                  
                     // Threshold condition
                     if(currentE < fThresholdEnergy)
                          fill=false;
                     else
                          fill=true;

                    if(fill)
                    {
                        // Add a line
                        physTableWithSelection.push_back(std::vector<ullint>());
                        // Fill the line
                        physTableWithSelection[notDuplicatedLine].push_back(eventNum);
                        physTableWithSelection[notDuplicatedLine].push_back(strand);
                        physTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                        physTableWithSelection[notDuplicatedLine].push_back(0);
                        physTableWithSelection[notDuplicatedLine].push_back(0);
                        physTableWithSelection[notDuplicatedLine].push_back((ullint)(currentE*1000000) );
                        physTableWithSelection[notDuplicatedLine].push_back(0);
                        ++notDuplicatedLine;
                    }
                }
            }
        }
        // *******************************************
        // Print the "physTableWithSelection" table for the current chromosome
        // *******************************************
        std::cout << "### Phys SB for chromosome "<<chrom<<" : " << physTableWithSelection.size() << " ###" << std::endl;
        if (physTableWithSelection.size()>0) fphysSlectedTables.insert({chrom,physTableWithSelection});
    }
    fphysTables.clear();//Free memory
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::SortChemTableWithSelection()
{
    for (const auto& [chromo,chemTable] : fchemTables)
    {
      	// Final table
        Table chemTableWithSelection;
        
        int notDuplicatedLine = 0;

        for(unsigned int line=0, eline=chemTable.size(); line<eline; ++line)
        {
            ullint eventNum = chemTable[line][0];
            ullint strand = chemTable[line][1];
            ullint copyNumber = chemTable[line][2];
            ullint base = chemTable[line][3];
            ullint time = chemTable[line][4];
        
     		// Random number between 0 and <1
			if (base==1) {
                // Add a line
                chemTableWithSelection.push_back(std::vector<ullint>());
                // Fill the line
                chemTableWithSelection[notDuplicatedLine].push_back(eventNum);
                chemTableWithSelection[notDuplicatedLine].push_back(strand);
                chemTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                chemTableWithSelection[notDuplicatedLine].push_back(base);
                chemTableWithSelection[notDuplicatedLine].push_back(time);
                chemTableWithSelection[notDuplicatedLine].push_back(0);
                chemTableWithSelection[notDuplicatedLine].push_back(1);              
                ++notDuplicatedLine;
			}
			else { 
			    //double r = double(std::rand() ) / RAND_MAX;
                double r = gRandomGen.Rndm(); 
                //std::cout<<r<<std::endl;
			    if(r <= fProbabilityForIndirectSB){
                    // Add a line
                    chemTableWithSelection.push_back(std::vector<ullint>());
                    // Fill the line
                    chemTableWithSelection[notDuplicatedLine].push_back(eventNum);
                    chemTableWithSelection[notDuplicatedLine].push_back(strand);
                    chemTableWithSelection[notDuplicatedLine].push_back(copyNumber);
                    chemTableWithSelection[notDuplicatedLine].push_back(base);
                    chemTableWithSelection[notDuplicatedLine].push_back(time);
                    chemTableWithSelection[notDuplicatedLine].push_back(0);
                    chemTableWithSelection[notDuplicatedLine].push_back(1);                          
                    ++notDuplicatedLine;
                }
	        }
      	}   	
        if (chemTableWithSelection.size() > 0) fchemSlectedTables.insert({chromo,chemTableWithSelection});
	}
    fchemTables.clear();//free memory
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::MergeDamageFromPhysChem()
{
    for(int  chromo=0; chromo<46; chromo++)
    {
        // *****************************
        // Add one table after the other
        // *****************************

        // MergedTable to be built
        Table mergedTable;

        auto itr = fphysSlectedTables.find(chromo);
        if(itr != fphysSlectedTables.end())
        {
            Table physTable= itr->second;
            for(auto const& vec : physTable) mergedTable.push_back(vec);
        }
        itr = fchemSlectedTables.find(chromo);
        if(itr != fchemSlectedTables.end())
        {
            Table chemTable = itr->second;
            for(auto const& vec : chemTable) mergedTable.push_back(vec);
        }

        // *******************************************
        // Sort the merged table to put the event in the correct order
        // *******************************************
        Trier tri;
        Table::iterator it = mergedTable.begin();
        Table::iterator ite = mergedTable.end();
        std::sort(it, ite, tri);

        // *******************************************
        // Delete duplicate SB
        // *******************************************

        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>> removeDuplicateValueMap;

        // Put all the values of the mergedTable in the map created just above.
        // Loop on all the lines of the mergedTable
        for(unsigned int line=0; line<mergedTable.size(); ++line)
        {
            ullint eventNum = mergedTable[line][0];
            ullint strand = mergedTable[line][1];
            ullint copyNumber = mergedTable[line][2];
            ullint base = mergedTable[line][3];
            std::vector<ullint> lineV;
            // If more elements are presents, add them here
            int lineSize = mergedTable[line].size();
            if(lineSize > 4)
            {
                for(int i=4; i<lineSize; i++)
                {
                    lineV.push_back(mergedTable[line][i]);
                }
                removeDuplicateValueMap[eventNum][strand][copyNumber][base]=lineV;
            }
            else
            {
                lineV.push_back(0);
                removeDuplicateValueMap[eventNum][strand][copyNumber][base]=lineV;
            }
        }

        // *******************************************
        // Create the "mergedTableWithoutDuplicatedSB" table
        // *******************************************

        // At this point, we have created a map named "removeDuplicateValueMap" that organized all the mergedTable content
        // AND that does not contain any duplicate because of the override caracteristic of a map.
        // Indeed, doing map[2] = "hello" followed by map[2]  = "bye" will put map[2] value to "bye" because it overrided the first "hello".
        // This is a cheap way to remove duplicates by overriding them.

        // The next part is dedicated to the creation of the "mergedTableWithoutDuplicatedSB" table from the "removeDuplicateValueMap" map.
        // Only the event, strand and copynumber and  will be put in this final table.
        Table mergedTableWithoutDuplicatedSB;
        Table mergedTableWithoutDuplicatedSBandbases;
        int notDuplicatedLine = 0;
        int notDuplicatedLine2= 0;
        // Loop on all the events
        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>>::iterator 
            iit = removeDuplicateValueMap.begin();
        std::map<ullint,std::map<ullint,std::map<ullint,std::map<ullint,std::vector<ullint>>>>>::iterator 
            iite = removeDuplicateValueMap.end();
        for(; iit!=iite;++iit)
        {
            ullint eventNum = iit->first;
            // Loop on all the strands
            std::map<ullint, std::map<ullint, std::map<ullint, std::vector<ullint > > > >::iterator 
                itt = iit->second.begin();
            std::map<ullint, std::map<ullint, std::map<ullint, std::vector<ullint > > > >::iterator 
                itte = iit->second.end();
            for(; itt!=itte;++itt)
            {
                ullint strand = itt->first;
                // Loop on all the copy numbers
                std::map<ullint, std::map<ullint, std::vector<ullint > > >::iterator ittt = itt->second.begin();
                std::map<ullint, std::map<ullint, std::vector<ullint > > >::iterator ittte = itt->second.end();
                for(; ittt!=ittte;++ittt)
                {
                    ullint copyNumber = ittt->first;                    
		            // Loop on all the base flag
                    std::map<ullint, std::vector<ullint > >::iterator itttt = ittt->second.begin();
                    std::map<ullint, std::vector<ullint > >::iterator itttte = ittt->second.end();
                    for(; itttt!=itttte;++itttt)
                    {
                        ullint base = itttt->first;
                        // Fill the table
                        if (strand>0 && strand<3)
                        {
                            // Add a line
                            mergedTableWithoutDuplicatedSB.push_back(std::vector<ullint>());                           
                            // Fill the line
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(eventNum);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(strand);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(copyNumber);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(base);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(
                                removeDuplicateValueMap[eventNum][strand][copyNumber][base][0]);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(
                                removeDuplicateValueMap[eventNum][strand][copyNumber][base][1]);
                            mergedTableWithoutDuplicatedSB[notDuplicatedLine].push_back(
                                removeDuplicateValueMap[eventNum][strand][copyNumber][base][2]);
                            ++notDuplicatedLine;                           
                            if (base==0)
                            {
                                // Add a line
                                mergedTableWithoutDuplicatedSBandbases.push_back(std::vector<ullint>());                              
                                // Fill the line
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(eventNum);
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(strand);
                                mergedTableWithoutDuplicatedSBandbases[notDuplicatedLine2].push_back(copyNumber);                               
                                ++notDuplicatedLine2;
                            }
                        }
                    }
                }
            }
        }

        // *******************************************
        // Print the "mergedTableWithoutDuplicatedSB" table for the current chromosome
        // *******************************************
        if (mergedTableWithoutDuplicatedSB.size() > 0) {
            //PrintTable(fMergeFolder + "/chromo_"+std::to_string(chromo)+".dat", 
            //mergedTableWithoutDuplicatedSB, "eventNum, strand, copyNumber, isbase, 
            //time(ns*1000000), edep, phy:0 chem:1");
            fMergedTables.insert({chromo,mergedTableWithoutDuplicatedSB});
        }
        mergedTable.clear();
        mergedTableWithoutDuplicatedSBandbases.clear();
        mergedTableWithoutDuplicatedSB.clear();
    }
    //free memory:
    fphysSlectedTables.clear();
    fchemSlectedTables.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ScanDamage::ReadCellandVoxelDefFilePaths()
{
    fs::path thisP = fs::current_path();
    for (const auto entry : fs::directory_iterator(thisP)){
        if (entry.path().filename() == "imp.info") {
            std::ifstream file(entry.path().c_str());
            if(!file.good() ){
                std::cerr<<"**** Fatal Error *****"<<std::endl;
                std::cerr<<"ScanDamage::ReadCellandVoxelDefFilePaths(): File corupted: "
                <<entry.path()<<std::endl;
                std::cerr<<"*************** *****"<<std::endl;
                exit(EXIT_FAILURE);
            }

            std::string line;
            while(std::getline(file, line) ){
                std::istringstream iss(line);
                std::string flag;
                iss >> flag;
                if ( flag == "_geovolxelpath") {
                    std::string voxname;
                    iss >> voxname;
                    fVoxelDefFilesList.insert(voxname);
                }
                if ( flag == "_geocellpath") {
                    std::string cellpname;
                    iss >> cellpname;
                    fCellDefFilePath = cellpname;
                }
                if ( flag == "_numberOfBasepairs") {
                    iss >> fTotalNbBpPlacedInGeo;
                }
                if ( flag == "_numberOfHistones") {
                    iss >> fTotalNbHistonePlacedInGeo;
                }
                if ( flag == "_nucleusVolume") {
                    iss >> fNucleusVolume;
                }
                if ( flag == "_nucleusMassDensity") {
                    iss >> fNucleusMassDensity;
                }
                if ( flag == "_nucleusMass") {
                    iss >> fNucleusMass;
                }
            }
            file.close();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....