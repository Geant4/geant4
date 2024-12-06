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
/// \file Analysis.cc
/// \brief Implementation of the Analysis class
/// \file Analysis.cc
/// \brief Implementation of the Analysis class

#include "Analysis.hh"
#include "DetectorConstruction.hh"
#include "G4AutoDelete.hh"
#include "G4Filesystem.hh"
#include "G4DNAMolecule.hh"


#include <sstream>
#include <fstream>
#ifdef USE_MPI
#include "G4MPImanager.hh"
#endif

G4ThreadLocal Analysis* the_analysis = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Analysis* Analysis::GetAnalysis()
{
    if (!the_analysis) {
        the_analysis = new Analysis();
        G4AutoDelete::Register(the_analysis);
    }

    return the_analysis;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::OpenFile(const G4String outFolder)
{
    G4AnalysisManager* anManager = G4AnalysisManager::Instance();
    anManager->SetDefaultFileType("root");
    G4String fullFileName = fFileName;
    if (gRunMode == RunningMode::Chem) {
        G4String slash = "";
#if defined(_WIN32) || defined(WIN32)
        slash= "\\";
#else
        G4String slashu = "/";
        slash = slashu;
#endif
        fullFileName = outFolder+slash+fFileName;
    }

    // open output file
    G4bool fileOpen = anManager->OpenFile(fullFileName.c_str());
    if (!fileOpen) {
        G4cout << "\n---> HistoManager::book(): cannot open " << fFileName
            << G4endl;
        return;
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::Save()
{
    G4AnalysisManager* anManager = G4AnalysisManager::Instance();
    anManager->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::Close(G4bool reset)
{
    G4AnalysisManager* anManager = G4AnalysisManager::Instance();
    anManager->CloseFile(reset);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::Book()
{   
    G4AnalysisManager* anManager = G4AnalysisManager::Instance();
    anManager->SetVerboseLevel(1);
    if (anManager->GetFirstNtupleId() != 1) anManager->SetFirstNtupleId(1);
    
    if (gRunMode == RunningMode::Phys) {     
#ifdef G4MULTITHREADED
        // MT ntuple merging
        anManager->SetNtupleMerging(true);//for future use MT+MPI
#endif
        anManager->SetNtupleDirectoryName("ntuple");
        // create ntuple
        anManager->CreateNtuple("ntuple_1","physical_stage");
        anManager->CreateNtupleIColumn("flagParticle");
        anManager->CreateNtupleIColumn("flagParentID");
        anManager->CreateNtupleIColumn("flagProcess");
        anManager->CreateNtupleDColumn("x");
        anManager->CreateNtupleDColumn("y");
        anManager->CreateNtupleDColumn("z");
        anManager->CreateNtupleDColumn("edep");
        anManager->CreateNtupleIColumn("eventNumber");
        anManager->CreateNtupleIColumn("volumeName");
        anManager->CreateNtupleIColumn("copyNumber");
        anManager->CreateNtupleIColumn("lastMetVoxelCopyNum");
        anManager->FinishNtuple(1);

        // For total edep
        anManager->CreateNtuple("ntuple_3","total_edep");
        anManager->CreateNtupleIColumn("eventNumber");
        anManager->CreateNtupleDColumn("edep");
        anManager->FinishNtuple(2);
    }

    if (gRunMode == RunningMode::Chem) {
        // Create directories
        const G4String directoryName = "ntuple";
        if (anManager->GetNtupleDirectoryName() != directoryName) {
            anManager->SetNtupleDirectoryName(directoryName);
        }
        
        // DBScan

        anManager->CreateNtuple("ntuple_2","DB_chemical_stage");
        anManager->CreateNtupleIColumn(1,"strand");
        anManager->CreateNtupleIColumn(1,"copyNumber");
        anManager->CreateNtupleDColumn(1,"xp");
        anManager->CreateNtupleDColumn(1,"yp");
        anManager->CreateNtupleDColumn(1,"zp");
        anManager->CreateNtupleDColumn(1,"time");
        anManager->CreateNtupleIColumn(1,"base");
        anManager->FinishNtuple(1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AnalysisManager* Analysis::GetAnalysisManager()
{
    return G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::AddInfoForChemGeo(InfoForChemGeo b)
{
    fInfoForChemGeoVector.push_back(b);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::ClearVector()
{
    fInfoForChemGeoVector.clear();
    fInfoInPhysStageVector.clear();
    fOutputFiles.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::AddInfoInPhysStage(InfoInPhysStage b)
{
    fInfoInPhysStageVector.push_back(b);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::UpdateChemInputDataAndFillNtuple()
{
    //if (fInfoForChemGeoVector.size() == 0) return;

    // We will loop on all the "element" of the vector and create several files:
    // - New file for each event/voxel couple
    // Create a map to define all the output file

    for (auto const& binfo : fInfoForChemGeoVector) {
        if( binfo.fVolume == 161  // voxelStraight
            || binfo.fVolume == 162 // voxelRight
            || binfo.fVolume == 163 // voxelLeft
            || binfo.fVolume == 164 // voxelUp
            || binfo.fVolume == 165 // voxelDown
            || binfo.fVolume == 261 // voxelStraight2
            || binfo.fVolume == 262 // voxelRight2
            || binfo.fVolume == 263 // voxelLeft2
            || binfo.fVolume == 264 // voxelUp2
            || binfo.fVolume == 265) // voxelDown2
        {
            // We are in a voxel

            std::string voxelName("");

            if(binfo.fVolume==161) voxelName = "VoxelStraight";
            else if(binfo.fVolume==162) voxelName = "VoxelRight";
            else if(binfo.fVolume==163) voxelName = "VoxelLeft";
            else if(binfo.fVolume==164) voxelName = "VoxelUp";
            else if(binfo.fVolume==165) voxelName = "VoxelDown";
            else if(binfo.fVolume==261) voxelName = "VoxelStraight2";
            else if(binfo.fVolume==262) voxelName = "VoxelRight2";
            else if(binfo.fVolume==263) voxelName = "VoxelLeft2";
            else if(binfo.fVolume==264) voxelName = "VoxelUp2";
            else if(binfo.fVolume==265) voxelName = "VoxelDown2";
            
            // Get the event number
            int eventNum = int(binfo.fEventNumber);
            // Here we have all the information for one ntuple line
            // We want to know if we have to create a new output file (new eventNumber/voxel couple)
            // or if we already have one.

            // Check if the event has already been registered
            if(fOutputFiles.find(eventNum)==fOutputFiles.end() )
            {
                // If not then create the event and voxel case
                fOutputFiles[eventNum][binfo.fVolumeCopyNumber] = 
                CreateChemInputFile(eventNum, int(binfo.fVolumeCopyNumber), voxelName);
            }
            // If the event has been registered then we need to check that the current voxel has also been
            // registered.
            else
            {
                if(fOutputFiles[eventNum].find(binfo.fVolumeCopyNumber)==fOutputFiles[eventNum].end())
                {
                    // We register the event-voxel couple
                    fOutputFiles[eventNum][binfo.fVolumeCopyNumber] = 
                    CreateChemInputFile(eventNum, int(binfo.fVolumeCopyNumber), voxelName);
                }
            }

            // Write in the file the information needed by the chemistry simulation to build
            // the input water molecule or solvated electron.
            UpdatingChemInputFile(binfo);
        }
    }

    if (fInfoInPhysStageVector.size() == 0) return;
    for (auto const& binfo : fInfoInPhysStageVector) {
        // Check if the process is an ionisation
        // Only ionisation should trigger the removal of a DNA molecule from the chemical step
            if(binfo.fFlagProcess == 13
            || binfo.fFlagProcess == 113
            || binfo.fFlagProcess == 18
            || binfo.fFlagProcess == 21
            || binfo.fFlagProcess == 24
            || binfo.fFlagProcess == 27
            || binfo.fFlagProcess == 31) {
                // Check the interaction happened in a dna molecule or its hydration shell
                if( binfo.fVolumeName == 1 // d1
                    || binfo.fVolumeName == 11 // p1
                    || binfo.fVolumeName == 2 // d2
                    || binfo.fVolumeName == 22 // p2
                    || binfo.fVolumeName == 3 // cyto
                    || binfo.fVolumeName == 4 // gua
                    || binfo.fVolumeName == 5 // thy
                    || binfo.fVolumeName == 6 // ade
                    || binfo.fVolumeName == 7 // d1_w
                    || binfo.fVolumeName == 71 // p1_w
                    || binfo.fVolumeName == 8 // d2_w
                    || binfo.fVolumeName == 81 // p2_w
                    || binfo.fVolumeName == 9 // ade_w
                    || binfo.fVolumeName == 10 // gua_w
                    || binfo.fVolumeName == 13 // cyto_w
                    || binfo.fVolumeName == 12) // thy_w
                {
                    // Retrieve the voxel copy number
                    double voxelCopyNumber  = binfo.fLastMetVoxelCopyNum;
                    // Get the event number
                    double eventNum = binfo.fEventNumber;

                    // Check if the event has already been registered
                    if(fOutputFiles.find(eventNum) != fOutputFiles.end() )
                    {
                        // Check if the volume number has been registered
                        if(fOutputFiles.at(eventNum).find(voxelCopyNumber) 
                            != fOutputFiles.at(eventNum).end() )
                        {
                            // If we are here then the event and volume couple has a already generated 
                            //file in which we should add a dna molecule to be removed
                            UpdatingChemInputFile(binfo);
                        }
                    }
                }
        }

        // fill ntuple
        auto analysisManager =G4AnalysisManager::Instance();
        analysisManager->FillNtupleIColumn(1, 0, binfo.fFlagParticle);
        analysisManager->FillNtupleIColumn(1, 1, binfo.fFlagParentID);
        analysisManager->FillNtupleIColumn(1, 2, binfo.fFlagProcess);
        analysisManager->FillNtupleDColumn(1, 3, binfo.fX);
        analysisManager->FillNtupleDColumn(1, 4, binfo.fY);
        analysisManager->FillNtupleDColumn(1, 5, binfo.fZ);
        analysisManager->FillNtupleDColumn(1, 6, binfo.fEdep);
        analysisManager->FillNtupleIColumn(1, 7, binfo.fEventNumber);
        analysisManager->FillNtupleIColumn(1, 8, binfo.fVolumeName);
        analysisManager->FillNtupleIColumn(1, 9, binfo.fCopyNumber);
        analysisManager->FillNtupleIColumn(1, 10, binfo.fLastMetVoxelCopyNum);
        analysisManager->AddNtupleRow(1);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Analysis::CreateChemInputFile(G4int eventNum, G4int volumeCopyNumber, const G4String &voxelName)
{
    std::stringstream sstream;
    sstream<<"./"<<fChemInputFolderName<<"/event_"<<eventNum<<"_voxel_"<<volumeCopyNumber<<".dat";
    std::ofstream oFile;
    oFile.open(sstream.str().c_str() );

    oFile<<"_eventNum"<<"\t\t"<<eventNum<<std::endl;
    oFile<<"_voxelType"<<"\t\t"<<voxelName<<std::endl;
    oFile<<"_voxelCopyNumber"<<"\t\t"<<volumeCopyNumber<<std::endl;
    oFile<<std::endl;

    oFile
            <<"# Chemistry input informations\n"
            <<"# "<<"_input, type, state, electronicLevel, x, y, z, parentTrackID"<<std::endl;
    oFile<<"# type=1 -> water molecule && type=2 -> solvated electron"<<std::endl;
    oFile<<std::endl;

    oFile.close();

    return sstream.str().c_str();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::UpdatingChemInputFile(InfoForChemGeo b)
{
    std::ofstream out;
    out.open( fOutputFiles[b.fEventNumber][b.fVolumeCopyNumber].c_str(), std::ios::app); // Open the file and go at the end

    if (b.fType==1)
    {
    out<<"_input"<<"\t"
                                            <<b.fType<<"\t\t"
                                            <<b.fState<<"\t\t"
                                        <<4-b.fElectronicLevel<<"\t\t"
                                        <<b.fRelX<<"\t\t"
                                        <<b.fRelY<<"\t\t"
                                        <<b.fRelZ<<"\t\t"
                                    <<b.fParentTrackID<<"\t\t"
                                <<"\n";
    } 
    if (b.fType==2)
    {
    out<<"_input"<<"\t"
                                            <<b.fType<<"\t\t"
                                            <<b.fState<<"\t\t"
                                        <<b.fElectronicLevel<<"\t\t"
// If we are here then the event and volume couple has a already 
//generated file in which we should add a dna molecule to be removed
                                        <<b.fRelX<<"\t\t"
                                        <<b.fRelY<<"\t\t"
                                        <<b.fRelZ<<"\t\t"
                                    <<b.fParentTrackID<<"\t\t"
                                <<"\n";
    }                              
    out.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::UpdatingChemInputFile(InfoInPhysStage b)
{
    G4String name;
    G4double volumeName =  b.fVolumeName;
    if(volumeName == 1 || volumeName == 2 || volumeName == 7 || 
            volumeName == 8) name =  G4Deoxyribose::Definition()->GetName();
    else if(volumeName == 11 || volumeName == 22 || 
            volumeName == 71 || volumeName == 81) name =  G4Phosphate::Definition()->GetName();
    else if(volumeName == 6 || volumeName == 9) name =  G4Adenine::Definition()->GetName();
    else if(volumeName == 4 || volumeName == 10) name =  G4Guanine::Definition()->GetName();
    else if(volumeName == 5 || volumeName == 12) name =  G4Thymine::Definition()->GetName();
    else if(volumeName == 3 || volumeName == 13) name =  G4Cytosine::Definition()->GetName();
    else
    {
        G4ExceptionDescription msg;
        msg <<"Volume number "<<volumeName<<" not registered.";
        G4Exception("Analysis::UpdatingChemInputFile", 
                    "", FatalException, msg);
    }

    G4double strand (-1);

    // Determine the strand
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
    std::ofstream out;
    out.open(fOutputFiles[b.fEventNumber][b.fLastMetVoxelCopyNum].c_str(), std::ios::app); // Open the file and go at the end
    out<<"_remove"<<"\t"<<name<<"\t\t"<<b.fCopyNumber<<"\t\t"<<strand<<"\t\t"<<std::endl; 
    out.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::WritePhysGeo()
{
    //Create a file containing geopath and Physlist in "FS build directory" for chemStage and anlysis
    std::string fname = "imp.info";
    std::ofstream fout(fname);
    fout<<"====> Auto_generated file. Do not delete me !!!\n";
    fout<<"====> The file conveys some information for chem_geo and Analysis modules!!!\n";
    fout<<"_geocellpath "<<fCellDefFilePath<<"\n";
    for (const auto & entry : fVoxelDefFilesList) {
        fout<<"_geovolxelpath "<<std::string(entry)<<"\n";
    }
    fout<<"_numberOfBasepairs "<<fTotalNbBpPlacedInGeo<<"\n";
    fout<<"_numberOfHistones "<<fTotalNbHistonePlacedInGeo<<"\n";
    fout<<"_nucleusVolume "<<fNucleusVolume<<"\n"; // m3
    fout<<"_nucleusMassDensity "<<fNucleusMassDensity<<"\n"; // kg/m3
    fout<<"_nucleusMass "<<fNucleusVolume*fNucleusMassDensity; //kg
    fout.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Analysis::DefineCommands()
{
    fMessenger = std::make_unique<G4GenericMessenger>(this,
                             "/dsbandrepair/output/",
                             "cmd controld");
    auto & fChemOutFolderCmd = fMessenger->DeclareProperty ("folderForChemOut",
                                                 fChemOutFolderName);
    fChemOutFolderCmd.SetParameterName("outFolderForChem",true);
    fChemOutFolderCmd.SetDefaultValue("chem_output");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Analysis::CheckAndCreateNewFolderInChemStage()
{
    G4fs::file_status myFs = G4fs::file_status{};
    auto folderPath = G4fs::path{fChemOutFolderName.c_str()};
    auto isExist = G4fs::status_known(myFs) ? G4fs::exists(myFs) : G4fs::exists(folderPath);
    if (isExist && !G4fs::is_empty(folderPath)) {
        G4ExceptionDescription msg;
        msg <<"==>> Chem output folder "<<fChemOutFolderName
            <<" is already existing and not empty!!!\n";
        msg<<"==>> Please delete or rename it, or use other name for  "
        <<"Chem output folder in macro file!!!\n";
        G4Exception("Analysis::CheckAndCreateNewFolderInChemStage()","",FatalException,msg);
    }
    G4fs::create_directory(folderPath);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Analysis::CheckAndCreateNewFolderInPhysStage()
{
    G4fs::file_status myFs = G4fs::file_status{};
    const G4fs::path folderPath{fPhysOutFolderName.c_str()};
    auto isExist = G4fs::status_known(myFs) ? G4fs::exists(myFs) : G4fs::exists(folderPath);
    if (isExist && !G4fs::is_empty(folderPath)) {
        G4fs::remove_all(folderPath);
    }
    G4fs::create_directory(folderPath);

    const G4fs::path folderPathPC{fChemInputFolderName.c_str()};
    isExist = G4fs::status_known(myFs) ? G4fs::exists(myFs) : G4fs::exists(folderPathPC);
    if (isExist && !G4fs::is_empty(folderPathPC)) {
        G4fs::remove_all(folderPathPC);
    }
    G4fs::create_directory(folderPathPC);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....