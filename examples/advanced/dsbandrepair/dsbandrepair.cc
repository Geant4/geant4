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
//  author: Le Tuan Anh, 20/10/2023
/// \file main.cc
/// \brief Main program of the dsbandrepair

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIExecutive.hh"

#include "G4RunManagerFactory.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4Timer.hh"
#include "G4ExceptionSeverity.hh"
#include "G4DNAChemistryManager.hh"
#include "G4VisExecutive.hh"
#include "G4Filesystem.hh"

#include "PhysActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "InformationKeeper.hh"

#include "ChemActionInitialization.hh"
#include "ChemPhysicsList.hh"
#include "ChemNtupleManager.hh"

#ifdef USE_MPI
#include "G4MPImanager.hh"
#include "G4MPIsession.hh"
#include "G4MPIextraWorker.hh"
#endif
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CheckingSomeFilesAndFolders();
G4String ExtractChemListNameFromMacroFile(G4String);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv)
{
#ifdef USE_MPI
        G4MPImanager* g4MPI = new G4MPImanager(argc, argv, 0);
        g4MPI->SetVerbose(1);
        G4MPIsession* session = g4MPI-> GetMPIsession();
        G4String prompt = "[40;01;33m";
        prompt += "G4MPI";
        prompt += "[40;31m(%s)[40;36m[%/][00;30m:";
        session-> SetPrompt(prompt);
#else 
        G4UIExecutive* ui = nullptr;
        if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }
#endif // USE_MPI
    if (argc < 2) {
        G4cerr<<"====>> Wrong input. To run Physgeo, type : ./dsbandrepair macrofile\n"
            <<"To run Chem_geo, type : ./dsbandrepair macrofile chem"<<G4endl;
#ifdef USE_MPI
        delete g4MPI;
#endif // USE_MPI
        return EXIT_FAILURE;
    }
    G4String macrofileName = argv[1];
    if (argc > 2) {
        const G4String rmode = argv[2];
        if (rmode == "phys") gRunMode = RunningMode::Phys;
        if (rmode == "chem") gRunMode = RunningMode::Chem;
    } 
    // Choose the Random engine
    time_t timeStart;
    time(&timeStart);
    unsigned long seed = timeStart;
#ifdef USE_MPI
    // Le Tuan Anh: add rankID to get different seeds for multi parallel processes
    seed += 1987*g4MPI->GetRank(); 
#endif // USE_MPI
    G4cout<<"Initial Seed for random engine: "<<seed<<G4endl;
    CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
    CLHEP::HepRandom::setTheSeed(seed);
    if (gRunMode == RunningMode::Phys) {
#ifdef USE_MPI
        if (g4MPI->GetRank() == 0 ){
            CheckingSomeFilesAndFolders();
        }
#else   
        CheckingSomeFilesAndFolders();
#endif // USE_MPI
        G4Random::setTheSeed(seed);
        auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
        G4int threadNumber= 1;
        runManager-> SetNumberOfThreads(threadNumber);
#endif  // G4MULTITHREADED
        
        DetectorConstruction* detector = new DetectorConstruction(1.,0,false);
        runManager->SetUserInitialization(detector);

        PhysicsList* physList = new PhysicsList;
        runManager->SetUserInitialization(physList);

        PhysActionInitialization* actionIni = new PhysActionInitialization();
        runManager->SetUserInitialization(actionIni);
#ifdef USE_MPI 
        // extra worker (for collecting ntuple data)
        if ( g4MPI->IsExtraWorker() ) {
            G4cout << "Set extra worker" << G4endl;
            G4UserRunAction* runAction = const_cast<G4UserRunAction*>(runManager->GetUserRunAction());
            g4MPI->SetExtraWorker(new G4MPIextraWorker(runAction));
        }
        session-> SessionStart();
        

        if (g4MPI->GetRank() == 0 ){
            InformationKeeper::Instance()->WritePhysGeo();
        }
        delete g4MPI;
#else 
        // Get the pointer to the User Interface manager
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        // Process macro or start UI session
        if ( ! ui ) {
            // batch mode
            G4String command = "/control/execute ";
            UImanager->ApplyCommand(command+macrofileName);
        }
        InformationKeeper::Instance()->WritePhysGeo();
#endif // USE_MPI
        delete runManager;
    } else if (gRunMode == RunningMode::Chem) {
        std::string inputFileorFolder = "chem_input";
        if (argc == 4) inputFileorFolder = argv[3];
        G4fs::path p{inputFileorFolder};
        G4String outputFileName = "test";
        std::vector<G4String> totalNumberofFilesVector, numberOfFilesTobeProcessedVector;
        if (G4fs::is_directory(p)) {
            for (const auto& entry : G4fs::directory_iterator(p)) {
                if (entry.path().extension() == ".dat") {
                    totalNumberofFilesVector.push_back(entry.path().string());
                }
            }
            std::sort(totalNumberofFilesVector.begin(),totalNumberofFilesVector.end());
#ifdef USE_MPI
            G4int numberofRanks = g4MPI->GetActiveSize();
            size_t filesTobeProcessedSlave = (size_t)(
                std::floor(G4double(totalNumberofFilesVector.size())/G4double(numberofRanks)));
            // note: should not use "std::ceil"
            size_t filesTobeProcessedMaster = 
                totalNumberofFilesVector.size() - (numberofRanks-1)*filesTobeProcessedSlave;
            if (g4MPI->IsMaster()) {
                for (size_t ii=0; ii< filesTobeProcessedMaster; ii++) {
                    numberOfFilesTobeProcessedVector.push_back(totalNumberofFilesVector.at(ii));
                }
            } else {
                for (size_t ii=0; ii< filesTobeProcessedSlave; ii++) {
                    auto rankID = g4MPI->GetRank();
                    size_t kk = filesTobeProcessedMaster + (rankID-1)*filesTobeProcessedSlave + ii;
                    numberOfFilesTobeProcessedVector.push_back(totalNumberofFilesVector.at(kk));
                }
            }
            G4cout<<"-----> "<<numberOfFilesTobeProcessedVector.size()
                    <<" files will be processed on rank #"<<g4MPI->GetRank()<<G4endl;
#else 
            numberOfFilesTobeProcessedVector = totalNumberofFilesVector;
#endif
            if (totalNumberofFilesVector.size() == 0) {
                G4cout<<"===>> There is no files found in "<<inputFileorFolder
                <<". You have to run Phys_geo first!!!"<<G4endl;
                //return EXIT_SUCCESS;
            } else {
                G4cout<<"===>> Total files found in "<<inputFileorFolder
                <<" : "<<totalNumberofFilesVector.size()<<G4endl;
            }
        } else if (G4fs::is_regular_file(p)) {
            numberOfFilesTobeProcessedVector.push_back(inputFileorFolder);
            if (p.has_stem()) {
                outputFileName = p.stem().string();
            } else outputFileName = inputFileorFolder;
        } 
        else G4cout<<"===>>dsbandrepair: "<<p.string()<<" is Not Directory or file !!!"<<G4endl;
        
        //------------------------------------------
        // Initialization classes
        //------------------------------------------

        auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
        
        G4DNAChemistryManager::Instance()->SetChemistryActivation(true);

        ChemNtupleManager ntupleManager;

        DetectorConstruction* detector = new DetectorConstruction();
        G4String firstFileForInit="";
        if (numberOfFilesTobeProcessedVector.size()>0) {
            firstFileForInit=numberOfFilesTobeProcessedVector.at(0);
            detector->ParseGeoFileForChemMode(firstFileForInit); // read to build voxel
        }
        runManager->SetUserInitialization(detector);
        G4String chemListName = ExtractChemListNameFromMacroFile(macrofileName);
        ChemPhysicsList* physList;
        if (chemListName != "") {
            physList = new ChemPhysicsList(chemListName);
        } else physList = new ChemPhysicsList();
        runManager->SetUserInitialization(physList);

        ChemActionInitialization* actionIni = new ChemActionInitialization(&ntupleManager,physList);
        runManager->SetUserInitialization(actionIni);
        
        //get the pointer to the User Interface manager
        G4UImanager* UI = G4UImanager::GetUIpointer();
        G4String command = "/control/execute ";
        UI->ApplyCommand(command+macrofileName);
        runManager->Initialize();
        if (numberOfFilesTobeProcessedVector.size()>0) {
            size_t nprocessedfiles{0}, ncounts{1};
            if (numberOfFilesTobeProcessedVector.size()>=100) ncounts=10;
            if (numberOfFilesTobeProcessedVector.size()>=10000) ncounts=100;
            if (numberOfFilesTobeProcessedVector.size()>=1000000) ncounts=500;
            for (auto const &fileInput : numberOfFilesTobeProcessedVector) {
                G4fs::path aP{std::string(fileInput)};
                if (aP.has_stem()) {
                    outputFileName = aP.stem().string();
                } else outputFileName = fileInput;
                ntupleManager.SetFileName(outputFileName);
                if (fileInput != firstFileForInit) detector->ParseGeoFileForChemMode(fileInput);
                detector->InsertMoleculeInWorld();
                UI->ApplyCommand("/run/beamOn 1");
                nprocessedfiles++;
                if (nprocessedfiles == 1 || 
                    nprocessedfiles == numberOfFilesTobeProcessedVector.size() ||
                    0 == (nprocessedfiles % ncounts)) {
                    G4cout<<"=====> Processed file: "<<nprocessedfiles<<"-th/("
                          <<numberOfFilesTobeProcessedVector.size()<<" files)"
#ifdef USE_MPI
                          <<" in rank #"<<g4MPI->GetRank()
#endif
                          <<"!!!"<<G4endl;
                }
            }   
        } else {
            UI->ApplyCommand("/run/beamOn 1");
        }
#ifdef USE_MPI
        delete g4MPI;
#endif // USE_MPI
        delete runManager;
    } else {
        G4cout<<"Undefined Running Mode; dsbansrepair will quit now. See you!\n";
    }

    return EXIT_SUCCESS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CheckingSomeFilesAndFolders()
{
    const G4fs::path curPath{G4fs::current_path()};
    // create output folder for containing Phys output
    G4fs::file_status fst = G4fs::file_status{};
    const G4fs::path outFolderPhysP{InformationKeeper::Instance()->GetPhysOutFolderName().c_str()};
    auto isExist = G4fs::status_known(fst) ? G4fs::exists(fst) : G4fs::exists(outFolderPhysP);
    if (! isExist) {
        G4fs::create_directory(outFolderPhysP);
    } else {
        G4fs::remove_all(outFolderPhysP);//delete old folder
        G4fs::create_directory(outFolderPhysP);
    }
    // create output folder for containing chem_input
    const G4fs::path outFolderP{InformationKeeper::Instance()->GetChemInputFolderName().c_str()};
    isExist = G4fs::status_known(fst) ? G4fs::exists(fst) : G4fs::exists(outFolderP);
    if (! isExist) {
        G4fs::create_directory(outFolderP);
    } else {
        G4fs::remove_all(outFolderP);//delete old folder
        G4fs::create_directory(outFolderP);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ExtractChemListNameFromMacroFile(G4String fileName)
{
    G4String out="";
    std::ifstream file;
    file.open(fileName.c_str());
    if (!file.is_open()) {
        G4String msg = "Error in openning file " + fileName;
        G4Exception("G4String ExtractChemListNameFromMacroFile()", "", FatalException, msg);
    } else {
        G4String line;
        while(std::getline(file, line))
        {
            std::istringstream iss(line);
            G4String flag;
            iss >> flag;
            G4String tvalue;
            iss >> tvalue;
            if (flag == "/dsbandrepair/chem/chemList") out = tvalue;
        }
    }
    file.close();
    return out;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
