////////////////////////////////////////////////////////////////////////////////
//
#include "MLAnalysisManager.hh"
#include "MLFluenceAnalyser.hh"
#include "MLDoseAnalyser.hh"
#include "MLNIELAnalyser.hh"
#include "MLPHSAnalyser.hh"
#include "MLAnalysisMessenger.hh"
#include "MLNIELFunction.hh"
#include "MLGeometryConstruction.hh"
#include "MLPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#ifdef USEHBOOK
 #include "HbookHistogram.hh"
#else
 #include "MLHisto1D.hh"
#endif

#include "G4RunManager.hh"
#include <iomanip>  
#include "G4UnitsTable.hh"

#include <strstream>

MLAnalysisManager* MLAnalysisManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
MLAnalysisManager::MLAnalysisManager (int ,char** )
  :nFactor(1./cm2), filenamePrefix("mulassis")
{
//
//
// Get the geometry and physics list.
//
  G4RunManager *runManager = G4RunManager::GetRunManager();
  geometry =
    (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());
//
//
// Define the messenger and the analysis system.
//
  analysisMessenger = new MLAnalysisMessenger(this);
  fluenceAnalyser   = new MLFluenceAnalyser(this, geometry);
  doseAnalyser      = new MLDoseAnalyser(this, geometry);
  NIELAnalyser      = new MLNIELAnalyser(this, geometry);
  PHSAnalyser       = new MLPHSAnalyser(this, geometry);
}
////////////////////////////////////////////////////////////////////////////////
//
MLAnalysisManager::~MLAnalysisManager ()
{
  delete analysisMessenger;
  delete fluenceAnalyser;
  delete doseAnalyser;
  delete NIELAnalyser;
  delete PHSAnalyser;
  G4cout <<"Analyis Manager deleting... " <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
MLAnalysisManager* MLAnalysisManager::getInstance (int argc,char** argv)
{
  if (instance == 0) instance = new MLAnalysisManager(argc,argv);
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
/*
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void MLAnalysisManager::BeginOfRunAction ()
{
#ifdef USEHBOOK
  G4String hbookfilename = filenamePrefix + ".hbook";
  theHbookManager.SetFilename(hbookfilename);
  //  theHbookManager.SetFilename(strcat(filenamePrefix,".hbook"));
#endif
//
//
// Note that the order in which these member functions are executed is
// important since we need to keep track of the data-block number for the
// CSV file, which is calculated here.
//
  G4int i = 0;
  doseAnalyser->BeginOfRunAction(nFactor);
  i = doseAnalyser->GetAddedNBlocks(i);
  PHSAnalyser->BeginOfRunAction(nFactor);
  i = PHSAnalyser->GetAddedNBlocks(i);
  NIELAnalyser->BeginOfRunAction(nFactor);
  i = NIELAnalyser->GetAddedNBlocks(i);
  fluenceAnalyser->BeginOfRunAction(nFactor);
  i = fluenceAnalyser->GetAddedNBlocks(i);
}
////////////////////////////////////////////////////////////////////////////////
//
//   This member is called at the end of each run

void MLAnalysisManager::EndOfRunAction (G4double nevent)
{
#ifdef USEHBOOK
  theHbookManager.WriteOut();
#else
  //  if (nevent > 999) PrintHistograms(nevent);
  NormaliseOutput(nevent);
  PrintOutput(nevent);
#endif

}
////////////////////////////////////////////////////////////////////////////////
//
//   This member is called at the end of selected event via event action

void MLAnalysisManager::EndOfEventAction (G4double nevent)
{
#ifdef USEHBOOK
  theHbookManager.Backup();
#else
  NormaliseOutput(nevent);
  PrintOutput(nevent);
#endif
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
void MLAnalysisManager::NormaliseOutput (G4double nevent)
{
  fluenceAnalyser->NormaliseOutput(nevent);
  NIELAnalyser->NormaliseOutput(nevent);
  doseAnalyser->NormaliseOutput(nevent);
  PHSAnalyser->NormaliseOutput(nevent);
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4ParticleTable.hh"
#include "MLPhysicsList.hh"
#include "MLRunManager.hh"
#include "MLVersion.hh"

void MLAnalysisManager::PrintOutput (G4double nevent)
{
  // Open a CSV format file named fileName1.
  //
  MLVersion codeVersion = MLVersion();
  G4RunManager* runManager   = G4RunManager::GetRunManager();
  const MLPhysicsList* physicsList = dynamic_cast <const MLPhysicsList*>
    (runManager->GetUserPhysicsList());
  const MLPrimaryGeneratorAction *primaryGenerator = dynamic_cast
    <const MLPrimaryGeneratorAction*>
    (runManager->GetUserPrimaryGeneratorAction());

  G4String CSVfilename = filenamePrefix + ".csv";

  //  CSVFile.open(CSVfilename, ios::out);
  CSVFile.open(CSVfilename);
  if (!CSVFile.rdbuf()->is_open()) {

    G4cerr <<G4endl <<"CSV output file cannot be opened !" <<G4endl;
    G4cerr <<"No CSV output is being produced." <<G4endl;
  } else {
//
//
// Note that the order in which these member functions are executed is
// important since we need to keep track of the data-block number for the
// CSV file.
//
    CSVFile <<*fluenceAnalyser;
    CSVFile <<*NIELAnalyser;
    CSVFile <<*PHSAnalyser;
    CSVFile <<*doseAnalyser;
    CSVFile.close();
  }
//
//
// Open the report output file and output.
//
  G4String RPTfilename = filenamePrefix + ".rpt";
  
  //  RPTFile.open(RPTfilename, ios::out);
  RPTFile.open(RPTfilename);
  if (!RPTFile.rdbuf()->is_open()) {

    G4cerr <<G4endl <<"Report output file cannot be opened !" <<G4endl;
    G4cerr <<"No RPT output is being produced." <<G4endl;
  } else {
    RPTFile <<codeVersion;
    RPTFile <<*geometry;
    RPTFile <<*physicsList;

    RPTFile <<G4endl;
    RPTFile <<"-------------------------------------------------------------"
            <<G4endl;
    RPTFile <<"Incident Particle Definition:" <<G4endl;
    RPTFile <<"-------------------------------------------------------------"
            <<G4endl;
    RPTFile <<G4endl;
    RPTFile <<"The incident particle type:       "
            <<primaryGenerator->GetParticleGun()
              ->GetParticleDefinition()->GetParticleName() <<G4endl;
    RPTFile <<"The incident position is:         "
            <<primaryGenerator->GetParticleGun()
              ->GetParticlePosition() <<G4endl;
    RPTFile <<"The number of particles fired is: " <<nevent <<G4endl;
    RPTFile <<"The incident particle fluence:    " <<nFactor*cm2
            <<" particles/cm2" <<G4endl <<G4endl;
    
    RPTFile <<*fluenceAnalyser;
    RPTFile <<*NIELAnalyser;
    RPTFile <<*PHSAnalyser;
    RPTFile <<*doseAnalyser;
    RPTFile <<"-------------------------------------------------------------" 
	    << G4endl;
    RPTFile <<"Simulation time: "
            <<(static_cast <MLRunManager*>(runManager))->GetTimeUsed()
            <<" seconds." <<G4endl;

    RPTFile.close();
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif
