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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MicrobeamHistoManager.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamHistoManager::MicrobeamHistoManager()
:af(0),tree(0),factoryOn(false)
{
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  af = AIDA_createAnalysisFactory();
  if(!af) {
    G4cout << " MicrobeamHistoManager::MicrobeamHistoManager() :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }
#endif 
 
  fileName[0] = "microbeam";
  fileType    = "root";
  fileOption  = "";
  ntupl0=0;
  ntupl1=0;
  ntupl2=0;
  ntupl3=0;
  ntupl4=0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamHistoManager::~MicrobeamHistoManager()
{
#ifdef G4ANALYSIS_USE  
  delete af;
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamHistoManager::book()
{
#ifdef G4ANALYSIS_USE
  if(!af) return;

  // Creating a tree mapped to an hbook file.
  fileName[1] = fileName[0] + "." + fileType;
  G4bool readOnly  = false;
  G4bool createNew = true;
  AIDA::ITreeFactory* tf  = af->createTreeFactory();
  tree = tf->create(fileName[1], fileType, readOnly, createNew, fileOption);
  delete tf;
  if(!tree) {
    G4cout << "MicrobeamHistoManager::book() :" 
           << " problem creating the AIDA tree with "
           << " storeName = " << fileName[1]
           << " storeType = " << fileType
           << " readOnly = "  << readOnly
           << " createNew = " << createNew
           << " options = "   << fileOption
           << G4endl;
    return;
  }

  // Creating a histogram & ntuplr factory
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
  AIDA::ITupleFactory* ntf = af->createTupleFactory(*tree);
 
  ntupl0 = ntf->create( "ntuple0", "Stopping power", "double e, sp");
  ntupl1 = ntf->create( "ntuple1", "Beam position", "double x, y");
  ntupl2 = ntf->create( "ntuple2", "Range", "double x, y, z");
  ntupl3 = ntf->create( "ntuple3", "Doses", "double doseN, doseC");
  ntupl4 = ntf->create( "ntuple4", "3D", "double x, y, z, doseV");
  factoryOn = true;

  delete hf;
  delete ntf;
      
  if (factoryOn) 
     G4cout << "\n----> Histogram Tree is opened in " << fileName[1] << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamHistoManager::save()
{
#ifdef G4ANALYSIS_USE
  if (factoryOn) {
  
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;

    delete tree;
    tree = 0;
    factoryOn = false;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamHistoManager::FillNtuple(G4int nt, G4int column, G4double value)
{
  if (nt >= MaxNtupl) {
    G4cout << "---> warning from MicrobeamHistoManager::FillNtuple() : Ntuple " << nt
           << " does not exist " << column << value << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(nt==0) ntupl0->fill(column, value);
  if(nt==1) ntupl1->fill(column, value);
  if(nt==2) ntupl2->fill(column, value);
  if(nt==3) ntupl3->fill(column, value);
  if(nt==4) ntupl4->fill(column, value);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrobeamHistoManager::AddRowNtuple(G4int nt)
{
  if (nt >= MaxNtupl) {
    G4cout << "---> warning from MicrobeamHistoManager::AddRowNtuple() : Ntuple " << nt
           << " do not exist" << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(nt==0) ntupl0->addRow();
  if(nt==1) ntupl1->addRow();
  if(nt==2) ntupl2->addRow();
  if(nt==3) ntupl3->addRow();
  if(nt==4) ntupl4->addRow();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
