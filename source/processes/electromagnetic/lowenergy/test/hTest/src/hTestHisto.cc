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
//---------------------------------------------------------------------------
//
// ClassName:   hTestHisto
//  
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestHisto.hh"
#include <iomanip>

#include <memory> // for the auto_ptr(T>

#include "AIDA/AIDA.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto* hTestHisto::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto* hTestHisto::GetPointer()
{
  if(!fManager) {
    fManager = new hTestHisto();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto::hTestHisto() 
{
  verbose = 1;
  histName = G4String("histo.hbook");
  ntup = 0;
  nHisto = 1;
  nAbsSaved = 0;
  maxEnergy = 0.0;
  nTuple = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto::~hTestHisto() 
{
  histo.clear(); 
  G4cout << "hTestHisto: Histograms are deleted for " << theName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::BeginOfHisto(G4int num)
{  
  if(0 < verbose) G4cout << "hTestHisto # " << num << " started " << G4endl;
  zend     = 0.0;
  zend2    = 0.0;
  zEvt     = 0.0;
  
  if(0 < nHisto) bookHisto();

  if(verbose > 0) {
    G4cout << "hTestHisto: Histograms are booked and run has been started" 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::EndOfHisto()
{

  G4cout << "hTestHisto: End of run actions are started" << G4endl;

  // Zend average
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    zend2 -= zend*zend;
    G4double sig = 0.0;
    if(zend2 > 0.) sig = std::sqrt(zend2);
    zend2 = sig / std::sqrt(zEvt);
    G4cout<<"========================================================"<<G4endl;
    G4cout << setprecision(4) << "Range(mm)= " << zend/mm 
           << "; Stragling(mm)= " << sig/mm 
           << setprecision(2) << " +- " << zend2/mm 
           << "    " << zEvt << " events for range" << G4endl;
    G4cout<<"========================================================"<<G4endl;
  }  

   // Write histogram file
  if(0 < nHisto || ntup) {
    tree->commit();
    std::cout << "Closing the tree..." << std::endl;
    tree->close();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveEvent()
{
  if(ntup) {
    ntup->addRow();
  }                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveToTuple(const G4String& parname, G4double val)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveToTuple(const G4String& parname,G4double val,G4double defval)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::bookHisto()
{
  G4double zmax = (AbsorberThickness + gap) * NumberOfAbsorbers / mm;
  G4cout << "hTestHisto: Histograms will be saved to the file <" 
         << histName << ">"
         << " AbsThick(mm)= " << AbsorberThickness/mm
         << " Nabs= " << NumberOfAbsorbers
         << " zmax= " << zmax
         << " emax= " << maxEnergy
         << " nHisto= " << nHisto
         << G4endl;

  // Creating the analysis factory
  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

  // Creating the tree factory
  std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new hbook file.
  tree = tf->create(histName,false,false,"hbook");
  G4cout << "Tree store : " << tree->storeName() << G4endl;
 
  histo.resize(nHisto);

  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr< AIDA::IHistogramFactory > hf(af->createHistogramFactory( *tree ));

  // Creating an 1-dimensional histograms in the root directory of the tree

  if(0 < nHisto) histo[0] = hf->createHistogram1D("10",
    "Energy deposit (MeV) in absorber (mm)",NumberOfAbsorbers,0.0,zmax);

  if(1 < nHisto) histo[1] = hf->createHistogram1D("11",
    "Energy (MeV) of delta-electrons",50,0.0,maxEnergy/MeV);

  if(2 < nHisto) histo[2] = hf->createHistogram1D("12",
    "Theta (degrees) of delta-electrons",36,0.0,180.);

  if(3 < nHisto) histo[3] = hf->createHistogram1D("13",
    "Energy (MeV) of secondary gamma",50,0.0,maxEnergy/MeV);

  if(4 < nHisto) histo[4] = hf->createHistogram1D("14",
    "Theta (degrees) of secondary gamma",36,0.0,180.);

  // Creating a tuple factory, whose tuples will be handled by the tree
  std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

  // If using Anaphe HBOOK implementation, there is a limitation on the 
  // length of the variable names in a ntuple
  if(nTuple) ntup = tpf->createITuple( "100", "Range/Energy", 
  "float xend, yend, zend, ltpk, tend, teta, loss, dedx, back, leak, edep" );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddEnergy(G4double edep, G4double z)
{
  if(0 < nHisto) histo[0]->fill((float)z/mm, (float)edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddEndPoint(G4double z)
{
  zend  += z;
  zend2 += z*z;
  zEvt  += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddDeltaElectron(const G4DynamicParticle* elec)
{
  if(1 < nHisto) histo[1]->fill(elec->GetKineticEnergy()/MeV,1.0);
  if(2 < nHisto)
     histo[2]->fill((elec->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddPhoton(const G4DynamicParticle* ph)
{
  if(3 < nHisto) histo[3]->fill(ph->GetKineticEnergy()/MeV,1.0);
  if(4 < nHisto)
     histo[4]->fill((ph->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


