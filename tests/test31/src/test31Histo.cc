//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//---------------------------------------------------------------------------
//
// ClassName:   test31Histo
//  
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31Histo.hh"
#include "g4std/iomanip"

#include <memory> // for the auto_ptr(T>

#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"

#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::GetPointer()
{
  if(!fManager) {
    fManager = new test31Histo();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::test31Histo() 
{
  verbose = 1;
  histName = G4String("histo.hbook");
  ntup = 0;
  nHisto = 1;
  maxEnergy = 0.0;
  nTuple = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::~test31Histo() 
{
  histo.clear(); 
  G4cout << "test31Histo: Histograms are deleted for " << theName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::BeginOfHisto(G4int num)
{  
  if(0 < verbose) G4cout << "test31Histo # " << num << " started " << G4endl;
  zend     = 0.0;
  zend2    = 0.0;
  zEvt     = 0.0;
  
  if(0 < nHisto) bookHisto();

  if(verbose > 0) {
    G4cout << "test31Histo: Histograms are booked and run has been started" 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::EndOfHisto()
{

  G4cout << "test31Histo: End of run actions are started" << G4endl;

  // Zend average
  G4cout<<"========================================================"<<G4endl;
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    zend2 -= zend*zend;
    G4double sig = 0.0;
    if(zend2 > 0.) sig = sqrt(zend2);
    zend2 = sig / sqrt(zEvt);
    G4cout << setprecision(4) << "Range(mm)= " << zend/mm 
           << "; Stragling(mm)= " << sig/mm 
           << setprecision(2) << " +- " << zend2/mm 
           << "    " << zEvt << " events for range" << G4endl;
  }  
  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;
  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*(G4double)n_step;
  G4cout                    << "Number of events          " << n_evt <<G4endl; 
  G4cout << setprecision(4) << "Average number of e-      " << xe << G4endl; 
  G4cout << setprecision(4) << "Average number of gamma   " << xg << G4endl; 
  G4cout << setprecision(4) << "Average number of e+      " << xp << G4endl; 
  G4cout << setprecision(4) << "Average number of steps   " << xs << G4endl; 
  G4cout<<"========================================================"<<G4endl;


   // Write histogram file
  if(0 < nHisto) {
    tree->commit();
    G4std::cout << "Closing the tree..." << G4std::endl;
    tree->close();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveEvent()
{
  if(ntup) {
    ntup->addRow();
  }                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname, G4double val)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname,G4double val, G4double)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::bookHisto()
{
  G4double zmax = (AbsorberThickness + gap) * NumberOfAbsorbers / mm;
  G4cout << "test31Histo: Histograms will be saved to the file <" 
         << histName << ">"
         << " AbsThick(mm)= " << AbsorberThickness/mm
         << " Nabs= " << NumberOfAbsorbers
         << " zmax= " << zmax
         << " nHisto= " << nHisto
         << G4endl;

  n_evt  = 0;
  n_elec = 0;
  n_posit= 0;
  n_gam  = 0;
  n_step = 0;

  // Creating the analysis factory
  G4std::auto_ptr< IAnalysisFactory > af( AIDA_createAnalysisFactory() );

  // Creating the tree factory
  G4std::auto_ptr< ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new hbook file.
  tree = tf->create(histName,false,false,"hbook");
  G4cout << "Tree store : " << tree->storeName() << G4endl;
 
  histo.resize(nHisto);

  // Creating a histogram factory, whose histograms will be handled by the tree
  G4std::auto_ptr< IHistogramFactory > hf(af->createHistogramFactory( *tree ));

  // Creating an 1-dimensional histograms in the root directory of the tree

  if(0 < nHisto) histo[0] = hf->create1D("10",
    "Energy deposit (MeV) in absorber (mm)",NumberOfAbsorbers,0.0,zmax);

  if(1 < nHisto) histo[1] = hf->create1D("11",
    "Energy (MeV) of secondary electrons",50,0.0,maxEnergy/MeV);

  if(2 < nHisto) histo[2] = hf->create1D("12",
    "Theta (degrees) of delta-electrons",36,0.0,180.);

  if(3 < nHisto) histo[3] = hf->create1D("13",
    "Energy (MeV) of secondary gamma",50,0.0,maxEnergy/MeV);

  if(4 < nHisto) histo[4] = hf->create1D("14",
    "Theta (degrees) of secondary gamma",36,0.0,180.);

  // Creating a tuple factory, whose tuples will be handled by the tree
  G4std::auto_ptr< ITupleFactory > tpf( af->createTupleFactory( *tree ) );

  // If using Anaphe HBOOK implementation, there is a limitation on the 
  // length of the variable names in a ntuple
  if(nTuple) ntup = tpf->create( "100", "Range/Energy", 
  "float xend, yend, zend, ltpk, tend, teta, loss, dedx, back, leak, edep" );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEnergy(G4double edep, G4double z)
{
  if(0 < nHisto) histo[0]->fill((float)z/mm, (float)edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEndPoint(G4double z)
{
  zend  += z;
  zend2 += z*z;
  zEvt  += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddDeltaElectron(const G4DynamicParticle* elec)
{
  n_elec++;
  if(1 < nHisto) histo[1]->fill((float)elec->GetKineticEnergy()/MeV,1.0);
  if(2 < nHisto)
     histo[2]->fill((float)(elec->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddPhoton(const G4DynamicParticle* ph)
{
  n_gam++;
  if(3 < nHisto) histo[3]->fill((float)ph->GetKineticEnergy()/MeV,1.0);
  if(4 < nHisto)
     histo[4]->fill((float)(ph->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

