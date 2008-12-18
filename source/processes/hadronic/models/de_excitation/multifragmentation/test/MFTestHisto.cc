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
// Igor Pshenichnov 10.05.2008
// This main{} used to run standalone tests for multifragmentation of excited nuclei
// ROOT installation is reqired
  

#include "G4StatMF.hh"
#include "G4NucleiProperties.hh"
#include "MFTestHistoManager.hh"
#include "Randomize.hh"

#ifdef G4ANALYSIS_USE_ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#endif

#include "MFTestHistoManager.cc"

int main()
{
  // The user will be asked for the charge and mass of the hot nuclear system to
  // simulate decay. 

  MFTestHistoManager histoManager;

  // Histograms will be booked now.
  histoManager.BookHisto();

  G4int A = histoManager.GetA();
  G4int Z = histoManager.GetZ();
    
  G4StatMF model;

  G4LorentzVector p4null(0.0,0.0,0.0,0.0);


  G4double totalNumFragments = 0.;
  G4int thisEventNumFragments = 0;

  G4int iterations = histoManager.GetIterations();

  std::cout.setf(std::ios::scientific,std::ios::floatfield);

  for (G4int i = 0; i < iterations; ++i) {

  G4double energy = (histoManager.GetUpEn() - histoManager.GetLowEn())*G4UniformRand() 
                                                             + histoManager.GetLowEn();

  G4double NuclearMass = G4NucleiProperties::GetNuclearMass(A,Z) + energy;
  G4LorentzVector p4(0.0,0.0,0.0,NuclearMass);
  G4Fragment aFragment(A,Z,p4);
                                                          
  G4LorentzVector p4cons;
  p4cons = p4null;
      
  G4FragmentVector * theResult = model.BreakItUp(aFragment);
  thisEventNumFragments = theResult->size();
  if ( i%100 == 0 ) std::cout << "### event " << i << " at " 
			      << energy/A << " MeV/nucleon" 
			      << " with "<< thisEventNumFragments <<" fragments #####\n";
  histoManager.fill(0,energy/A,thisEventNumFragments);

  totalNumFragments += thisEventNumFragments;
 
  for (G4FragmentVector::iterator i = theResult->begin(); i != theResult->end(); ++i)
    {
      p4cons += (*i)->GetMomentum();

      G4double thisFragmentZ = (*i)->GetZ();
      G4double thisFragmentA = (*i)->GetA();
      histoManager.fill2(0,energy/A,thisFragmentZ,1.);
      histoManager.fill2(1,energy/A,thisFragmentA,1.);
     
      histoManager.fill(3,thisFragmentA,thisFragmentZ);
      histoManager.fill(4,thisFragmentA,1.);
       
      if ( (int)thisFragmentZ == 6 && (int)thisFragmentA == 10) {
	histoManager.fill(1,energy/A,1.);
      }
      if ( (int)thisFragmentZ == 6 && (int)thisFragmentA == 11) {
	histoManager.fill(2,energy/A,1.);
      }
                      
      delete (*i);
    }
  delete theResult;
  }

  std::cout << " ------------------------------------------------------------------  " << '\n';
  std::cout << " Average fragment multiplicty for the whole energy range: " 
	    << totalNumFragments/iterations << '\n';
  
  histoManager.NormalizeHisto();
  histoManager.CleanHisto();

  return 0;
}


      
