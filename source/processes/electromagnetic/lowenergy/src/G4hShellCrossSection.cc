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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hShellCrossSection
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20.10.01
//
//
// -------------------------------------------------------------------

#include "G4hShellCrossSection.hh"
#include "G4AtomicTransitionManager.hh"

G4hShellCrossSection::G4hShellCrossSection():G4VhShellCrossSection()
{}


G4hShellCrossSection::~G4hShellCrossSection() 
{}


G4std::vector<G4double>* G4hShellCrossSection::Probabilities(G4int Z, 
                                               G4double kineticEnergy) const 
{
  G4AtomicTransitionManager* transitionManager = 
                             G4AtomicTransitionManager::Instance();
  G4int nShells = transitionManager->NumberOfShells(Z);
  G4std::vector<G4double>* p = new G4std::vector<G4double>;
  G4double norm = 0.0;
  G4double x;
  p->resize(nShells);

  for(G4int shell=0; shell<nShells; shell++) {

    G4double bindingE = transitionManager->Shell(Z, shell)->BindingEnergy();

    // screen function to be replace
    x  = 1./(bindingE + kineticEnergy);
    (*p)[shell] = x*x;
    norm += x;
  }
  if(norm) norm = 1.0/norm;
  for(G4int i=0; i<nShells; i++) {
    x = (*p)[i] * norm;
    (*p)[i] = x;
  }

  return p;
}


