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
// File name:     G4VhShellCrossSection
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20.10.01
//
//
// -------------------------------------------------------------------

#include "G4VhShellCrossSection.hh"


G4VhShellCrossSection::G4VhShellCrossSection()
{}


G4VhShellCrossSection::~G4VhShellCrossSection() 
{}


G4int G4VhShellCrossSection::SelectRandomShell(G4int Z, 
                                               G4double kineticEnergy) const 
{
  G4std::vector<G4double>* p = Probabilities(Z, kineticEnergy);
  G4int shell = 0;
  if(p) {
    G4int nShells = p->size();
    G4double q = G4UniformRand();
    for (shell=0; shell<nShells; shell++) {
       
      if((*p)[shell] >= q) break;
      q -= (*p)[shell];

    } 
    delete p;
  }
  return shell;
}


