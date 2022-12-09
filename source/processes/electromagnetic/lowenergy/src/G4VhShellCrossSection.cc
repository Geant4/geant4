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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VhShellCrossSection
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// History:
// -----------
// 20 Oct 2001 V.Ivanchenko   1st implementation
// 24 Oct 2001 MGP            Minor clean-up
// 29 Oct 2001 VI             Add delta energy
//
// -------------------------------------------------------------------

#include "G4VhShellCrossSection.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VhShellCrossSection::G4VhShellCrossSection(const G4String& xname)
 :name(xname)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VhShellCrossSection::~G4VhShellCrossSection() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VhShellCrossSection :: SetTotalCS(G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VhShellCrossSection::SelectRandomShell(G4int Z, 
                                               G4double incidentEnergy,
					       G4double mass, 
					       G4double deltaEnergy,
					       const G4Material* mat)
//  returns the shell ionized if the shell exists. If the shell is 
// not counted, it returns -1
{
  std::vector<G4double> p = 
    Probabilities(Z,incidentEnergy,mass,deltaEnergy,mat);
  G4int shell = -1;
  G4int nShells = (G4int)p.size();
  G4double q = G4UniformRand();
  for (G4int i=0; i<nShells; ++i) {
    
    if (p[i] >= q) {
      shell = i;
      break;
    }
    q -= p[i];
  }
  return shell;
}
