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

#include "globals.hh"

#include "G4ios.hh"

#include "G4CompetitiveFission.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "../src/G4FissionParameters.cc"

int main()
{

  G4int A, Z;

  G4cout << "Please enter Z and A" << G4endl;
  G4cin >> Z >> A;

  G4int iter;
  G4cout << "Please enter number of iterations " << G4endl;
  G4cin >> iter;
  
  G4double energy;
  G4cout << "Please enter the excitation energy" << G4endl;
  G4cin >> energy;

  G4CompetitiveFission theFission;

  G4int i=0;
  for (i=0; i<iter; i++)
  {
    G4LorentzVector p4(0.,0.,0.,G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( Z, A )+energy);
    G4Fragment aFragment(A,Z,p4);
    G4FragmentVector * theResult = theFission.BreakUp(aFragment);
    int ii=0;
    for(ii=0; ii<theResult->entries(); ii++)
    {
      G4cout <<theResult->at(ii)->GetA() << " "<<theResult->at(ii)->GetZ() << G4endl;
      delete theResult->at(ii);
    }
    delete theResult;
  }
}



















