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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4FermiDecayProbability.hh"
#include "G4FermiFragment.hh"
#include "G4PhysicalConstants.hh"

G4FermiDecayProbability::G4FermiDecayProbability()
{}

G4double 
G4FermiDecayProbability::ComputeProbability(G4int, G4int A, G4int spin,
					    G4double etot, 
					    const G4FermiFragment* f1,
					    const G4FermiFragment* f2) const
{
  G4double prob = 0.0;
  G4double mass1 = f1->GetTotalEnergy();
  G4double mass2 = f2->GetTotalEnergy();
  G4double bCouloumb = f1->GetCoulombBarrier(f2->GetA(), f2->GetZ(), 0.0);
  if(etot <= mass1 + mass2 + bCouloumb) { return prob; }
  
  //G4cout << "ComputeProbability M1= " << mass1 << " M2= " << mass2 << G4endl;
  G4double ekin = etot - mass1 - mass2;

  // mass factors
  G4double x = mass1*mass2/(mass1 + mass2);
  G4double massFactor = x*std::sqrt(x);

  // Spin factor S_n
  G4double S_n = 1.0;
  if(spin >= 0) {
    G4int spin1 = f1->GetSpin();
    G4int spin2 = f2->GetSpin();
    if(spin1 >= 0 && spin2 >= 0) {
      S_n = (spin1+1)*(spin2+1); 
    }
  }
    
  // Permutation Factor G_n
  // search for identical fragments
  G4double G_n = (f1 == f2) ? 0.5 : 1.0;

  prob = A*massFactor*S_n*G_n*std::sqrt(ekin);

  //G4cout << "prob= " << prob << " Coeff= " << Coeff << G4endl;
  return prob; 
}
