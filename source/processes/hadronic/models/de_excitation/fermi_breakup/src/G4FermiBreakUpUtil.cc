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

#include "G4FermiBreakUpUtil.hh"
#include "G4FermiFragment.hh"
#include "G4NuclearRadii.hh"
#include "G4PhysicalConstants.hh"

namespace G4FermiBreakUpUtil {

  const G4double deltaR = 0.6*CLHEP::fermi;
  const G4double coeff = 0.9;

  // Coulomb barrier
  G4double CoulombBarrier(const G4int Z1, const G4int A1,
                          const G4int Z2, const G4int A2, const G4double exc) {
    const G4double r1 = G4NuclearRadii::RadiusCB(Z1, A1);
    const G4double r2 = G4NuclearRadii::RadiusCB(Z2, A2);
    G4double CB = coeff*CLHEP::elm_coupling*(Z1*Z2)/(r1 + r2 - deltaR);
    if(exc > 0.0) { CB /= (1.0 + std::sqrt(exc/((2*(A1 + A2))*CLHEP::MeV))); }
    return CB;
  }

  // 2-body decay probability
  G4double Probability(const G4int A,
                       const G4FermiFragment* f1, const G4FermiFragment* f2,
                       const G4double mass, const G4double exc) {
    G4double prob = 0.0;
    const G4double mass1 = f1->GetTotalEnergy();
    const G4double mass2 = f2->GetTotalEnergy();
    const G4double bCouloumb =
      CoulombBarrier(f1->GetZ(), f1->GetA(), f2->GetZ(), f2->GetA(), exc);
    if (mass < mass1 + mass2 + bCouloumb)
      return prob;

    // free energy
    const G4double e = mass - mass1 - mass2;

    // Spin factor S_n
    const G4int S_n = (std::abs(f1->TwoSpinParity())+1)
      *(std::abs(f2->TwoSpinParity())+1); 

    // mass factors
    const G4double x = mass1*mass2/(mass1 + mass2);
    G4double massFactor = x*std::sqrt(x);

    // Permutation Factor G_n - search for identical fragments
    G4double G_n = (f1 == f2) ? 0.5 : 1.0;

    prob = (A*S_n) * massFactor*G_n*std::sqrt(e);
    //G4cout << "prob= " << prob << " Coeff= " << Coeff << G4endl;
    return prob;
  }

  G4bool CheckSpinParity(const G4FermiFragment* f1, const G4FermiFragment* f2,
			 const G4FermiFragment* f3)
  {
    // check parity
    G4int spin1 = f1->TwoSpinParity();
    G4int spin2 = f2->TwoSpinParity();
    G4int spin3 = f3->TwoSpinParity();
    if ((spin3 > 0 && spin1*spin2 < 0) || (spin3 < 0 && spin1*spin2 > 0))
      return false;

    return true;
  }
}
