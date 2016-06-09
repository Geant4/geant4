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
#include "G4HETCHe3.hh"

G4double G4HETCHe3::K(const G4Fragment & aFragment)
{
  if (GetStage() != 1) return 1.0;
  // Number of protons in projectile
  G4double Pa = static_cast<G4int>(aFragment.GetParticleDefinition()->GetPDGCharge());
  // Number of neutrons in projectile 
  G4double Na = aFragment.GetParticleDefinition()->GetBaryonNumber();
  G4double TargetA = aFragment.GetA() - Na;
  G4double TargetZ = aFragment.GetZ() - Pa;
  Na -= Pa;
  G4double r = TargetZ/TargetA;


  
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();

  G4double result = 0.0;
  if (P > 2)
    {
      result = 3.0/(P*(P-1.0)*(P-2.0)) * 
	(H*(H-1.0)*(H-2.0)*r*r*(r-1.0) +
	 H*(H-1.0)*(2.0*Na*r*(1.0-r)+Pa*r*r) +
	 H*(Pa*(Pa-1.0)*(r-1.0)+2.0*Na*Pa*r) +
	 Pa*Na*(Pa-1.0));

      result /= 3.0*pow(TargetZ/TargetA, 2.0) * ((TargetA-TargetZ)/TargetA);
    }

  return std::max(0.0,result);

}
