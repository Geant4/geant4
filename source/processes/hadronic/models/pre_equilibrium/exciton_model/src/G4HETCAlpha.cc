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
#include "G4HETCAlpha.hh"

G4double G4HETCAlpha::K(const G4Fragment & aFragment)
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
  if (P > 3)
    {
      result = 3.0/(P*(P-1.0)*(P-2.0)*(P-3.0)) * 
	(H*(H-1.0)*(H-2.0)*(H-3.0)*r*r*(r-1.0)*(r-1.0) +
	 2.0*H*(H-1.0)*(H-2.0)*(Pa*r*(1.0-r)*(1.0-r)+Na*r*r*(1.0-r)) +
	 H*(H-1.0)*(Pa*(Pa-1.0)*(1.0-r)*(1.0-r)+4.0*Na*Pa*r*(1.0-r)+Na*(Na-1.0)*r*r) +
	 2*H*(Pa*Na*(Na-1.0)*r+Pa*(Pa-1.0)*Na*(1.0-r)) +
	 Pa*(Pa-1.0)*Na*(Na-1.0));

      result /= 6.0*std::pow((TargetZ/TargetA)*((TargetA-TargetZ)/TargetA),2.0);
    }

  return std::max(0.0,result);

}
