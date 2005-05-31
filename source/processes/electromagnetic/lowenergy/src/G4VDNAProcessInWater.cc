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
// $Id: G4VDNAProcessInWater.cc,v 1.1 2005-05-31 09:58:40 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VDNAProcessInWater.hh"

G4VParticleChange *                      G4VDNAProcessInWater :: PostStepDoIt(const G4Track & aTrack, const G4Step & aStep)
{
 ValidateInWater(aTrack);
 
 return G4VLowEnergyTestableDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void                                     G4VDNAProcessInWater :: ValidateInWater(const G4Track & aTrack) const
{
 G4Material * theMaterial(aTrack.GetMaterial());
 
 // Water validation consists in checking that density is 1. and that Oxigen and Hydrogen are the only two elements of the material

 if (fabs(theMaterial->GetDensity()*cm3/g-1.)<0.1)
  if (theMaterial->GetNumberOfElements()==2)
  {
   const G4Element * element1(theMaterial->GetElement(0));
   const G4Element * element2(theMaterial->GetElement(1));

   if ((element1->GetZ()==1. && element2->GetZ()==8.) ||
       (element1->GetZ()==8. && element2->GetZ()==1.))
    return;
  }
  
 G4String message;
 message="G4VDNAProcessInWater::ValidateInWater - Process ";
 message+=GetProcessName();
 message+=" happened in material ";
 message+=theMaterial->GetName();
  
 G4Exception(message);
}

