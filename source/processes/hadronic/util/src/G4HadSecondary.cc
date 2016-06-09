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
#include "G4HadSecondary.hh"
#include "G4HadronicException.hh"

G4HadSecondary::G4HadSecondary(G4DynamicParticle * aT, G4double aWeight) :
    theP(aT), theWeight(aWeight), theTime(-1)
{
  if(aT->GetKineticEnergy()<0)
  {
    throw G4HadronicException(__FILE__, __LINE__, 
    "ATTEMPTING TO CREATE A SECONDARY WITH NEGATIVE KINETIC ENERGY.");
  }
}
