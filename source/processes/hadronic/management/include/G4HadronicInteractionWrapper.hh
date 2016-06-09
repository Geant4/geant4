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
#ifndef G4HadronicInteractionWrapper_h
#define G4HadronicInteractionWrapper_h 1

#include "G4HadronicWhiteBoard.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicInteraction.hh"

class G4HadronicInteractionWrapper
{
  public:
  G4HadFinalState * ApplyInteraction(G4HadProjectile & thePro, 
                                     G4Nucleus & targetNucleus,
                                     G4HadronicInteraction * theInteraction);
};

#endif
