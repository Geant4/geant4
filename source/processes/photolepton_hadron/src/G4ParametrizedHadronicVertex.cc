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
// $Id: G4ParametrizedHadronicVertex.cc,v 1.5 2001/10/24 17:47:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// --------------------------------------------------------------
#include "G4ParametrizedHadronicVertex.hh"

G4VParticleChange * G4ParametrizedHadronicVertex::
ApplyYourself(G4Nucleus & theTarget, const G4Track &thePhoton)
{   
    G4double theKineticEnergy = thePhoton.GetKineticEnergy();
    if(RandBit::shootBit())
    {
      if(theKineticEnergy<20*GeV) return theLowEPionMinus.ApplyYourself(thePhoton, theTarget);
      return theHighEPionMinus.ApplyYourself(thePhoton, theTarget);
    }
    else
    {
      if(theKineticEnergy<20*GeV) return theLowEPionPlus.ApplyYourself(thePhoton, theTarget);
      return theHighEPionPlus.ApplyYourself(thePhoton, theTarget);
    }
    return 0;
}
