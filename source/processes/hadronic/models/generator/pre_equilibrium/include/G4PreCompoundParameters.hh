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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PreCompoundParameters.hh,v 1.8 2002/06/06 17:10:07 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// by V. Lara


#ifndef G4PreCompoundParameters_h
#define G4PreCompoundParameters_h 1

#include "globals.hh"

class G4PreCompoundParameters
{
private:
    static G4PreCompoundParameters thePreCompoundParameters;

    // default constructor
    G4PreCompoundParameters() : theLevelDensity(0.125/MeV),
      r0(1.5*fermi),Transitions_r0(0.6*fermi),FermiEnergy(35.0*MeV) 
	{}

public:

    ~G4PreCompoundParameters() {};
 
    static G4PreCompoundParameters * GetAddress();

    G4double GetLevelDensity()
	{ return theLevelDensity; }
 

    G4double Getr0()
	{ return r0; }

    G4double GetTransitionsr0()
	{ return Transitions_r0; }


    G4double GetFermiEnergy()
	{ return FermiEnergy; }

private:
    // Level density parameter
    const G4double theLevelDensity;

    // Nuclear radius r0
    const G4double r0;
	
    // Nuclear radius r0 for transitions
    const G4double Transitions_r0;

    // Fermi energy level
    const G4double FermiEnergy;
  
};

#endif
