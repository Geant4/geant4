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
// $Id: G4PreCompoundParameters.hh,v 1.3 2006/06/29 20:58:32 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
