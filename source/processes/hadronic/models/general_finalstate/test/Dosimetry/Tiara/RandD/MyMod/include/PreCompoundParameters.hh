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
// $Id: PreCompoundParameters.hh,v 1.1 2003-10-08 12:32:12 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by V. Lara


#ifndef PreCompoundParameters_h
#define PreCompoundParameters_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class PreCompoundParameters
{
private:
  static PreCompoundParameters thePreCompoundParameters;

  // default constructor
  PreCompoundParameters() :r0(1.5*fermi),Transitions_r0(0.6*fermi){};

public:

    ~PreCompoundParameters() {};
 
    static PreCompoundParameters * GetAddress();

  G4double GetLevelDensity(const G4Fragment& fragm);
 

  G4double Getr0()
  { return r0; }

  G4double GetTransitionsr0()
  { return Transitions_r0; }


  G4double GetFermiEnergy(const G4Fragment& fragm);

private:
    // Nuclear radius r0
    const G4double r0;
	
    // Nuclear radius r0 for transitions
    const G4double Transitions_r0;
};

#endif
