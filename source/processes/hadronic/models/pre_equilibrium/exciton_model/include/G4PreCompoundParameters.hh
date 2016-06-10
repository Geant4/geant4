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
// $Id: G4PreCompoundParameters.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// by V. Lara
//
// 01.04.2008 J. M. Quesada Level density set to A/10 at preequilibrium
// 18.08.2010 V.Ivanchenko make this class as a standard singleton
//

#ifndef G4PreCompoundParameters_h
#define G4PreCompoundParameters_h 1

#include "globals.hh"

class G4PreCompoundParameters
{
public:

  G4PreCompoundParameters();

  ~G4PreCompoundParameters();
 
  inline G4double GetLevelDensity();

  inline G4double Getr0();

  inline G4double GetTransitionsr0();

  inline G4double GetFermiEnergy();

private:


  // Level density parameter
  G4double fLevelDensity;

  // Nuclear radius r0
  G4double fR0;
	
  // Nuclear radius r0 for transitions
  G4double fTransitions_r0;

  // Fermi energy level
  G4double fFermiEnergy;
};

inline G4double G4PreCompoundParameters::GetLevelDensity()
{ 
  return fLevelDensity; 
}
 
inline G4double G4PreCompoundParameters::Getr0()
{ 
  return fR0; 
}

inline G4double G4PreCompoundParameters::GetTransitionsr0()
{ 
  return fTransitions_r0; 
}

inline G4double G4PreCompoundParameters::GetFermiEnergy()
{ 
  return fFermiEnergy; 
}

#endif
