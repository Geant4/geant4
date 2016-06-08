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
// $Id: G4VFermiBreakUp.hh,v 1.5 2001/08/01 17:04:56 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4VFermiBreakUp_h
#define G4VFermiBreakUp_h 1

#include "globals.hh"
#include "G4FragmentVector.hh"

class G4VFermiBreakUp 
{
public:
  G4VFermiBreakUp();
  virtual ~G4VFermiBreakUp();
  
private:
  G4VFermiBreakUp(const G4VFermiBreakUp &right);
  
  const G4VFermiBreakUp & operator=(const G4VFermiBreakUp &right);
  G4bool operator==(const G4VFermiBreakUp &right) const;
  G4bool operator!=(const G4VFermiBreakUp &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;
};


#endif


