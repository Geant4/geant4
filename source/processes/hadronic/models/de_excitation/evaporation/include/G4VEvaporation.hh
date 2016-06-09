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
// $Id: G4VEvaporation.hh,v 1.1 2003/08/26 18:30:51 lara Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) written from G4Evaporation.hh (May 1998)
//



#ifndef G4VEvaporation_h
#define G4VEvaporation_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VEvaporation 
{
public:
  G4VEvaporation() {};
  virtual ~G4VEvaporation() {}; // *

private:  
  G4VEvaporation(const G4VEvaporation &right);

  const G4VEvaporation & operator=(const G4VEvaporation &right);
  G4bool operator==(const G4VEvaporation &right) const;
  G4bool operator!=(const G4VEvaporation &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;


};


#endif
