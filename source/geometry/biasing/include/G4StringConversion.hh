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
// $Id: G4StringConversion.hh,v 1.1 2002-10-16 14:30:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Declaration
//
// Declaration description:
//
// This fucntions convert some non character objects into 
// G4Strings

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4StringConversion_hh
#define G4StringConversion_hh G4StringConversion_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Nsplit_Weight;


G4String str(G4int i);
  // convert an integer into a G4String 

G4String str(G4double d);
  // convert a double into a G4String

G4String str(const G4ThreeVector &v);
  // convert a G4ThreeVector into a G4String


G4std::ostream& operator<<(G4std::ostream &out, 
			   const G4Nsplit_Weight &nw);

G4String str(const G4Nsplit_Weight &nw);

#endif
