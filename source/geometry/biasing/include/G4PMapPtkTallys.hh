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
// $Id: G4PMapPtkTallys.hh,v 1.7 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PMapPtkTallys
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PMapPtkTallys_hh 
#define G4PMapPtkTallys_hh

#include "g4std/map"
#include "globals.hh"
#include "g4std/iostream"
#include "G4PTouchableKey.hh"
#include "G4Sigma.hh"

typedef G4std::map<const char *, G4Sigma> G4PMapNameTally;

typedef G4std::map<G4PTouchableKey, G4PMapNameTally, G4PTkComp> G4PMapPtkTallys; 

G4std::ostream& operator<<(G4std::ostream &out, const G4PMapNameTally &tally);
G4std::ostream& operator<<(G4std::ostream &out, const G4PMapPtkTallys &ptktally);

#endif
