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
// $Id: G4ios.hh,v 1.6 2001-07-11 10:00:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4ios.hh
//
// ---------------------------------------------------------------
#ifndef included_G4ios
#define included_G4ios

#include "G4Types.hh"

#if defined(OO_DDL_TRANSLATION)
/*
 * stdlib needs to be included before iostream.h
 * during oddlx runs to work around a parser
 * problem of AIX ooddlx v4.0.2 with some versions of
 * AIX system header files.
 */
#include <stdlib.h>
#endif

#include "g4std/iostream"

extern G4std::ostream G4cout;
extern G4std::ostream G4cerr;

#define G4cin G4std::cin
#define G4endl G4std::endl

#endif
