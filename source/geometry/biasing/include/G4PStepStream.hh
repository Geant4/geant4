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
// $Id: G4PStepStream.hh,v 1.5 2002-04-10 13:13:06 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Declarations
//
// Declaration description:
//
// declarations of streams for G4PTouchableKey and G4PStep

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PStepStream_hh
#define G4PStepStream_hh G4PStepStream_hh

#include "G4PTouchableKey.hh"
#include "G4PStep.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4PTouchableKey &tk);
G4std::ostream& operator<<(G4std::ostream &out, const G4PStep &ps);

#endif
