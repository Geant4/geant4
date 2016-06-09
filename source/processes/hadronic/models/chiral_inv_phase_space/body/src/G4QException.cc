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
// $Id: G4QException.cc,v 1.7 2004/03/25 10:44:58 gunter Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
//      ---------------- G4QException ----------------
//           by Mikhail Kossov, November 2003.
//      class for an Exception used by the CHIPS Model
// ------------------------------------------------------------

#include "G4QException.hh"

G4QException::G4QException(G4String message) {theMessage = message;}

G4QException::~G4QException() throw(){}

const G4String& G4QException::GetMessage() {return theMessage;}

