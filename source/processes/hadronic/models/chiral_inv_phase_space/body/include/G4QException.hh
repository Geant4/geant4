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
// $Id: G4QException.hh,v 1.2 2003-11-13 14:40:45 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QException ----------------
//           by Mikhail Kossov, November 2003.
//  header of the class for an Exception used by the CHIPS Model
// ---------------------------------------------------------------

#ifndef G4QException_h
#define G4QException_h

#include <exception>
#include <iostream>
#include "globals.hh"

class G4QException : public std::exception
{
  public:
  G4QException(G4String message);
  const G4String& GetMessage();
  private:
  G4String theMessage;
};

#endif
