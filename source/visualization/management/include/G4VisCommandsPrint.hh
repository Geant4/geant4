// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsPrint.hh,v 1.3 1999-12-15 14:54:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/print/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSPRINT_HH
#define G4VISCOMMANDSPRINT_HH

#include "globals.hh"

////////////////////////////////////////////////////  /vis~/print/...  ////
//vis \hline
//vis /vis~/print/ &&
//vis ...menu of printing possibilities, also controlled by ``Verbose''
//vis parameter. \\%
class G4VisCommandPrint {
public:
  G4String GetCommandName () const {return "/vis~/print/";}
  G4String GetGuidance () const {
    return "...menu of printing possibilities, also controlled by \"Verbose\" "
      "parameter.";
  }
};

#endif
