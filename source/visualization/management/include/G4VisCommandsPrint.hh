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
// $Id: G4VisCommandsPrint.hh,v 1.3.4.1 2001/06/28 19:16:09 gunter Exp $
// GEANT4 tag $Name:  $
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
