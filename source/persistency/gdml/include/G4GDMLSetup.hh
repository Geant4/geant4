//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GDMLSetup.hh,v 1.10 2007-12-05 15:44:48 ztorzsok Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLSetup
//
// Class description:
//
// GDML class for geometry setup of physical volumes.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLSETUP_INCLUDED_
#define _G4GDMLSETUP_INCLUDED_

#include <map>

#include "G4String.hh"

#include "G4GDMLSolids.hh"

class G4GDMLSetup : public G4GDMLSolids {
private:
   std::map<G4String,G4String> setupMap;
   void setupRead(const xercesc::DOMElement* const element);
public:
   G4String getSetup(const G4String&);
};

#endif
