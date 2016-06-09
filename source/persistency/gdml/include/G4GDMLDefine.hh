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
// $Id: G4GDMLDefine.hh,v 1.11 2007/11/30 11:58:46 ztorzsok Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
// class G4GDMLDefine
//
// Class description:
//
// GDML class for positionings and transformations according to
// specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLDEFINE_INCLUDED_
#define _G4GDMLDEFINE_INCLUDED_

#include <map>

#include "G4ThreeVector.hh"
#include "G4GDMLBase.hh"

class G4GDMLDefine : public G4GDMLBase {
private:
   std::map<G4String,G4ThreeVector*> positionMap;
   std::map<G4String,G4ThreeVector*> rotationMap;
   std::map<G4String,G4ThreeVector*> scaleMap;

   void constantRead(const xercesc::DOMElement* const); 
   void positionRead(const xercesc::DOMElement* const);
   void rotationRead(const xercesc::DOMElement* const);
   void scaleRead(const xercesc::DOMElement* const);
   void variableRead(const xercesc::DOMElement* const); 
   void defineRead(const xercesc::DOMElement* const);
protected:
   G4ThreeVector* getPosition(const G4String&);
   G4ThreeVector* getRotation(const G4String&);
   G4ThreeVector* getScale(const G4String&);
};

#endif
