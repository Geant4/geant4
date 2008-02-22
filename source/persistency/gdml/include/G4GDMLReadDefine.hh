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
// class G4GDMLReadDefine
//
// Class description:
//
// GDML class for positionings and transformations according to
// specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLREADDEFINE_INCLUDED_
#define _G4GDMLREADDEFINE_INCLUDED_

#include <map>

#include "G4ThreeVector.hh"
#include "G4GDMLRead.hh"

class G4GDMLReadDefine : public G4GDMLRead {

   std::map<G4String,G4ThreeVector*> positionMap;
   std::map<G4String,G4ThreeVector*> rotationMap;
   std::map<G4String,G4ThreeVector*> scaleMap;
   std::map<G4String,G4double> quantityMap;

   void constantRead(const xercesc::DOMElement* const); 
   void matrixRead(const xercesc::DOMElement* const);
   void positionRead(const xercesc::DOMElement* const);
   void rotationRead(const xercesc::DOMElement* const);
   void scaleRead(const xercesc::DOMElement* const);
   void variableRead(const xercesc::DOMElement* const); 
   void quantityRead(const xercesc::DOMElement* const); 
   void defineRead(const xercesc::DOMElement* const);
protected:
   void vectorRead(const xercesc::DOMElement* const,G4ThreeVector&);
   G4String refRead(const xercesc::DOMElement* const);
   G4ThreeVector* getPosition(const G4String&);
   G4ThreeVector* getRotation(const G4String&);
   G4ThreeVector* getScale(const G4String&);
   G4double getQuantity(const G4String&);
   G4RotationMatrix getRotationMatrix(const G4ThreeVector&);
};

#endif
