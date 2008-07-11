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

class G4GDMLMatrix
{
   G4double *m;
   size_t rows,cols;

 public:

   G4GDMLMatrix();
   G4GDMLMatrix(size_t rows0,size_t cols0);
   ~G4GDMLMatrix();

   void set(size_t r,size_t c,G4double a);
   G4double get(size_t r,size_t c) const;
   size_t getRows() const;
   size_t getCols() const;
};

class G4GDMLReadDefine : public G4GDMLRead
{
 private:

   std::map<G4String,G4double> quantityMap;
   std::map<G4String,G4ThreeVector> positionMap;
   std::map<G4String,G4ThreeVector> rotationMap;
   std::map<G4String,G4ThreeVector> scaleMap;
   std::map<G4String,G4GDMLMatrix> matrixMap;

   void constantRead(const xercesc::DOMElement* const); 
   void matrixRead(const xercesc::DOMElement* const);
   void positionRead(const xercesc::DOMElement* const);
   void rotationRead(const xercesc::DOMElement* const);
   void scaleRead(const xercesc::DOMElement* const);
   void variableRead(const xercesc::DOMElement* const); 
   void quantityRead(const xercesc::DOMElement* const); 
   void defineRead(const xercesc::DOMElement* const);

 protected:

   G4RotationMatrix getRotationMatrix(const G4ThreeVector&);
   void vectorRead(const xercesc::DOMElement* const,G4ThreeVector&);
   G4String refRead(const xercesc::DOMElement* const);

 public:

   G4double getConstant(const G4String&);
   G4double getVariable(const G4String&);
   G4double getQuantity(const G4String&);
   G4ThreeVector getPosition(const G4String&);
   G4ThreeVector getRotation(const G4String&);
   G4ThreeVector getScale(const G4String&);
   G4GDMLMatrix getMatrix(const G4String&);
};

#endif
