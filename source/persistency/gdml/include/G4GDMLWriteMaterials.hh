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
// $Id: G4GDMLWriteMaterials.hh,v 1.12 2008/07/16 15:46:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// class G4GDMLWriteMaterials
//
// Class description:
//
// GDML class for writing materials definitions.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLWRITEMATERIALS_INCLUDED_
#define _G4GDMLWRITEMATERIALS_INCLUDED_

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include "G4GDMLWriteDefine.hh"

class G4GDMLWriteMaterials : public G4GDMLWriteDefine
{

 protected:

   void AddIsotope(const G4Isotope* const);
   void AddElement(const G4Element* const);
   void AddMaterial(const G4Material* const);

 private:

   void AtomWrite(xercesc::DOMElement*,const G4double&);
   void DWrite(xercesc::DOMElement*,const G4double&);
   void PWrite(xercesc::DOMElement*,const G4double&);
   void TWrite(xercesc::DOMElement*,const G4double&);
   void IsotopeWrite(const G4Isotope* const);
   void ElementWrite(const G4Element* const);
   void MaterialWrite(const G4Material* const);
   void MaterialsWrite(xercesc::DOMElement*);

 private:

   std::vector<const G4Isotope*> isotopeList;
   std::vector<const G4Element*> elementList;
   std::vector<const G4Material*> materialList;
   xercesc::DOMElement* materialsElement;
};

#endif
