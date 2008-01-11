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
// $Id: G4GDMLMaterials.hh,v 1.15 2008-01-11 12:21:21 ztorzsok Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLMaterials
//
// Class description:
//
// GDML class for loading isotopes, elements and materials according to
// specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLMATERIALS_INCLUDED_
#define _G4GDMLMATERIALS_INCLUDED_

#include "G4GDMLReadDefine.hh"

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

class G4GDMLMaterials : public G4GDMLReadDefine {
   G4double atomRead(const xercesc::DOMElement* const);
   G4int compositeRead(const xercesc::DOMElement* const,G4String&);
   G4double DRead(const xercesc::DOMElement* const);
   void elementRead(const xercesc::DOMElement* const);
   G4double fractionRead(const xercesc::DOMElement* const,G4String&);
   void isotopeRead(const xercesc::DOMElement* const);
   void materialRead(const xercesc::DOMElement* const);
   void mixtureRead(const xercesc::DOMElement* const,G4Element*);
   void mixtureRead(const xercesc::DOMElement* const,G4Material*);
   void opticalsurfaceRead(const xercesc::DOMElement* const);
   void materialsRead(const xercesc::DOMElement* const);
   G4Element* getElement(const G4String&) const;
   G4Isotope* getIsotope(const G4String&) const;
protected:
   G4Material* getMaterial(const G4String&) const;
};

#endif
