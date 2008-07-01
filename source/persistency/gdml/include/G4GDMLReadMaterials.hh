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
// $Id: G4GDMLReadMaterials.hh,v 1.6 2008-07-01 08:12:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLReadMaterials
//
// Class description:
//
// GDML class for loading isotopes, elements and materials according to
// specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLREADMATERIALS_INCLUDED_
#define _G4GDMLREADMATERIALS_INCLUDED_

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GDMLReadDefine.hh"

class G4GDMLReadMaterials : public G4GDMLReadDefine {
   G4double atomRead(const xercesc::DOMElement* const);
   G4int compositeRead(const xercesc::DOMElement* const,G4String&);
   G4double DRead(const xercesc::DOMElement* const);
   G4double PRead(const xercesc::DOMElement* const);
   G4double TRead(const xercesc::DOMElement* const);
   void elementRead(const xercesc::DOMElement* const);
   G4double fractionRead(const xercesc::DOMElement* const,G4String&);
   void isotopeRead(const xercesc::DOMElement* const);
   void materialRead(const xercesc::DOMElement* const);
   void mixtureRead(const xercesc::DOMElement* const,G4Element*);
   void mixtureRead(const xercesc::DOMElement* const,G4Material*);
   void propertyRead(const xercesc::DOMElement* const,G4Material*);
   void materialsRead(const xercesc::DOMElement* const);
protected:
   G4Element* getElement(const G4String&,bool verbose=true) const;
   G4Isotope* getIsotope(const G4String&,bool verbose=true) const;
   G4Material* getMaterial(const G4String&,bool verbose=true) const;
};

#endif
