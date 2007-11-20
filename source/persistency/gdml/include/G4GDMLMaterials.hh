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
// $Id: G4GDMLMaterials.hh,v 1.6 2007-11-20 09:31:44 gcosmo Exp $
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

#include <xercesc/dom/DOM.hpp>

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include "G4GDMLEvaluator.hh"

class G4GDMLMaterials {

   G4String prename;

   G4GDMLEvaluator* evaluator;

   G4bool atomRead     (const xercesc::DOMElement* const,G4double& _value);
   G4bool DRead        (const xercesc::DOMElement* const,G4double& _value);
   G4bool elementRead  (const xercesc::DOMElement* const);
   G4bool fractionRead (const xercesc::DOMElement* const,G4double& _n,G4String& ref);
   G4bool isotopeRead  (const xercesc::DOMElement* const);
   G4bool materialRead (const xercesc::DOMElement* const);
   G4bool mixtureRead  (const xercesc::DOMElement* const,G4Element*);
   G4bool mixtureRead  (const xercesc::DOMElement* const,G4Material*);
public:
   G4bool Read(const xercesc::DOMElement* const element,G4GDMLEvaluator*,const G4String&);
   G4Material* getMaterial(const G4String&) const;
};

#endif
