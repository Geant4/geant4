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
// $Id: G4GDMLStructure.hh,v 1.10 2007-11-21 13:23:53 ztorzsok Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLStructure
//
// Class description:
//
// GDML class for loading physical volumes according to various
// specifications in Geant4.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLSTRUCTURE_INCLUDED_
#define _G4GDMLSTRUCTURE_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SolidStore.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GDMLEvaluator.hh"
#include "G4GDMLMaterials.hh"
#include "G4GDMLSolids.hh"
#include "G4GDMLSetup.hh"

class G4GDMLStructure {

   G4String file,loop;

   G4GDMLEvaluator* evaluator;

   xercesc::XercesDOMParser* parser;

   G4bool directionRead           (const xercesc::DOMElement* const,EAxis&);
   G4bool divisionvolRead         (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool fileRead                (const xercesc::DOMElement* const,G4LogicalVolume**);
   G4bool loopRead                (const xercesc::DOMElement* const);
   G4bool paramvolRead            (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool physvolRead             (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool positionRead            (const xercesc::DOMElement* const,G4ThreeVector&);
   G4bool quantityRead            (const xercesc::DOMElement* const,G4double&);
   G4bool refRead                 (const xercesc::DOMElement* const,G4String&);
   G4bool replicate_along_axisRead(const xercesc::DOMElement* const,G4double&,G4double&,EAxis&);
   G4bool replicavolRead          (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool volumeRead              (const xercesc::DOMElement* const);
   G4bool volume_contentRead      (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool volume_content_loopRead (const xercesc::DOMElement* const,G4LogicalVolume*);
   G4bool Read                    (const xercesc::DOMElement* const,const G4String&,const G4String&);
public:
   G4GDMLMaterials materials;
   G4GDMLSolids solids;
   G4GDMLSetup setup;

   G4GDMLStructure();

   G4bool gdmlRead(const G4String&,xercesc::XercesDOMParser*);
   G4LogicalVolume* getVolume(const G4String&) const;
};

#endif
