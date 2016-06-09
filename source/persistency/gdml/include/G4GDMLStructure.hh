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
// $Id: G4GDMLStructure.hh,v 1.16 2007/11/30 14:51:20 ztorzsok Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SolidStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ReflectionFactory.hh"

#include "G4GDMLSetup.hh"

class G4GDMLStructure : public G4GDMLSetup {
private:
   EAxis directionRead(const xercesc::DOMElement* const);
   void divisionvolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   G4LogicalVolume* fileRead(const xercesc::DOMElement* const);
   void loopRead(const xercesc::DOMElement* const);
   void paramvolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   void physvolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   G4double quantityRead(const xercesc::DOMElement* const);
   G4String refRead(const xercesc::DOMElement* const);
   void replicate_along_axisRead(const xercesc::DOMElement* const,G4double&,G4double&,EAxis&);
   void replicavolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   void volumeRead(const xercesc::DOMElement* const);
   void volume_contentRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   void volume_loopRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   void structureRead(const xercesc::DOMElement* const);
protected:
   G4LogicalVolume* getVolume(const G4String&) const;
};

#endif
