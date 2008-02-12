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
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLREADSTRUCTURE_INCLUDED_
#define _G4GDMLREADSTRUCTURE_INCLUDED_

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SolidStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ReflectionFactory.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4GDMLReadParamvol.hh"

class G4GDMLReadStructure : public G4GDMLReadParamvol {
public:
   typedef std::pair<G4String,G4String> AuxPairType;
   typedef std::vector<AuxPairType> AuxListType;
   typedef std::map<const G4LogicalVolume* const,AuxListType> AuxMapType;
private:
   AuxMapType auxMap;

   G4LogicalVolume *pMotherLogical;

   AuxPairType auxiliaryRead(const xercesc::DOMElement* const);
   EAxis directionRead(const xercesc::DOMElement* const);
   void divisionvolRead(const xercesc::DOMElement* const);
   G4LogicalVolume* fileRead(const xercesc::DOMElement* const);
   void physvolRead(const xercesc::DOMElement* const);
   G4double quantityRead(const xercesc::DOMElement* const);
   void replicate_along_axisRead(const xercesc::DOMElement* const,G4double&,G4double&,EAxis&);
   void replicavolRead(const xercesc::DOMElement* const);
   void volumeRead(const xercesc::DOMElement* const);
   void volume_contentRead(const xercesc::DOMElement* const);
   void structureRead(const xercesc::DOMElement* const);
public:
   G4LogicalVolume* getVolume(const G4String&) const;
   AuxListType getVolumeAuxiliaryInformation(const G4LogicalVolume* const);
   const AuxMapType* getAuxiliaryMap();
};

#endif
