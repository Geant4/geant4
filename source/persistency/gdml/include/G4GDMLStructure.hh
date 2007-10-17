#ifndef _G4GDMLSTRUCTURE_INCLUDED_
#define _G4GDMLSTRUCTURE_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include <iostream>
#include <string>

#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SolidStore.hh"

#include "G4GDMLMaterials.hh"
#include "G4GDMLSolids.hh"

class G4GDMLStructure {

   bool physvolRead(const xercesc::DOMElement* const,G4LogicalVolume*);
   bool volumeRead(const xercesc::DOMElement* const);
public:
   G4GDMLMaterials materials;
   G4GDMLSolids solids;

   bool Read(const xercesc::DOMElement* const);
};

#endif
