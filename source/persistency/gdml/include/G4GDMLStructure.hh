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

   G4String module;

   G4GDMLEvaluator* evaluator;

   xercesc::XercesDOMParser* parser;

   bool directionRead           (const xercesc::DOMElement* const,EAxis&);
   bool divisionvolRead         (const xercesc::DOMElement* const,G4LogicalVolume*);
   bool fileRead                (const xercesc::DOMElement* const,G4String&);
   bool physvolRead             (const xercesc::DOMElement* const,G4LogicalVolume*);
   bool quantityRead            (const xercesc::DOMElement* const,G4double&);
   bool refRead                 (const xercesc::DOMElement* const,G4String&);
   bool replicate_along_axisRead(const xercesc::DOMElement* const,G4double&,G4double&,EAxis&);
   bool replicavolRead          (const xercesc::DOMElement* const,G4LogicalVolume*);
   bool volumeRead              (const xercesc::DOMElement* const);
   bool Read(const xercesc::DOMElement* const element,const G4String& newModule);
public:
   G4GDMLMaterials materials;
   G4GDMLSolids solids;
   G4GDMLSetup setup;

   G4GDMLStructure();

   bool gdmlRead(const G4String&,xercesc::XercesDOMParser*);
};

#endif
