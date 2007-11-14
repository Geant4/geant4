#ifndef _G4GDMLSETUP_INCLUDED_
#define _G4GDMLSETUP_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include <map>

#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

class G4GDMLSetup {
   
   G4String module;

   std::map<G4String,G4String> setupMap;
public:
   G4VPhysicalVolume *Get(const G4String&);
   G4String GetS(const G4String&);

   bool Read(const xercesc::DOMElement* const element,const G4String& newModule);
};

#endif
