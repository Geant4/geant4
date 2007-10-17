#ifndef _G4GDMLSETUP_INCLUDED_
#define _G4GDMLSETUP_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include <iostream>
#include <string>
#include <map>

#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

class G4GDMLSetup {
   std::map<std::string,std::string> setupMap;
public:
   G4VPhysicalVolume *Get(const std::string&);
   bool Read(const xercesc::DOMElement* const);
};

#endif
