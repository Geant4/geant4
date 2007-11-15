#ifndef _G4GDMLPARSER_INCLUDED_
#define _G4GDMLPARSER_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>

#include "G4GDMLStructure.hh"

class G4GDMLParser {
   xercesc::XercesDOMParser *parser;
   G4GDMLStructure structure;
public:
   G4GDMLParser();
   ~G4GDMLParser();

   bool Read(const G4String &fileName);

   G4VPhysicalVolume* GetWorldVolume(const G4String& setupName="Default");
};

#endif
