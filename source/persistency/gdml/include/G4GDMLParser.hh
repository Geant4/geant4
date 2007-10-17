#ifndef _G4GDMLPARSER_INCLUDED_
#define _G4GDMLPARSER_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>

#include "G4GDMLStructure.hh"
#include "G4GDMLSetup.hh"

class G4GDMLParser {

   xercesc::XercesDOMParser *Parser;

   G4GDMLStructure structure;
   G4GDMLSetup setup;

public:
   G4GDMLParser();
   ~G4GDMLParser();

   bool Read(const std::string &fileName);

   G4VPhysicalVolume *GetWorldVolume(const std::string &setupName="Default");
};

#endif
