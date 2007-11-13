#include "G4GDMLParser.hh"

G4GDMLParser::G4GDMLParser() {

   try {

      xercesc::XMLPlatformUtils::Initialize();
   }
   catch(xercesc::XMLException& e) {

      char* message = xercesc::XMLString::transcode(e.getMessage());
      std::cerr << "XML toolkit initialization error: " << message << std::endl;
      xercesc::XMLString::release(&message);
   }

   parser = new xercesc::XercesDOMParser;

   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);
   parser->setValidationSchemaFullChecking(true);
}

G4GDMLParser::~G4GDMLParser() {

   if (parser) delete parser;

   try {

      xercesc::XMLPlatformUtils::Terminate();
   }
   catch(xercesc::XMLException& e) {
    
      char* message = xercesc::XMLString::transcode(e.getMessage());
      std::cerr << "XML toolkit termination error: " << message << std::endl;
      xercesc::XMLString::release(&message);
   }
}

bool G4GDMLParser::Read(const G4String& fileName) {

     return structure.gdmlRead(fileName,parser);
}

G4VPhysicalVolume *G4GDMLParser::GetWorldVolume(const G4String &setupName) {

   return structure.setup.Get(setupName);
}
