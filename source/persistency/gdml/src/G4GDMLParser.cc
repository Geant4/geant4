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

   Parser = new xercesc::XercesDOMParser;

   Parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   Parser->setDoNamespaces(true);
   Parser->setDoSchema(true);
   Parser->setValidationSchemaFullChecking(true);
}

G4GDMLParser::~G4GDMLParser() {

   if (Parser) delete Parser;

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

   try {

      Parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      std::cerr << "XML: " << message << std::endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      std::cerr << "DOM: " << message << std::endl;
      xercesc::XMLString::release(&message);
   }
   catch (...) {

      std::cerr << "Unexpected exception!" << std::endl;
   }

   xercesc::DOMDocument* doc = Parser->getDocument();

   if (!doc) {
   
      std::cout << "Unable to open document '" << fileName << "'!" << std::endl;
      return false;
   }

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) {
   
      std::cout << "Empty document!" << std::endl;
      return false;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) { if (!structure.solids.define.Read(child)) return false; } else
      if (tag=="materials") { if (!structure.materials.Read(child)    ) return false; } else
      if (tag=="solids"   ) { if (!structure.solids.Read(child)       ) return false; } else
      if (tag=="structure") { if (!structure.Read(child)              ) return false; } else
      if (tag=="setup"    ) { if (!setup.Read(child)                  ) return false; }
   }

   return true;
}

G4VPhysicalVolume *G4GDMLParser::GetWorldVolume(const G4String &setupName) {

   return setup.Get(setupName);
}
