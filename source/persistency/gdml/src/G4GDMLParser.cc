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

bool G4GDMLParser::Read(const std::string &fileName) {

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

   xercesc::DOMNodeList* children = element->getChildNodes();
   const  XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue; 
      
      element = dynamic_cast<xercesc::DOMElement*>(node);   
      
      std::string tag = xercesc::XMLString::transcode(element->getTagName());

      if (tag=="define"   ) { if (!structure.solids.define.Read(element)) return false; } else
      if (tag=="materials") { if (!structure.materials.Read(element)    ) return false; } else
      if (tag=="solids"   ) { if (!structure.solids.Read(element)       ) return false; } else
      if (tag=="structure") { if (!structure.Read(element)              ) return false; } else
      if (tag=="setup"    ) { if (!setup.Read(element)                  ) return false; }
   }

   return true;
}

G4VPhysicalVolume *G4GDMLParser::GetWorldVolume(const std::string &setupName) {

   return setup.Get(setupName);
}
