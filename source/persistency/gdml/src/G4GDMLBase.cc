#include "G4GDMLBase.hh"

G4GDMLBase::G4GDMLBase() {

   try {

      xercesc::XMLPlatformUtils::Initialize();
   }
   catch(xercesc::XMLException& e) {

      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cerr << "XML toolkit initialization error: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   parser = new xercesc::XercesDOMParser;

   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);
   parser->setValidationSchemaFullChecking(true);

   evaluator = new G4GDMLEvaluator();
}

G4GDMLBase::~G4GDMLBase() {

   if (evaluator) delete evaluator;

   if (parser) delete parser;

   try {

      xercesc::XMLPlatformUtils::Terminate();
   }
   catch(xercesc::XMLException& e) {
    
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cerr << "XML toolkit termination error: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
}

std::string G4GDMLBase::nameProcess(const std::string& in) {

   std::string out(prename);
   
   std::string::size_type open = in.find("[",0);

   out.append(in,0,open);
   
   while (open != std::string::npos) {
   
      std::string::size_type close = in.find("]",open);

      if (close == std::string::npos) G4Exception("Bracket mismatch in loop!");
   
      std::string expr = in.substr(open+1,close-open-1);

      std::stringstream stream;
      
      stream << "[" << evaluator->EvaluateInteger(expr) << "]";
   
      out.append(stream.str());

      open = in.find("[",close);
   }

   return out;
}

void G4GDMLBase::gdmlRead(const G4String& fileName) {

   prename = fileName + "_";

   try {

      parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "XML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "DOM: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) G4Exception("GDML: Unable to open document: "+fileName);

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) G4Exception("GDML: Empty document!");

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) defineRead(child); else
      if (tag=="materials") materialsRead(child); else
      if (tag=="solids"   ) solidsRead(child); else
      if (tag=="setup"    ) setupRead(child); else
      if (tag=="structure") structureRead(child);
   }
}
