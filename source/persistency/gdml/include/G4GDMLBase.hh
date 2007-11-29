#ifndef _G4GDMLBASE_INCLUDED_
#define _G4GDMLBASE_INCLUDED_

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/dom/DOM.hpp>

#include "G4GDMLEvaluator.hh"

class G4GDMLBase {
private:
   xercesc::XercesDOMParser* parser;
protected:
   G4GDMLEvaluator* evaluator;
   G4String prename;

   std::string nameProcess(const std::string&);
public:
   G4GDMLBase();
   ~G4GDMLBase();
   
   void gdmlRead(const G4String&);

   virtual void defineRead(const xercesc::DOMElement* const)=0;
   virtual void materialsRead(const xercesc::DOMElement* const)=0;
   virtual void solidsRead(const xercesc::DOMElement* const)=0;
   virtual void structureRead(const xercesc::DOMElement* const)=0;
   virtual void setupRead(const xercesc::DOMElement* const)=0;
};

#endif
