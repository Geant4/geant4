#ifndef _G4GDMLMATERIALS_INCLUDED_
#define _G4GDMLMATERIALS_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include "G4Material.hh"

#include "G4GDMLEvaluator.hh"

class G4GDMLMaterials {

   G4GDMLEvaluator* evaluator;

   bool compositeRead(const xercesc::DOMElement* const,G4Material*);
   bool elementRead(const xercesc::DOMElement* const);
   bool fractionRead(const xercesc::DOMElement* const,double&,std::string&);
   bool materialRead(const xercesc::DOMElement* const);
   bool valueRead(const xercesc::DOMElement* const,double&);
public:
   G4GDMLMaterials();

   bool Read(const xercesc::DOMElement* const);
};

#endif
