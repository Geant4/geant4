#ifndef _G4GDMLMATERIALS_INCLUDED_
#define _G4GDMLMATERIALS_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include "G4GDMLEvaluator.hh"

class G4GDMLMaterials {

   G4String prename;

   G4GDMLEvaluator* evaluator;

   bool atomRead     (const xercesc::DOMElement* const,double& _value);
   bool DRead        (const xercesc::DOMElement* const,double& _value);
   bool elementRead  (const xercesc::DOMElement* const);
   bool fractionRead (const xercesc::DOMElement* const,double& _n,std::string& ref);
   bool isotopeRead  (const xercesc::DOMElement* const);
   bool materialRead (const xercesc::DOMElement* const);
   bool mixtureRead  (const xercesc::DOMElement* const,G4Element*);
   bool mixtureRead  (const xercesc::DOMElement* const,G4Material*);
public:
   bool Read(const xercesc::DOMElement* const element,G4GDMLEvaluator*,const G4String&);
   G4Material* getMaterial(const G4String&) const;
};

#endif
