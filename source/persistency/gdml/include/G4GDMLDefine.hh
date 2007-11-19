#ifndef _G4GDMLDEFINE_INCLUDED_
#define _G4GDMLDEFINE_INCLUDED_

#include <xercesc/dom/DOM.hpp>

#include <map>

#include "G4ThreeVector.hh"

#include "G4GDMLEvaluator.hh"

class G4GDMLDefine {

   G4String prename;

   G4GDMLEvaluator* evaluator;

   std::map<G4String,G4ThreeVector*> positionMap;
   std::map<G4String,G4ThreeVector*> rotationMap;

   bool constantRead(const xercesc::DOMElement* const);
   bool positionRead(const xercesc::DOMElement* const);
   bool rotationRead(const xercesc::DOMElement* const);

public:
   G4GDMLDefine();
   ~G4GDMLDefine();

   bool Read(const xercesc::DOMElement* const element,G4GDMLEvaluator*,const G4String&);

   G4ThreeVector *getPosition(const G4String&);
   G4ThreeVector *getRotation(const G4String&);
};

#endif
