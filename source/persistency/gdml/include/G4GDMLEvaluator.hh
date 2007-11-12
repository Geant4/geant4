#ifndef _G4GDMLEVALUATOR_INCLUDED_
#define _G4GDMLEVALUATOR_INCLUDED_

#include "CLHEP/Evaluator/Evaluator.h"

#include "G4String.hh"
#include "G4Types.hh"

class G4GDMLEvaluator {
   HepTool::Evaluator eval;
public:
   G4GDMLEvaluator();
   ~G4GDMLEvaluator();

   static G4GDMLEvaluator *GetInstance();

   bool RegisterConstant(const G4String& name,G4double value);

   bool Evaluate(G4double& value,const G4String& expression,const G4String& unit="");

// If the expression is an empty string, the value of the expression is considered as zero
// The unit is an empty string by default, what means no unit or unit of one

// Do NOT change the default unit into "1" because a string is empty by default, so that if an empty
// string is passed as unit it means no unit or unit of one

};

#endif
