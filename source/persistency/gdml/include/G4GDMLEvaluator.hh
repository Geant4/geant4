#ifndef _G4GDMLEVALUATOR_INCLUDED_
#define _G4GDMLEVALUATOR_INCLUDED_

#include "CLHEP/Evaluator/Evaluator.h"

#include <iostream>
#include <string>

class G4GDMLEvaluator {
   HepTool::Evaluator eval;
public:
   G4GDMLEvaluator();
   ~G4GDMLEvaluator();

   static G4GDMLEvaluator *GetInstance();

   bool RegisterConstant(const std::string &name,double value);

   bool Evaluate(double& value,const std::string& expression,const std::string& unit="");

// If the expression is an empty string, the value of the expression is considered as zero
// The unit is an empty string by default, what means no unit or unit of one

// Do NOT change the default unit into "1" because a string is empty by default, so that if an empty
// string is passed as unit it means no unit or unit of one

};

#endif
