#include "G4GDMLEvaluator.hh"

G4GDMLEvaluator::G4GDMLEvaluator() {

   eval = new HepTool::Evaluator[8];
   max_eval = 8;
   index = 0;

   Init();
}

G4GDMLEvaluator::~G4GDMLEvaluator() {

   if (eval) delete [] eval;
}

G4GDMLEvaluator *G4GDMLEvaluator::GetInstance() {

   static G4GDMLEvaluator instance;

   return &instance;
}

bool G4GDMLEvaluator::RegisterConstant(const G4String &name,G4double value) {

   eval[index].setVariable(name.c_str(),value);

   return true;
}

bool G4GDMLEvaluator::Evaluate(G4double& value,const G4String& expression,const G4String& unit) {

   G4double _expression = 0.0;
   G4double _unit = 1.0;

   if (expression != "") {
   
      _expression = eval[index].evaluate(expression.c_str());

      if (eval[index].status() != HepTool::Evaluator::OK) {

         eval[index].print_error();
         return false;
      }
   }

   if (unit != "") {
   
      _unit = eval[index].evaluate(unit.c_str());

      if (eval[index].status() != HepTool::Evaluator::OK) {

         eval[index].print_error();
         return false;
      }
   }

   value = _expression*_unit;

   return true;
}

void G4GDMLEvaluator::Push() {

   if (index >= max_eval) return;

   index++;
   Init();
}

void G4GDMLEvaluator::Pop() {

   if (index <= 0) return;

   index--;
}

void G4GDMLEvaluator::Init() {

   eval[index].clear();
   eval[index].setStdMath();
   eval[index].setSystemOfUnits(1.e+3,1./1.60217733e-25,1.e+9,1./1.60217733e-10,1.0,1.0,1.0);
}
