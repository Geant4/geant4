#ifndef GDML_EXPRESSION_EVALUATOR_H
#define GDML_EXPRESSION_EVALUATOR_H 1

#include "CLHEP/Evaluator/Evaluator.h"

#include "StatusCode.hh"
#include "GDMLDefineTable.hh"

#include "string"

class GDMLExpressionEvaluator
{
public:
  GDMLExpressionEvaluator();
  ~GDMLExpressionEvaluator();
  
//   StatusCode RegisterConstant( const GDMLConstant* const c );
//   StatusCode RegisterPhysConstant( const GDMLPhysicalConstant* const physc );
//   StatusCode RegisterExpression( const GDMLExpression* e );
  StatusCode RegisterConstant( const define::constant* const c );
  StatusCode RegisterPhysConstant( const define::quantity* const physc );
  StatusCode RegisterExpression( const define::expression* e );
  double Eval( const std::string& expr );
  double Eval( const char* expr );
  
private:
  HepTool::Evaluator     fCalc;
  ConstantsTable         fCTable;
  PhysicalConstantsTable fPCTable;
};

#endif // GDML_EXPRESSION_EVALUATOR_H

