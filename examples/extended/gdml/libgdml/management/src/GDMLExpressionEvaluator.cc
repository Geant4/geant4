#include "GDMLExpressionEvaluator.hh"
#include "CLHEP/Units/PhysicalConstants.h"

GDMLExpressionEvaluator::GDMLExpressionEvaluator()
{
  fCalc.clear();
  fCalc.setStdMath();                 // set standard constants and functions
  //fCalc.setSystemOfUnits();           // set SI units
  // Set Geant4 system of units
  fCalc.setSystemOfUnits(1.e+3, 1./1.60217733e-25, 1.e+9, 1./1.60217733e-10,1.0, 1.0, 1.0);  
}

GDMLExpressionEvaluator::~GDMLExpressionEvaluator()
{
  fCTable.clear();
  fPCTable.clear();
  fCalc.clear();
}

//StatusCode GDMLExpressionEvaluator::RegisterConstant( const GDMLConstant* const c )
StatusCode GDMLExpressionEvaluator::RegisterConstant( const define::constant* const c )
{
  double value = fCalc.evaluate( c->get_value().c_str() );
  
  if( fCalc.status() != HepTool::Evaluator::OK )
  {
    std::cerr << "Expression evaluator:: Error registering constant " << c->get_name() << std::endl;
    fCalc.print_error();
    std::cout << std::endl;
    return StatusCode::eError;
  }
  
  std::cout << "Expression evaluator:: Registering constant "
            << c->get_name() << ": " << value << std::endl;
  
  //fCalc.setVariable( c->get_name().c_str(), c->get_value().c_str() );
  fCalc.setVariable( c->get_name().c_str(), value );
  return StatusCode::eOk;
}

//StatusCode GDMLExpressionEvaluator::RegisterPhysConstant( const GDMLPhysicalConstant* const physc )
StatusCode GDMLExpressionEvaluator::RegisterPhysConstant( const define::quantity* const physc )
{
  std::string expr = physc->get_value();
  expr += "*(";
  expr += physc->get_unit();
  expr += ")";
  
  //std::cout << "Expression evaluator:: evaluating string: " << expr << std::endl;
  
  double value      = fCalc.evaluate( expr.c_str() );
  double unit_value = fCalc.evaluate( physc->get_unit().c_str() );
  
  if( fCalc.status() != HepTool::Evaluator::OK )
  {
    std::cerr << "Expression evaluator:: Error registering quantity "
              << physc->get_name() << std::endl;
    fCalc.print_error();
    std::cout << std::endl;
    return StatusCode::eError;
  }
  
  std::cout << "Expression evaluator:: Registering quantity "
            << physc->get_name() << ": " << (value/unit_value) << physc->get_unit() << std::endl;
  
  //fCalc.setVariable( physc->get_name().c_str(), expr.c_str() );
  fCalc.setVariable( physc->get_name().c_str(), value );
  return StatusCode::eOk;
}

//StatusCode GDMLExpressionEvaluator::RegisterExpression( const GDMLExpression* e )
StatusCode GDMLExpressionEvaluator::RegisterExpression( const define::expression* e )
{
  std::string expr = "(";
  expr += e->get_text();
  expr += ")";
  double value = fCalc.evaluate( expr.c_str() );
  
  if( fCalc.status() != HepTool::Evaluator::OK )
  {
    std::cerr << "Expression evaluator:: Error registering expression " << e->get_name() << std::endl;
    fCalc.print_error();
    std::cout << std::endl;
    return StatusCode::eError;
  }
  
  std::cout << "Expression evaluator:: Registering expression "
            << e->get_name() << ": " << value << std::endl;
  
  //fCalc.setVariable( e->get_name().c_str(), e->get_text().c_str() );
  fCalc.setVariable( e->get_name().c_str(), value );
  return StatusCode::eOk;
}

double GDMLExpressionEvaluator::Eval( const std::string& expr )
{
  return Eval( expr.c_str() );
}

double GDMLExpressionEvaluator::Eval( const char* expr )
{
  double result = fCalc.evaluate( expr );
  if( fCalc.status() != HepTool::Evaluator::OK )
  {
    std::cerr << expr << std::endl;
    //std::cerr << "------";
    for (int i=0; i<fCalc.error_position(); i++)
    {
      std::cerr << "-";
    }
    std::cerr << "^\a" << std::endl;
    fCalc.print_error();
    std::cerr << std::endl;
  }
  return result;
}


