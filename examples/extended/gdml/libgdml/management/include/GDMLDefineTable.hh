#ifndef GDML_DEFINE_TABLE
#define GDML_DEFINE_TABLE 1

// #include "GDMLConstant.hh"
// #include "GDMLPhysicalConstant.hh"
// #include "GDMLExpression.hh"
// #include "GDMLCartesianVectorType.hh"

#include "define.hh"

#include <map>

// typedef std::map< std::string, GDMLConstant >         ConstantsTable;
// typedef std::map< std::string, GDMLPhysicalConstant > PhysicalConstantsTable;
// typedef std::map< std::string, GDMLExpression >       ExpressionsTable;
// typedef std::map< std::string, GDMLPosition >         PositionsTable;
// typedef std::map< std::string, GDMLRotation >         RotationsTable;
typedef std::map< std::string, define::constant >         ConstantsTable;
typedef std::map< std::string, define::quantity >         PhysicalConstantsTable;
typedef std::map< std::string, define::expression >       ExpressionsTable;
typedef std::map< std::string, define::position >         PositionsTable;
typedef std::map< std::string, define::rotation >         RotationsTable;

#endif // GDML_DEFINE_TABLE


