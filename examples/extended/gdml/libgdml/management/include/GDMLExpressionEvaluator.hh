//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: GDMLExpressionEvaluator.hh,v 1.2 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
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

