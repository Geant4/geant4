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
// $Id: GDMLDefineTable.hh,v 1.2 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
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


