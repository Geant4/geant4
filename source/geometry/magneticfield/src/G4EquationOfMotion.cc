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
// $Id: G4EquationOfMotion.cc,v 1.5 2001-11-09 20:39:07 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4EquationOfMotion.hh"

static const int G4maximum_number_of_field_components = 16;

G4EquationOfMotion::~G4EquationOfMotion()
{}

void 
G4EquationOfMotion::RightHandSide( const  G4double y[],
				   G4double dydx[]   ) const
{
     G4double Field[G4maximum_number_of_field_components];   
     G4double  PositionAndTime[4];

     // Position
     PositionAndTime[0] = y[0];
     PositionAndTime[1] = y[1];
     PositionAndTime[2] = y[2];
     // Global Time
     PositionAndTime[3] = y[7];  // See G4FieldTrack::LoadFromArray

     GetFieldValue(PositionAndTime, Field) ;
     EvaluateRhsGivenB( y, Field, dydx );
}

void 
G4EquationOfMotion::EvaluateRhsReturnB( const G4double y[],
				 G4double dydx[],
				 G4double  Field[]  ) const
{
     G4double  PositionAndTime[4];

     // Position
     PositionAndTime[0] = y[0];
     PositionAndTime[1] = y[1];
     PositionAndTime[2] = y[2];
     // Global Time
     PositionAndTime[3] = y[7];  // See G4FieldTrack::LoadFromArray

     GetFieldValue(PositionAndTime, Field) ;
     EvaluateRhsGivenB( y, Field, dydx );
}

#if  HELP_THE_COMPILER
void 
G4EquationOfMotion::doNothing()
{
}
#endif
