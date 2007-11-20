//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GDMLEvaluator.hh,v 1.5 2007-11-20 09:31:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLEvaluator
//
// Class description:
//
// GDML class for evaluation of expressions.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLEVALUATOR_INCLUDED_
#define _G4GDMLEVALUATOR_INCLUDED_

#include <CLHEP/Evaluator/Evaluator.h>

#include "G4String.hh"
#include "G4Types.hh"

class G4GDMLEvaluator {
   HepTool::Evaluator eval;
public:
   G4GDMLEvaluator();

   void Set(const G4GDMLEvaluator&);

   G4bool RegisterConstant(const G4String& name,G4double value);
   G4bool Evaluate(G4double& value,const G4String& expression,const G4String& unit="");

// If the expression is an empty string, the value of the expression is considered as zero
// The unit is an empty string by default, what means no unit or unit of one

// Do NOT change the default unit into "1" because a string is empty by default, so that if an empty
// string is passed as unit it means no unit or unit of one
};

#endif
