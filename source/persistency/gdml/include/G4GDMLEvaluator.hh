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
// $Id: G4GDMLEvaluator.hh,v 1.15 2008-07-11 07:50:07 gcosmo Exp $
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
#include <vector>

#include "globals.hh"

class G4GDMLEvaluator
{
   HepTool::Evaluator eval;
   std::vector<G4String> variableList;

 public:

   G4GDMLEvaluator();

   void defineConstant(const G4String&, G4double);
   void defineVariable(const G4String&, G4double);
   void defineMatrix(const G4String&, G4int, std::vector<G4double>);
   void setVariable(const G4String&, G4double);
   G4bool isVariable(const G4String&) const;
   G4String SolveBrackets(const G4String&);
   G4double Evaluate(const G4String&);
   G4int EvaluateInteger(const G4String&);
   G4double getConstant(const G4String&);
   G4double getVariable(const G4String&);
};

#endif
