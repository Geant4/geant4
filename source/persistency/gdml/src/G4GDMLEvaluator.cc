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
// $Id: G4GDMLEvaluator.cc,v 1.9 2007-11-26 14:31:32 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLEvaluator Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLEvaluator.hh"

G4GDMLEvaluator::G4GDMLEvaluator() {

   eval.clear();
   eval.setStdMath();
   eval.setSystemOfUnits(1.e+3,1./1.60217733e-25,1.e+9,1./1.60217733e-10,1.0,1.0,1.0);
}

void G4GDMLEvaluator::Set(const G4GDMLEvaluator &right) {

   memcpy(this,&right,sizeof(G4GDMLEvaluator));
}

void G4GDMLEvaluator::defineVariable(const G4String& name,G4double value) {

   if (eval.findVariable(name)) G4Exception("GDML: Variable '"+name+"' is already defined!");

   eval.setVariable(name.c_str(),value);
}

void G4GDMLEvaluator::setVariable(const G4String& name,G4double value) {

   if (!eval.findVariable(name)) G4Exception("GDML: Variable '"+name+"' was not found!");

   eval.setVariable(name.c_str(),value);
}

G4double G4GDMLEvaluator::Evaluate(const G4String& expression) {

   G4double value = 0.0;

   if (!expression.empty()) {
   
      value = eval.evaluate(expression.c_str());

      if (eval.status() != HepTool::Evaluator::OK) {

         eval.print_error();
         G4Exception("Error in evaluator!");
      }
   }
   
   return value;
}

