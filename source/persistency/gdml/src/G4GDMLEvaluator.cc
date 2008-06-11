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
// $Id: G4GDMLEvaluator.cc,v 1.16 2008-06-11 09:36:13 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLEvaluator Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include <sstream>
#include "G4GDMLEvaluator.hh"

G4GDMLEvaluator::G4GDMLEvaluator() {

   eval.clear();
   eval.setStdMath();
   eval.setSystemOfUnits(1.e+3,1./1.60217733e-25,1.e+9,1./1.60217733e-10,1.0,1.0,1.0);
}

void G4GDMLEvaluator::defineConstant(const G4String& name,G4double value) {

   if (eval.findVariable(name)) G4Exception("G4GDML: ERROR! Redefinition of constant or variable: "+name);
   eval.setVariable(name.c_str(),value);
}

void G4GDMLEvaluator::defineVariable(const G4String& name,G4double value) {

   if (eval.findVariable(name)) G4Exception("G4GDML: ERROR! Redefinition of constant or variable: "+name);
   eval.setVariable(name.c_str(),value);
   variableList.push_back(name);
}

void G4GDMLEvaluator::defineMatrix(const G4String& name,G4int coldim,std::vector<G4double> valueList) {

   const G4int size = valueList.size();

   if (size == 0) G4Exception("G4GDML: ERROR! Matrix '"+name+"' is empty!");
   if (size == 1) G4Exception("G4GDML: ERROR! Matrix '"+name+"' has only one element! Define a constant instead!");
   if (size % coldim != 0) G4Exception("G4GDML: ERROR! Matrix '"+name+"' is not filled correctly!");

   if ((size == coldim) || (coldim == 1)) { // Row- or column matrix
   
      for (G4int i=0;i<size;i++) {
   
         std::stringstream MatrixElementNameStream;
         MatrixElementNameStream << name << "_" << i;
         defineConstant(MatrixElementNameStream.str(),valueList[i]);
      }
   } else { // Normal matrix
   
      const G4int rowdim = size/coldim;

      for (G4int i=0;i<rowdim;i++)
      for (G4int j=0;j<coldim;j++) {
      
         std::stringstream MatrixElementNameStream;
         MatrixElementNameStream << name << "_" << i << "_" << j;
         defineConstant(MatrixElementNameStream.str(),valueList[coldim*i+j]);
      }
   }
}

void G4GDMLEvaluator::setVariable(const G4String& name,G4double value) {

   if (!isVariable(name)) G4Exception("G4GDML: ERROR! Variable '"+name+"' is not defined!");
   eval.setVariable(name.c_str(),value);
}

bool G4GDMLEvaluator::isVariable(const G4String& name) const {

   const size_t variableCount = variableList.size();

   for (size_t i=0;i<variableCount;i++) {

      if (variableList[i] == name) return true;
   }

   return false;
}

G4String G4GDMLEvaluator::SolveBrackets(const G4String& in) {

   const std::string::size_type open = in.find("[",0);
   const std::string::size_type close = in.find("]",0);

   if (open==close) return in;

   if (open>close) G4Exception("G4GDML: ERROR! Bracket mismatch: "+in);

   std::string::size_type begin = open;
   std::string::size_type end = 0;

   std::string out;
   out.append(in,0,open);
   
   do {

      end = in.find(",",begin+1);
      if (end==std::string::npos) end = close;

      std::stringstream indexStream;
      indexStream << "_" << EvaluateInteger(in.substr(begin+1,end-begin-1));

      out.append(indexStream.str());

      begin = end;

   } while (end<close);

   return out;
}

G4double G4GDMLEvaluator::Evaluate(const G4String& in) {

   G4String expression = SolveBrackets(in);

   G4double value = 0.0;

   if (!expression.empty()) {
   
      value = eval.evaluate(expression.c_str());

      if (eval.status() != HepTool::Evaluator::OK) {

         eval.print_error();
         G4Exception("G4GDML: Error in expression: "+expression);
      }
   }
   
   return value;
}

G4int G4GDMLEvaluator::EvaluateInteger(const G4String& expression) {

// This function is for evaluating integer expressions, like loop variables and matrix indices
// Complains if the evaluated expression has a fractional part different from zero

   G4double value = Evaluate(expression);

   G4int whole = (G4int)value;
   G4double frac = value - (G4double)whole;

   if (frac != 0.0) G4Exception("G4GDML: ERROR! Expression '"+expression+"' is expected to have an integer value!");

   return whole;
}

G4double G4GDMLEvaluator::getConstant(const G4String& name) {

   if (isVariable(name)) G4Exception("G4GDML: ERROR! Constant '"+name+"' is not defined! It is a variable!");
   if (!eval.findVariable(name)) G4Exception("G4GDML: ERROR! Constant '"+name+"' is not defined!"); 
   return Evaluate(name);
}

G4double G4GDMLEvaluator::getVariable(const G4String& name) {

   if (!isVariable(name)) G4Exception("G4GDML: ERROR! Variable '"+name+"' is not a defined!");
   return Evaluate(name);
}
