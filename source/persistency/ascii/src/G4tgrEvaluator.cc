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
// $Id: G4tgrEvaluator.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrEvaluator

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrEvaluator.hh"

#include <cmath>

// -------------------------------------------------------------------------
G4tgrEvaluator::G4tgrEvaluator()
{
  AddCommonFunctions();
}


// -------------------------------------------------------------------------
G4tgrEvaluator::~G4tgrEvaluator()
{
}


// -------------------------------------------------------------------------
void G4tgrEvaluator::print_error( G4int estatus ) const
{
  switch (estatus)
  {
    case ERROR_SYNTAX_ERROR:
      G4cerr << "G4tgrEvaluator: syntax error!" << G4endl;
      return;
    default:
      G4Evaluator::print_error();    
    return;
  }
} 
 
G4double fsin( G4double arg ){  return std::sin(arg); }
G4double fcos( G4double arg ){  return std::cos(arg); }
G4double ftan( G4double arg ){  return std::tan(arg); }
G4double fasin( G4double arg ){  return std::asin(arg); }
G4double facos( G4double arg ){  return std::acos(arg); }
G4double fatan( G4double arg ){  return std::atan(arg); }
G4double fatan2( G4double arg1, G4double arg2 ){ return std::atan2(arg1,arg2); }
G4double fsinh( G4double arg ){  return std::sinh(arg); }
G4double fcosh( G4double arg ){  return std::cosh(arg); }
G4double ftanh( G4double arg ){  return std::tanh(arg); }
// G4double fasinh( G4double arg ){  return std::asinh(arg); }
// G4double facosh( G4double arg ){  return std::acosh(arg); }
// G4double fatanh( G4double arg ){  return std::atanh(arg); }
G4double fsqrt( G4double arg ){  return std::sqrt(arg); }
G4double fexp( G4double arg ){  return std::exp(arg); }
G4double flog( G4double arg ){  return std::log(arg); }
G4double flog10( G4double arg ){  return std::log10(arg); }
G4double fpow( G4double arg1, G4double arg2 ){  return std::pow(arg1,arg2); }


//--------------------------------------------------------------------
void G4tgrEvaluator::AddCommonFunctions()
{
  setFunction("sin", (*fsin));
  setFunction("cos", (*fcos));
  setFunction("tan", (*ftan));
  setFunction("asin", (*fasin));
  setFunction("acos", (*facos));
  setFunction("atan", (*fatan));
  setFunction("atan2", (*fatan2));
  setFunction("sinh", (*fsinh));
  setFunction("cosh", (*fcosh));
  setFunction("tanh", (*ftanh));
//  setFunction("asinh", (*fasinh));
//  setFunction("acosh", (*facosh));
//  setFunction("atanh", (*fatanh));
  setFunction("sqrt", (*fsqrt));
  setFunction("exp", (*fexp));
  setFunction("log", (*flog));
  setFunction("log10", (*flog10));
  setFunction("pow", (*fpow));
}
