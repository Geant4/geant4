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
// $Id: G4Fcn.hh 72292 2013-07-15 12:01:43Z ihrivnac $

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisUtilities_h
#define G4AnalysisUtilities_h 1

#include "globals.hh"

#include <vector>
#include <memory>

// Enumeration for definition of available output types

enum class G4AnalysisOutput {
  kCsv,
  kRoot,
  kXml,
  kNone
};  

namespace G4Analysis 
{

// Constant expressions
//	
// constexpr G4int kX = 0;  // not yet supported on vc12
// constexpr G4int kY = 1;  // not yet supported on vc12
// constexpr G4int kZ = 2;  // not yet supported on vc12
const G4int kX = 0;  // not yet supported on vc12
const G4int kY = 1;  // not yet supported on vc12
const G4int kZ = 2;  // not yet supported on vc12

// Invalid object Id
// 
// constexpr G4int kInvalidId = -1;  // not yet supported on vc12
const G4int kInvalidId = -1;  // not yet supported on vc12

// Utility functions for checking input parameters
G4bool CheckNbins(G4int nbins);
G4bool CheckMinMax(G4double xmin, G4double xmax,
                   const G4String& fcnName = "none",
                   const G4String& binSchemeName = "linear");
G4bool CheckEdges(const std::vector<G4double>& edges);
G4bool CheckName(const G4String& name, const G4String& objectType);

// Get unit value with added handling of "none" 
G4double GetUnitValue(const G4String& unit);

// Add unit & fcn to the title
void UpdateTitle(G4String& title, 
                 const G4String& unitName, const G4String& fcnName);
                 
// Tokenizer with taking into account composed strings within ""               
void Tokenize(const G4String& line, std::vector<G4String>& tokens);

// make-unique utility, not yet availabla in C++11
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// get output type from name
G4AnalysisOutput GetOutput(const G4String& outputName);
G4String GetOutputName(G4AnalysisOutput outputType);

/*
// make possible to print enumerators in class enum as integer
template <typename Enumeration>
auto as_integer(Enumeration const value)
    -> typename std::underlying_type<Enumeration>::type
{
  return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}
*/

}

#endif




