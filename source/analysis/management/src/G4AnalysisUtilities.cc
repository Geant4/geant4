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
// $Id:$

// Author: Ivana Hrivnacova, 22/08/2013  (ivana@ipno.in2p3.fr)

#include "G4AnalysisUtilities.hh"
#include "G4BinScheme.hh"
#include "G4UnitsTable.hh"
#include "G4String.hh"

namespace {

//_____________________________________________________________________________
G4bool GetToken(const G4String& line, G4String& token, 
                std::string::size_type begIdx, std::string::size_type& endIdx)
{
  while ( line[begIdx] == ' ') ++begIdx; // Loop checking, 23.06.2015, I. Hrivnacova
  if ( line[begIdx] == '"' ) {
    endIdx = line.find('"', begIdx+1);
    if ( endIdx == std::string::npos ) endIdx = line.length();
    token = line.substr(begIdx+1, (endIdx-1)-begIdx);
    ++endIdx;
  }
  else {
    endIdx = line.find(' ', begIdx);
    if ( endIdx == std::string::npos ) endIdx = line.length();
    token = line.substr(begIdx, endIdx-begIdx);
  }
  return ( token.length() > 0 );
}              

}

namespace G4Analysis
{

//_____________________________________________________________________________
G4bool CheckNbins(G4int nbins)
{
  if ( nbins <= 0 ) {
    G4ExceptionDescription description;
    description 
      << "    Illegal value of number of bins: nbins <= 0" << G4endl;
      G4Exception("G4VAnalysisManager::CheckNbins",
                  "Analysis_W013", JustWarning, description);
    return false;
  }
  else
    return true;                   
}  


//_____________________________________________________________________________
G4bool CheckMinMax(G4double xmin, G4double xmax, 
                   const G4String& fcnName, const G4String& binSchemeName)
{
  auto result = true;
  
  if ( xmax <= xmin ) {
    G4ExceptionDescription description;
    description 
      << "    Illegal values of (xmin >= xmax)" << G4endl;
      G4Exception("G4VAnalysisManager::CheckMinMax",
                  "Analysis_W013", JustWarning, description);
                  
    result = false;
  }
  
  if ( ( fcnName != "none" ) && ( binSchemeName != "linear" ) ) {
    G4ExceptionDescription description;
    description 
      << "    Combining Function and Binning scheme is not supported." 
      << G4endl;
      G4Exception("G4VAnalysisManager::CheckMinMax",
                  "Analysis_W013", JustWarning, description);
                  
    result = false;
  }
  
  if ( ( GetBinScheme(binSchemeName) == G4BinScheme::kLog ||
         fcnName == "log" || fcnName == "log10" ) && ( xmin == 0 ) ) {
    G4ExceptionDescription description;
    description 
      << "    Illegal value of (xmin = 0) with logarithmic function or binning" 
      << G4endl;
      G4Exception("G4VAnalysisManager::CheckMinMax",
                  "Analysis_W013", JustWarning, description);
                  
    result = false;
  }
  
  return result;
}  

//_____________________________________________________________________________
G4bool CheckEdges(const std::vector<G4double>& edges)
{
  if ( edges.size() <= 1 ) {
    G4ExceptionDescription description;
    description 
      << "    Illegal edges vector (size <= 1)" << G4endl;
      G4Exception("G4VAnalysisManager::CheckEdges",
                  "Analysis_W013", JustWarning, description);
    return false;
  }
  else
    return true;                   

}

//_____________________________________________________________________________
G4bool CheckName(const G4String& name, const G4String& objectType)
{
  if ( ! name.size() ) {
    G4ExceptionDescription description;
    description 
      << "    Empty " << objectType << " name is not allowed." << G4endl
      << "    " << objectType << " was not created." << G4endl;
      G4Exception("G4VAnalysisManager::CheckName",
                  "Analysis_W013", JustWarning, description);
    return false;
  }
  else
    return true;                   
}

//_____________________________________________________________________________
G4double GetUnitValue(const G4String& unit)
{
   G4double value = 1.;
   if ( unit != "none" ) {
     value = G4UnitDefinition::GetValueOf(unit);
     if ( value == 0. ) value = 1.; 
   }  
   return value;
}   

//_____________________________________________________________________________
void UpdateTitle(G4String& title, 
                 const G4String& unitName, 
                 const G4String& fcnName)
{
  if ( fcnName != "none" )  { title += " "; title += fcnName; title += "("; }
  if ( unitName != "none" ) { title += " ["; title += unitName; title += "]";}
  if ( fcnName != "none" )  { title += ")"; }
}                                                            

//_____________________________________________________________________________
void  Tokenize(const G4String& line, std::vector<G4String>& tokens)
{
  // Define start values
  std::string::size_type begIdx = 0;
  std::string::size_type endIdx = 0; 
  G4String token;
  
  do {
    if ( GetToken(line, token, begIdx, endIdx) ) {
      //G4cout << "got token: '" << token << "'" << G4endl;
      //G4cout << "beg, end: " << begIdx << ", " << endIdx << G4endl;
      tokens.push_back(token);
    }  
    begIdx = endIdx + 1;
  }
  while ( endIdx < line.length() ); // Loop checking, 23.06.2015, I. Hrivnacova
}

//_____________________________________________________________________________
G4AnalysisOutput GetOutput(const G4String& outputName) {
  if      ( outputName == "csv"  )  { return G4AnalysisOutput::kCsv;  }
  else if ( outputName == "root" )  { return G4AnalysisOutput::kRoot; }
  else if ( outputName == "xml"  )  { return G4AnalysisOutput::kXml;  }
  else if ( outputName == "none" )  { return G4AnalysisOutput::kNone; }
  else {
    G4ExceptionDescription description;
    description 
      << "    \"" << outputName << "\" output type is not supported." << G4endl
      << "    " << "Root type will be used.";
    G4Exception("G4Analysis::GetOutputType",
                "Analysis_W013", JustWarning, description);
    return G4AnalysisOutput::kNone; 
  }
}

//_____________________________________________________________________________
G4String GetOutputName(G4AnalysisOutput output) {
  switch ( output ) {
    case G4AnalysisOutput::kCsv:
      return "csv";
      break;
    case G4AnalysisOutput::kRoot: 
      return "root";
      break;
    case G4AnalysisOutput::kXml:
      return "xml";
      break;
    case G4AnalysisOutput::kNone:
      return "none";
      break;
  }
  // should never reach this line
  G4ExceptionDescription description;
  description
    << "    \"" << static_cast<int>(output) << "\" is not handled." << G4endl
    << "    " << "none type will be used.";
  G4Exception("G4Analysis::GetOutputName",
              "Analysis_W013", JustWarning, description);
  return "none";
}

}
