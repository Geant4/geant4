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

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisUtilities_h
#define G4AnalysisUtilities_h 1

#include "G4Exception.hh"
#include "globals.hh"

#include <vector>
#include <memory>
#include <string_view>

// Enumeration for definition of available output types

enum class G4AnalysisOutput {
  kCsv,
  kHdf5,
  kRoot,
  kXml,
  kNone
};

namespace G4Analysis
{

// Constant expressions
//
constexpr G4int kX { 0 };
constexpr G4int kY { 1 };
constexpr G4int kZ { 2 };
constexpr G4int kInvalidId { -1 };
constexpr G4int kVL0 { 0 };
constexpr G4int kVL1 { 1 };
constexpr G4int kVL2 { 2 };
constexpr G4int kVL3 { 3 };
constexpr G4int kVL4 { 4 };
constexpr unsigned int kDefaultBasketSize = 32000;
constexpr unsigned int kDefaultBasketEntries = 4000;
constexpr std::string_view kNamespaceName { "G4Analysis" };

// Warning
//
void Warn(const G4String& message,
          const std::string_view inClass,
          const std::string_view inFunction);

// Get unit value with added handling of "none"
G4double GetUnitValue(const G4String& unit);

// Tokenizer with taking into account composed strings within ""
void Tokenize(const G4String& line, std::vector<G4String>& tokens);

// Get output type from name
G4AnalysisOutput GetOutput(const G4String& outputName, G4bool warn = true);
size_t GetOutputId(const G4String& outputName, G4bool warn = true);
G4String GetOutputName(G4AnalysisOutput outputType);

// Get short hnType from the tools object
template <typename HT>
G4String GetHnType()
{
  // tools::histo::h1d etc.
  G4String hnTypeLong = HT::s_class();

  // tools::histo::h1d -> h1 etc.
  return hnTypeLong.substr(14, 2);
}

template <typename HT>
G4bool IsProfile()
{
  // tools::histo::h1d etc.
  G4String hnTypeLong = HT::s_class();

  // tools::histo::h1d -> h1 etc.
  return hnTypeLong[14] == 'p';
}

// String conversion
template <typename T>
inline
std::string ToString(const T& value)
{ return std::to_string(value); }

template <>
inline
std::string ToString<std::string>(const std::string& value)
{ return value; }

// File names utilities

// Get file base name (without dot)
G4String GetBaseName(const G4String& fileName);

// Get file base extension (without dot)
G4String GetExtension(const G4String& fileName,
            const G4String& defaultExtension = "");

// Compose and return the histogram or profile specific file name:
// - add _hn_hnName suffix to the file base name
// - add file extension if not present
G4String GetHnFileName(
            const G4String& fileName,
            const G4String& fileType,
            const G4String& hnType,
            const G4String& hnName);

// Update Hn file name:
// - add _vN  suffix to the base namer if cycle > 0
G4String GetHnFileName(
            const G4String& fileName,
            const G4String& fileType,
            G4int cycle = 0);

// Compose and return the ntuple specific file name:
// - add _nt_ntupleName suffix to the file base name
// - add _vN  suffix if cycle > 0
// - add _tN suffix if called on thread worker
// - add file extension if not present
G4String GetNtupleFileName(
            const G4String& fileName,
            const G4String& fileType,
            const G4String& ntupleName,
            G4int cycle = 0);

// Compose and return the ntuple specific file name:
// - add _mFN suffix to the file base name where FN = ntupleFileNumber
// - add _vN  suffix if cycle > 0
// - add file extension if not present
G4String GetNtupleFileName(
            const G4String& fileName,
            const G4String& fileType,
            G4int ntupleFileNumber,
            G4int cycle = 0);

// Update file base name with the thread suffix:
// - add _vN  suffix if cycle > 0
// - add _tN suffix if called on thread worker
// - add file extension if not present
G4String GetTnFileName(
            const G4String& fileName,
            const G4String& fileType,
            G4int cycle = 0);

// Generate plot file name for an output file name
G4String GetPlotFileName(const G4String& fileName);

}

/*
// make possible to print enumerators in class enum as integer
template <typename Enumeration>
auto as_integer(Enumeration const value)
    -> typename std::underlying_type<Enumeration>::type
{
  return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}
*/

#endif

