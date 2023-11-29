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

// Data class for the added Hn/Pn information (not available in g4tools).
//
// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4HnInformation_h
#define G4HnInformation_h 1

#include "G4AnalysisUtilities.hh"
#include "G4BinScheme.hh"
#include "G4Fcn.hh"
#include "globals.hh"

#include <utility>
#include <vector>

// The histogram input parameters per dimension
struct G4HnDimension
{
  G4HnDimension(
      G4int nbins,
      G4double minValue,
      G4double maxValue)
    : fNBins(nbins),
      fMinValue(minValue),
      fMaxValue(maxValue)
      {
        fEdges.clear();
      }

  G4HnDimension(const std::vector<G4double>& edges)
    : fNBins(0),
      fMinValue(0.),
      fMaxValue(0.),
      fEdges(edges)
      {}

  G4HnDimension() = default;
  G4HnDimension(const G4HnDimension& rhs) = default;
  G4HnDimension& operator=(const G4HnDimension& rhs) = default;

  void Print() const;

  G4int fNBins{0};
  G4double fMinValue{0.};
  G4double fMaxValue{0.};
  std::vector<G4double> fEdges;
};

// The additional histogram information per dimension
struct G4HnDimensionInformation
{
  G4HnDimensionInformation(
      G4String unitName,
      G4String fcnName,
      G4String binSchemeName = "linear")
    : fUnitName(std::move(unitName)),
      fFcnName(std::move(fcnName)),
      fBinSchemeName(std::move(binSchemeName)),
      fUnit(G4Analysis::GetUnitValue(fUnitName)),
      fFcn(G4Analysis::GetFunction(fFcnName)),
      fBinScheme(G4Analysis::GetBinScheme(fBinSchemeName))
      {}

  G4HnDimensionInformation()
    : G4HnDimensionInformation("none", "none", "linear") {}
  G4HnDimensionInformation(const G4HnDimensionInformation& rhs) = default;
  G4HnDimensionInformation& operator=(const G4HnDimensionInformation& rhs) = default;

  void Print() const;

  G4String fUnitName;
  G4String fFcnName;
  G4String fBinSchemeName;
  G4double fUnit;
  G4Fcn    fFcn;
  G4BinScheme fBinScheme;
};

// The additional histogram information
class G4HnInformation
{
  public:
    G4HnInformation(G4String name, G4int nofDimensions)
      : fName(std::move(name))
    { fHnDimensionInformations.reserve(nofDimensions); }

    // Deleted default constructor
    G4HnInformation() = delete;

    // Set methods
    void AddDimension(const G4HnDimensionInformation& hnDimensionInformation);
    void SetDimension(G4int dimension, const G4HnDimensionInformation& hnDimensionInformation);

    void SetIsLogAxis(G4int axis, G4bool isLog);
    void SetActivation(G4bool activation);
    void SetAscii(G4bool ascii);
    void SetPlotting(G4bool plotting);
    void SetFileName(const G4String& fileName);

    // Get methods
    G4String GetName() const;
    G4HnDimensionInformation* GetHnDimensionInformation(G4int dimension);
    const G4HnDimensionInformation& GetHnDimensionInformation(G4int dimension) const;
    G4bool  GetIsLogAxis(G4int axis) const;
    G4bool  GetActivation() const;
    G4bool  GetAscii() const;
    G4bool  GetPlotting() const;
    G4String GetFileName() const;

  private:
    // Data members
    G4String fName;
    std::vector<G4HnDimensionInformation> fHnDimensionInformations;
    std::vector<G4bool> fIsLogAxis { false, false, false };
    G4bool   fActivation { true };
    G4bool   fAscii { false };
    G4bool   fPlotting { false };
    G4String fFileName;
};

namespace G4Analysis
{

// Apply Hn information
void Update(G4double& value, const G4HnDimensionInformation& hnInfo);
void UpdateValues(G4HnDimension& bins, const G4HnDimensionInformation& hnInfo);
void Update(G4HnDimension& bins, const G4HnDimensionInformation& hnInfo);
void UpdateTitle(G4String& title, const G4HnDimensionInformation& hnInfo);

// Paremeters check
G4bool CheckMinMax(G4double min, G4double max);
G4bool CheckDimension(unsigned int idim,
         const G4HnDimension& dimension, const G4HnDimensionInformation& info);

template <unsigned int DIM>
G4bool CheckDimensions(
         const std::array<G4HnDimension, DIM>& bins,
         const std::array<G4HnDimensionInformation, DIM>& hnInfo,
         G4bool isProfile = false);
}

// inline functions

inline void G4HnInformation::AddDimension(
  const G4HnDimensionInformation& hnDimensionInformation)
{ fHnDimensionInformations.push_back(hnDimensionInformation); }

inline void G4HnInformation::SetDimension(
  G4int dimension, const G4HnDimensionInformation& hnDimensionInformation)
{
  auto info = GetHnDimensionInformation(dimension);
  (*info) = hnDimensionInformation;
}

inline void G4HnInformation::SetIsLogAxis(G4int axis, G4bool isLog)
{ fIsLogAxis[axis] = isLog; }

inline void G4HnInformation::SetActivation(G4bool activation)
{ fActivation = activation; }

inline void G4HnInformation::SetAscii(G4bool ascii)
{ fAscii = ascii; }

inline void G4HnInformation::SetPlotting(G4bool plotting)
{ fPlotting = plotting; }

inline void G4HnInformation::SetFileName(const G4String& fileName)
{ fFileName = fileName; }

inline G4String G4HnInformation::GetName() const
{ return fName; }

inline G4HnDimensionInformation* G4HnInformation::GetHnDimensionInformation(G4int dimension)
{ return &(fHnDimensionInformations[dimension]); }

inline const G4HnDimensionInformation& G4HnInformation::GetHnDimensionInformation(G4int dimension) const
{ return fHnDimensionInformations[dimension]; }

inline G4bool  G4HnInformation::GetIsLogAxis(G4int axis) const
{ return fIsLogAxis[axis]; }

inline G4bool  G4HnInformation::GetActivation() const
{ return fActivation; }

inline G4bool  G4HnInformation::GetAscii() const
{ return fAscii; }

inline G4bool  G4HnInformation::GetPlotting() const
{ return fPlotting; }

inline G4String G4HnInformation::GetFileName() const
{ return fFileName; }

template <unsigned int DIM>
inline G4bool G4Analysis::CheckDimensions(
         const std::array<G4HnDimension, DIM>& bins,
         const std::array<G4HnDimensionInformation, DIM>& hnInfo,
         G4bool isProfile)
{
  G4bool result = true;

  // Check bins parameters
  // (the last dimension has special meaninh for profiles)
  auto dimToCheck = (isProfile) ? DIM -1 : DIM ;
  for (unsigned int idim = 0; idim < dimToCheck; ++idim) {
    result &= CheckDimension(idim, bins[idim], hnInfo[idim]);
  }

  // Check profile min/max value
  if (isProfile) {
    result &= CheckMinMax(bins[DIM-1].fMinValue, bins[DIM-1].fMaxValue);
  }

  return result;
}

#endif
