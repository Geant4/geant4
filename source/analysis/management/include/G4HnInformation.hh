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

#include "globals.hh"
#include "G4Fcn.hh"
#include "G4BinScheme.hh"
#include "G4AnalysisUtilities.hh"

// The additional Hn information per dimension

struct G4HnDimensionInformation
{
  G4HnDimensionInformation(
      const G4String& unitName,
      const G4String& fcnName,
      G4double unit,
      G4Fcn fcn,
      G4BinScheme binScheme)
    : fUnitName(unitName),
      fFcnName(fcnName),
      fUnit(unit),
      fFcn(fcn),
      fBinScheme(binScheme)
      {}

  G4HnDimensionInformation() = default;
  G4HnDimensionInformation(const G4HnDimensionInformation& rhs) = default;
  G4HnDimensionInformation& operator=(const G4HnDimensionInformation& rhs) = default;

  //G4String fName;
  G4String fUnitName;
  G4String fFcnName;
  G4double fUnit;
  G4Fcn    fFcn;
  G4BinScheme fBinScheme;
};

class G4HnInformation
{
  public:
    G4HnInformation(const G4String& name, G4int nofDimensions)
      : fName(name)
    { fHnDimensionInformations.reserve(nofDimensions); }

    // Deleted default constructor
    G4HnInformation() = delete;

    // Set methods
    void AddHnDimensionInformation(
            const G4HnDimensionInformation& hnDimensionInformation);
    void AddDimension(
            const G4String& unitName, const G4String& fcnName, G4BinScheme binScheme);
    void SetDimension(G4int dimension,
            const G4String& unitName, const G4String& fcnName, G4BinScheme binScheme);
    void SetIsLogAxis(G4int axis, G4bool isLog);
    void SetActivation(G4bool activation);
    void SetAscii(G4bool ascii);
    void SetPlotting(G4bool plotting);
    void SetFileName(G4String fileName);

    // Get methods
    G4String GetName() const;
    G4HnDimensionInformation* GetHnDimensionInformation(G4int dimension);
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

// inline functions

inline void G4HnInformation::AddHnDimensionInformation(
  const G4HnDimensionInformation& hnDimensionInformation)
{ fHnDimensionInformations.push_back(hnDimensionInformation); }

inline void G4HnInformation::AddDimension(
  const G4String& unitName, const G4String& fcnName, G4BinScheme binScheme)
{
  auto unit = G4Analysis::GetUnitValue(unitName);
  auto fcn = G4Analysis::GetFunction(fcnName);
  fHnDimensionInformations.emplace_back(
    unitName, fcnName, unit, fcn, binScheme);
}

inline void G4HnInformation::SetDimension(G4int dimension,
  const G4String& unitName, const G4String& fcnName, G4BinScheme binScheme)
{
  auto info = GetHnDimensionInformation(dimension);
  auto unit = G4Analysis::GetUnitValue(unitName);
  auto fcn = G4Analysis::GetFunction(fcnName);
  info->fUnitName = unitName;
  info->fFcnName = fcnName;
  info->fUnit = unit;
  info->fFcn = fcn;
  info->fBinScheme = binScheme;
}

inline void G4HnInformation::SetIsLogAxis(G4int axis, G4bool isLog)
{ fIsLogAxis[axis] = isLog; }

inline void G4HnInformation::SetActivation(G4bool activation)
{ fActivation = activation; }

inline void G4HnInformation::SetAscii(G4bool ascii)
{ fAscii = ascii; }

inline void G4HnInformation::SetPlotting(G4bool plotting)
{ fPlotting = plotting; }

inline void G4HnInformation::SetFileName(G4String fileName)
{ fFileName = fileName; }

inline G4String G4HnInformation::GetName() const
{ return fName; }

inline G4HnDimensionInformation* G4HnInformation::GetHnDimensionInformation(G4int dimension)
{ return &(fHnDimensionInformations[dimension]); }

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

#endif
