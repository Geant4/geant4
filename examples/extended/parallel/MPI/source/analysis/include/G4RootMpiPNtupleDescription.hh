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

// Structure containing the information related to Root MPI pntuple.
// This class is temporarily provided with g4mpi,
// it will be integrated in Geant4 analysis category in future.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4RootMpiPNtupleDescription_h
#define G4RootMpiPNtupleDescription_h 1

#include "tools/impi"
#include "tools/ntuple_booking"
#include "tools/wroot/impi_ntuple"
#include "tools/wroot/ntuple"

#include "G4RootFileDef.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

namespace tools
{
namespace wroot
{
class branch;
class base_pntuple;
}  // namespace wroot
}  // namespace tools

using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;

class G4RootMpiPNtupleDescription
{
  public:
    G4RootMpiPNtupleDescription(G4NtupleBooking* g4NtupleBooking) : fDescription(g4NtupleBooking) {}

    ~G4RootMpiPNtupleDescription()
    {
      if (fDescription.GetIsNtupleOwner()) delete fNtuple;
    }

    // Set methods
    void SetNtuple(tools::wroot::impi_ntuple* intuple);
    void SetBasePNtuple(tools::wroot::base_pntuple* basePNtuple);
    void SetMainNtupleRank(G4int value);
    void SetImpi(tools::impi* value);
    void Reset();

    // Get methods
    RootNtupleDescription& GetDescription();
    tools::wroot::impi_ntuple* GetNtuple() const;
    tools::wroot::base_pntuple* GetBasePNtuple() const;
    std::vector<tools::wroot::branch*>& GetMainBranches();
    G4int GetMainNtupleRank() const;

  private:
    RootNtupleDescription fDescription;
    tools::wroot::impi_ntuple* fNtuple{nullptr};
    tools::wroot::base_pntuple* fBasePNtuple{nullptr};
    std::vector<tools::wroot::branch*> fMainBranches;
    G4int fMainNtupleRank{0};
    tools::impi* fImpi{nullptr};
};

// inline function

inline void G4RootMpiPNtupleDescription::SetNtuple(tools::wroot::impi_ntuple* intuple)
{
  fNtuple = intuple;
}

inline void G4RootMpiPNtupleDescription::SetBasePNtuple(tools::wroot::base_pntuple* basePNtuple)
{
  fBasePNtuple = basePNtuple;
}

inline void G4RootMpiPNtupleDescription::SetMainNtupleRank(G4int value)
{
  fMainNtupleRank = value;
}

inline void G4RootMpiPNtupleDescription::SetImpi(tools::impi* value)
{
  fImpi = value;
}

inline void G4RootMpiPNtupleDescription::Reset()
{
  if (fDescription.GetIsNtupleOwner()) delete fNtuple;
  fNtuple = nullptr;
}

inline RootNtupleDescription& G4RootMpiPNtupleDescription::GetDescription()
{
  return fDescription;
}

inline tools::wroot::impi_ntuple* G4RootMpiPNtupleDescription::GetNtuple() const
{
  return fNtuple;
}

inline tools::wroot::base_pntuple* G4RootMpiPNtupleDescription::GetBasePNtuple() const
{
  return fBasePNtuple;
}

inline std::vector<tools::wroot::branch*>& G4RootMpiPNtupleDescription::GetMainBranches()
{
  return fMainBranches;
}

inline G4int G4RootMpiPNtupleDescription::GetMainNtupleRank() const
{
  return fMainNtupleRank;
}

#endif
