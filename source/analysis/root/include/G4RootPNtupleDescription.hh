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

// Structure containing the information related to Root pntuple.
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#ifndef G4RootPNtupleDescription_h
#define G4RootPNtupleDescription_h 1

#include "G4TNtupleDescription.hh"
#include "G4RootFileDef.hh"
#include "globals.hh"

#include "tools/ntuple_booking"
#include "tools/wroot/mt_ntuple_row_wise"
#include "tools/wroot/mt_ntuple_column_wise"

#include <fstream>

namespace tools {
namespace wroot {
class branch;
class ntuple;
}
}

using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;

class G4RootPNtupleDescription
{
  public:
    G4RootPNtupleDescription(G4NtupleBooking* g4NtupleBooking)
      :  fDescription(g4NtupleBooking) {}

    ~G4RootPNtupleDescription()
      {
        if ( fDescription.GetIsNtupleOwner() ) delete fNtuple;
      }

    // Set methods
    void SetNtuple(tools::wroot::imt_ntuple* intuple);
    void SetBasePNtuple(tools::wroot::base_pntuple* basePNtuple);
    void Reset();

    // Get methods
    RootNtupleDescription& GetDescription();
    tools::wroot::imt_ntuple* GetNtuple() const;
    tools::wroot::base_pntuple* GetBasePNtuple() const;
    std::vector<tools::wroot::branch*>& GetMainBranches();

  private:
    RootNtupleDescription fDescription;
    tools::wroot::imt_ntuple* fNtuple { nullptr };
    tools::wroot::base_pntuple* fBasePNtuple { nullptr };
    std::vector<tools::wroot::branch*> fMainBranches;
};

// inline function

inline void G4RootPNtupleDescription::SetNtuple(
  tools::wroot::imt_ntuple* intuple)
{ fNtuple = intuple; }

inline void G4RootPNtupleDescription::SetBasePNtuple(
  tools::wroot::base_pntuple* basePNtuple)
{ fBasePNtuple = basePNtuple; }

inline void G4RootPNtupleDescription::Reset()
{
  if ( fDescription.GetIsNtupleOwner() ) delete fNtuple;
  fNtuple = nullptr;
}

inline RootNtupleDescription&
G4RootPNtupleDescription::GetDescription()
{ return fDescription; }

inline tools::wroot::imt_ntuple*
G4RootPNtupleDescription::GetNtuple() const
{ return fNtuple; }

inline tools::wroot::base_pntuple*
G4RootPNtupleDescription::GetBasePNtuple() const
{ return fBasePNtuple; }

inline std::vector<tools::wroot::branch*>&
G4RootPNtupleDescription::GetMainBranches()
{ return fMainBranches; }

#endif
