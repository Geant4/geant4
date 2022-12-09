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

// Structure templated containing the information related to an ntuple
// for all output types.
//
// Author: Ivana Hrivnacova, 19/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4TNtupleDescription_h
#define G4TNtupleDescription_h 1

#include "G4NtupleBookingManager.hh"
#include "globals.hh"

#include "tools/ntuple_booking"

#include <fstream>

template <typename NT, typename FT>
class G4TNtupleDescription
{
  public:
    G4TNtupleDescription(G4NtupleBooking* g4NtupleBooking)
      :  fNtupleBooking(g4NtupleBooking->fNtupleBooking),
         fFileName(g4NtupleBooking->fFileName),
         fActivation(g4NtupleBooking->fActivation)
      {}

    G4TNtupleDescription() = delete;
    ~G4TNtupleDescription()
      {
        if ( fIsNtupleOwner ) delete fNtuple;
      }

    // Set methods
    void SetFile(std::shared_ptr<FT> file);
    void SetNtuple(NT* ntuple);
    void SetFileName(const G4String& fileName);
    void SetActivation(G4bool activation);
    void SetIsNtupleOwner(G4bool isNtupleOwner);
    void SetHasFill(G4bool hasFill);
    void Reset();

    // Get methods
    std::shared_ptr<FT> GetFile() const;
    NT* GetNtuple() const;
    const tools::ntuple_booking& GetNtupleBooking() const;
    G4String GetFileName() const;
    G4bool GetActivation() const;
    G4bool GetIsNtupleOwner() const;
    G4bool GetHasFill() const;

  private:
    std::shared_ptr<FT> fFile { nullptr };
    NT* fNtuple { nullptr };
    tools::ntuple_booking fNtupleBooking;
    G4String fFileName;
    G4bool fActivation { true };
    G4bool fIsNtupleOwner { true };
    G4bool fHasFill { false };
};

// inline functions

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetFile(std::shared_ptr<FT> file)
{ fFile = file; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetNtuple(NT* ntuple)
{ fNtuple = ntuple; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetFileName(const G4String& fileName)
{ fFileName = fileName; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetActivation(G4bool activation)
{ fActivation = activation; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetIsNtupleOwner(G4bool isNtupleOwner)
{ fIsNtupleOwner = isNtupleOwner; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::SetHasFill(G4bool hasFill)
{ fHasFill = hasFill; }

template <typename NT, typename FT>
void G4TNtupleDescription<NT, FT>::Reset()
{
  if (fIsNtupleOwner) delete fNtuple;
  fNtuple = nullptr;
}

template <typename NT, typename FT>
std::shared_ptr<FT> G4TNtupleDescription<NT, FT>::GetFile() const
{ return fFile; }

template <typename NT, typename FT>
NT* G4TNtupleDescription<NT, FT>::GetNtuple() const
{ return fNtuple; }

template <typename NT, typename FT>
const tools::ntuple_booking&
G4TNtupleDescription<NT, FT>::GetNtupleBooking() const
{ return fNtupleBooking; }

template <typename NT, typename FT>
G4String G4TNtupleDescription<NT, FT>::GetFileName() const
{ return fFileName; }

template <typename NT, typename FT>
G4bool G4TNtupleDescription<NT, FT>::GetActivation() const
{ return fActivation; }

template <typename NT, typename FT>
G4bool G4TNtupleDescription<NT, FT>::GetIsNtupleOwner() const
{ return fIsNtupleOwner; }

template <typename NT, typename FT>
G4bool G4TNtupleDescription<NT, FT>::GetHasFill() const
{ return fHasFill; }

#endif
