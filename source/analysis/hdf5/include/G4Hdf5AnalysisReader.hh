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

// The main manager for Hdf5 analysis reader.
// It delegates most of functions to the object specific managers.

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5AnalysisReader_h
#define G4Hdf5AnalysisReader_h 1

#include "G4ToolsAnalysisReader.hh"
#include "globals.hh"

#include "toolx/hdf5/ntuple"

#include <memory>
#include <string_view>

class G4Hdf5AnalysisReader;
class G4Hdf5RFileManager;
class G4Hdf5RNtupleManager;
template <class T>
class G4ThreadLocalSingleton;

class G4Hdf5AnalysisReader : public G4ToolsAnalysisReader
{
  friend class G4ThreadLocalSingleton<G4Hdf5AnalysisReader>;

  public:
    ~G4Hdf5AnalysisReader() override;

    // Static methods
    static G4Hdf5AnalysisReader* Instance();

    // Access methods
    toolx::hdf5::ntuple* GetNtuple() const;
    toolx::hdf5::ntuple* GetNtuple(G4int ntupleId) const;
    using G4VAnalysisReader::GetNtuple;

  protected:
    // Virtual methods from base class
    G4bool CloseFilesImpl(G4bool reset) final;

  private:
    G4Hdf5AnalysisReader();

    // Methods
    G4bool Reset();

    // Static data members
    inline static G4Hdf5AnalysisReader* fgMasterInstance { nullptr };
    static constexpr std::string_view fkClass { "G4Hdf5AnalysisReader" };

    // Data members
    std::shared_ptr<G4Hdf5RNtupleManager> fNtupleManager { nullptr };
    std::shared_ptr<G4Hdf5RFileManager>   fFileManager { nullptr };
};

#include "G4Hdf5AnalysisReader.icc"

#endif

