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

// Class for Root main ntuple management.
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#ifndef G4RootMainNtupleManager_h
#define G4RootMainNtupleManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "G4TNtupleDescription.hh"
#include "G4RootNtupleManager.hh"
#include "globals.hh"

#include <vector>

class G4RootNtupleManager;

namespace tools {
namespace wroot {
class file;
class ntuple;
}
}

class G4RootMainNtupleManager : public G4BaseAnalysisManager
{
  friend class G4RootPNtupleManager;
  friend class G4RootNtupleManager;

  public:
    explicit G4RootMainNtupleManager(G4RootNtupleManager* ntupleBuilder,
                                     G4bool rowWise,
                                     const G4AnalysisManagerState& state);
    ~G4RootMainNtupleManager();

  protected:
    // Types alias
    using NtupleType = tools::wroot::ntuple;
    using NtupleDescriptionType = G4TNtupleDescription<NtupleType>;

    // Methods to manipulate ntuples  
    void   CreateNtuple(const tools::ntuple_booking& ntupleBooking, G4bool warn = true);
    void   CreateNtuplesFromBooking();
    G4bool Merge();
    G4bool Reset(G4bool deleteNtuple);

    // Set methods
    void SetNtupleFile(std::shared_ptr<tools::wroot::file> rfile);
    void SetNtupleDirectory(tools::wroot::directory*  directory);
    std::shared_ptr<tools::wroot::file> GetNtupleFile() const;
    tools::wroot::directory*  GetNtupleDirectory() const;

    // Access functions
    const std::vector<NtupleDescriptionType*>& GetNtupleDescriptionVector() 
      { return fNtupleBuilder->GetNtupleDescriptionVector(); }
    const std::vector<tools::wroot::ntuple*>& GetNtupleVector() 
      { return fNtupleVector; }
    unsigned int GetBasketSize() const;

  private:
    // Data members
    G4RootNtupleManager*  fNtupleBuilder;
    G4bool  fRowWise;
    std::shared_ptr<tools::wroot::file>  fNtupleFile;
    tools::wroot::directory*  fNtupleDirectory;
    std::vector<tools::wroot::ntuple*>   fNtupleVector;
};

inline void G4RootMainNtupleManager::SetNtupleFile(std::shared_ptr<tools::wroot::file> rfile) 
{ fNtupleFile = rfile; }

inline void G4RootMainNtupleManager::SetNtupleDirectory(tools::wroot::directory*  directory)
{ fNtupleDirectory = directory; }

inline std::shared_ptr<tools::wroot::file> G4RootMainNtupleManager::GetNtupleFile() const 
{ return fNtupleFile; }

inline tools::wroot::directory*  G4RootMainNtupleManager::GetNtupleDirectory() const
{ return fNtupleDirectory; }

#endif

