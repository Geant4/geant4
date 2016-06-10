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
// $Id: G4RootFileManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// The manager for Root output file operations.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4RootFileManager_h
#define G4RootFileManager_h 1

#include "G4VFileManager.hh"
#include "globals.hh"

#include <map>
#include <memory>

namespace tools {
namespace wroot { 
class file; 
class directory;
}
}  

class G4RootFileManager : public G4VFileManager
{
  public:
    explicit G4RootFileManager(const G4AnalysisManagerState& state);
    virtual ~G4RootFileManager();

    // Methods to manipulate file
    virtual G4bool OpenFile(const G4String& fileName) final;
    virtual G4bool WriteFile() final;
    virtual G4bool CloseFile() final; 

    G4bool CreateHistoDirectory();
    G4bool CreateNtupleDirectory();
    
    // Get methods
    tools::wroot::directory* GetHistoDirectory() const;
    tools::wroot::directory* GetNtupleDirectory() const;

  private:
    // data members
    std::unique_ptr<tools::wroot::file>  fFile;
    tools::wroot::directory*  fHistoDirectory;
    tools::wroot::directory*  fNtupleDirectory;
};

// inline functions

inline tools::wroot::directory* G4RootFileManager::GetHistoDirectory() const
{ return fHistoDirectory; }

inline tools::wroot::directory* G4RootFileManager::GetNtupleDirectory() const
{ return fNtupleDirectory; }

#endif

