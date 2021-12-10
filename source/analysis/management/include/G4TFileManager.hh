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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4TFileManager_h
#define G4TFileManager_h 1

#include "G4TFileInformation.hh"

#include <vector>
#include <memory>
#include <string_view>

template <typename FT>
class G4TFileManager
{
  public:
    explicit G4TFileManager(const G4AnalysisManagerState& state);
    G4TFileManager() = delete;
    virtual ~G4TFileManager();

    // Per file methods
    // (implements all tests and warnings, fills the map)
    std::shared_ptr<FT> CreateTFile(const G4String& fileName);
    G4bool WriteTFile(const G4String& fileName);
    G4bool CloseTFile(const G4String& fileName);
    G4bool SetIsEmpty(const G4String& fileName, G4bool isEmpty);
    //
    std::shared_ptr<FT> GetTFile(const G4String& fileName, G4bool warn = true) const;

    // Methods applied to all registered files
    G4bool OpenFiles();
    G4bool WriteFiles();
    G4bool CloseFiles();
    G4bool DeleteEmptyFiles();

    // Clear all data
    void ClearData();

  protected:
    // Methods to be implemented per file type manipulate file
    virtual std::shared_ptr<FT> CreateFileImpl(const G4String& fileName) = 0;
    virtual G4bool WriteFileImpl(std::shared_ptr<FT> file) = 0;
    virtual G4bool CloseFileImpl(std::shared_ptr<FT> file) = 0;

  private:
    // Methods
    void FileNotFoundWarning(const G4String& fileName,
          std::string_view functionName) const;
    G4TFileInformation<FT>* GetFileInfoInFunction(const G4String& fileName,
      std::string_view functionName, G4bool warn = true) const;
    std::shared_ptr<FT> GetFileInFunction(const G4String& fileName,
      std::string_view functionName, G4bool warn = true) const;

    G4bool WriteTFile(std::shared_ptr<FT> file, const G4String& fileName);
    G4bool CloseTFile(std::shared_ptr<FT> file, const G4String& fileName);
    G4bool DeleteEmptyFile(const G4String& fileName);

    // Static data members
    static constexpr std::string_view fkClass { "G4TFileManager<FT>" };

    // Data members
    const G4AnalysisManagerState& fAMState;
    std::map<G4String, G4TFileInformation<FT>*> fFileMap;
 };

#include "G4TFileManager.icc"

#endif
