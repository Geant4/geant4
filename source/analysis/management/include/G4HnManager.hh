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

// Class for management of G4HnInformation.
// It implements functions handling the added H1/H2 information
// (not available in g4tools).
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4HnManager_h
#define G4HnManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "G4HnInformation.hh"
#include "globals.hh"

#include <string_view>
#include <utility>
#include <vector>

class G4VFileManager;
class G4HnMessenger;

class G4HnManager : public G4BaseAnalysisManager
{
  public:
    G4HnManager(G4String hnType, const G4AnalysisManagerState& state);
    G4HnManager() = delete;
    ~G4HnManager() override;

    // Methods to manipulate additional information

    G4HnInformation* AddHnInformation(const G4String& name, G4int nofDimensions);

    void ClearData();

    // Access methofd
    G4HnInformation* GetHnInformation(G4int id,
                          std::string_view functionName,
                          G4bool warn = true) const;

    G4HnDimensionInformation* GetHnDimensionInformation(G4int id,
                          G4int dimension,
                          std::string_view functionName,
                          G4bool warn = true) const;

    const std::vector<G4HnInformation*>& GetHnVector() const;
    G4int GetNofHns() const;
    G4int GetNofActiveHns() const;
    G4String GetHnType() const;

    // Activation option

    // Return false if activation is enabled and there is no object activated,
    // return true otherwise
    G4bool IsActive() const;

    // ASCII option

    // Return false if there is no object selected for ASCII output,
    // return true otherwise
    G4bool IsAscii() const;

    // Plotting option

    // Return false if there is no object selected for plotting,
    // return true otherwise
    G4bool IsPlotting() const;

    // Return false if there is no object with a specific file name
    G4bool IsFileName() const;

    // Function implementing public analysis manager interface
    //
    void  SetActivation(G4bool activation);
    void  SetActivation(G4int id, G4bool activation);
    void  SetAscii(G4int id, G4bool ascii);
    void  SetPlotting(G4int id, G4bool plotting);
    void  SetPlotting(G4bool plotting);
    void  SetFileName(G4int id, const G4String& fileName);
    void  SetFileName(const G4String& fileName);
    G4bool  SetAxisIsLog(unsigned int idim, G4int id, G4bool isLogAxis);


    // Access to Hn additional information
    G4String GetName(G4int id) const;
    G4double GetUnit(unsigned int idim, G4int id) const;
    G4bool   GetAxisIsLog(unsigned int idim, G4int id) const;
    G4bool   GetActivation(G4int id) const;
    G4bool   GetAscii(G4int id) const;
    G4bool   GetPlotting(G4int id) const;
    G4String GetFileName(G4int id) const;

    void SetFileManager(std::shared_ptr<G4VFileManager> fileManager);

  private:
    // Methods
    void  SetActivation(G4HnInformation* info, G4bool activation);
    void  SetPlotting(G4HnInformation* info, G4bool plotting);
    void  SetFileName(G4HnInformation* info, const G4String& fileName);

    // Static data members
    static constexpr std::string_view fkClass { "G4HnManager" };

    // Data members
    G4String  fHnType;
    G4int     fNofActiveObjects { 0 };
    G4int     fNofAsciiObjects { 0 };
    G4int     fNofPlottingObjects { 0 };
    G4int     fNofFileNameObjects { 0 };

    // Additional histograms/ntuple properties not included in tools
    std::vector<G4HnInformation*> fHnVector;
    std::shared_ptr<G4VFileManager> fFileManager { nullptr };

    // Messenger
    std::unique_ptr<G4HnMessenger> fMessenger;
};

inline G4int G4HnManager::GetNofHns() const
{ return G4int(fHnVector.size()); }

inline G4int G4HnManager::GetNofActiveHns() const
{ return fNofActiveObjects; }

inline G4String G4HnManager::GetHnType() const
{ return fHnType; }

inline const std::vector<G4HnInformation*>& G4HnManager::GetHnVector() const
{ return fHnVector; }

inline void G4HnManager::SetFileManager(std::shared_ptr<G4VFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

#endif

