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

// Common implementation for histograms messengers.
// It implements commands in /analysis/hn|pn directories
//
// Author: Ivana Hrivnacova, 26/08/2022  (ivana@ipno.in2p3.fr)

#ifndef G4THnMessenger_h
#define G4THnMessenger_h 1

#include "G4UImessenger.hh"
// #include "G4THnToolsManager.hh"
    // to avoid include recursion this messenger header is included in G4THnToolsManager.hh
    // after the G4THnToolsManager class definition
#include "globals.hh"

#include <memory>
#include <array>
#include <vector>
#include <string_view>

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;
struct G4HnDimension;
struct G4HnDimensionInformation;

template <unsigned int DIM, typename HT>
class G4THnMessenger : public G4UImessenger
{
  public:
    G4THnMessenger(G4THnToolsManager<DIM,HT>* manager);
    G4THnMessenger() = delete;
    ~G4THnMessenger() override = default;

    // Methods
    G4String GetCurrentValue(G4UIcommand* command) final;
    void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    // Helper functions
    G4String GetObjectType() const;
    G4bool IsProfileLastDimension(unsigned int idim) const;
    std::unique_ptr<G4UIcommand> CreateCommand(G4String name, G4String guideline);
    void CreateDimensionParameters(unsigned int idim, 
           std::vector<G4UIparameter*>& parameters) const;
    void AddIdParameter(G4UIcommand& command);
    G4String GetTAddress(G4int id) const;
    G4String GetTVectorAddress() const;

    // Functions to create commands
    void CreateDirectory() const;
    void CreateCmd();
    void SetCmd();
    std::unique_ptr<G4UIcommand> CreateSetBinsCommand(unsigned int ibin);
    void CreateSetTitleCommand();
    std::unique_ptr<G4UIcommand> CreateSetAxisCommand(unsigned int ibin);
    std::unique_ptr<G4UIcommand> CreateSetAxisLogCommand(unsigned int ibin);
    void CreateListCommand();
    void CreateGetCommand();
    void CreateGetVectorCommand();

    // Functions to retrieve data
    void GetBinData(unsigned int idim, G4int& counter, 
           const std::vector<G4String>& parameters, G4HnDimension& bins) const;
    void GetBinInfoData(unsigned int idim, G4int& counter, 
           const std::vector<G4String>& parameters, G4HnDimension& bins, 
           G4HnDimensionInformation& info) const;
    void GetData(G4int& counter, const std::vector<G4String>& parameters,
           std::array<G4HnDimension, DIM>& bins,
           std::array<G4HnDimensionInformation, DIM>& info) const;

    // constants
    static constexpr unsigned int kMaxDim{3};
    static constexpr std::string_view fkClass { "G4THnMessenger" };

    // Data members
    G4THnToolsManager<DIM,HT>* fManager { nullptr }; ///< TO DO Associated class  

    std::unique_ptr<G4UIcommand>  fCreateCmd;
    std::unique_ptr<G4UIcommand>  fSetCmd;
    std::array<std::unique_ptr<G4UIcommand>, DIM>  fSetDimensionCmd;
    std::unique_ptr<G4UIcommand>  fSetTitleCmd;
    std::array<std::unique_ptr<G4UIcommand>, DIM+1> fSetAxisCmd;
    std::array<std::unique_ptr<G4UIcommand>, DIM+1> fSetAxisLogCmd;
    std::unique_ptr<G4UIcommand>  fListCmd;
    std::unique_ptr<G4UIcommand>  fGetTCmd;
    std::unique_ptr<G4UIcommand>  fGetTVectorCmd;

    std::array<unsigned int, DIM> fTmpId;
    std::array<G4HnDimension, DIM> fTmpBins;
    std::array<G4HnDimensionInformation, DIM> fTmpInfo;

    G4String fTValue;
    G4String fTVectorValue;
};

// #include "G4THnMessenger.icc"
    // to avoid include recursion the implementation is included in G4THnToolsManager.hh

#endif
