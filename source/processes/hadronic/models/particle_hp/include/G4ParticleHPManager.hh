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

#ifndef G4ParticleHPManager_h
#define G4ParticleHPManager_h 1

// 121031 First implementation done by T. Koi (SLAC/PPA)
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#include "G4ParticleHPReactionWhiteBoard.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4ParticleDefinition;
class G4ParticleHPChannel;
class G4ParticleHPChannelList;
class G4ParticleHPMessenger;
class G4ParticleHPVector;
class G4PhysicsTable;
struct E_isoAng;
struct E_P_E_isoAng;

class G4ParticleHPManager
{
  public:
    static G4ParticleHPManager* GetInstance();
    ~G4ParticleHPManager();

    G4ParticleHPReactionWhiteBoard* GetReactionWhiteBoard();
    void OpenReactionWhiteBoard();
    void CloseReactionWhiteBoard();

    void GetDataStream(const G4String&, std::istringstream& iss);
    void GetDataStream2(const G4String&, std::istringstream& iss);
    void SetVerboseLevel(G4int i);
    G4int GetVerboseLevel() const { return verboseLevel; };

    void DumpDataSource();

    G4bool GetUseOnlyPhotoEvaporation() const { return USE_ONLY_PHOTONEVAPORATION; };
    G4bool GetSkipMissingIsotopes() const { return SKIP_MISSING_ISOTOPES; };
    G4bool GetNeglectDoppler() const { return NEGLECT_DOPPLER; };
    G4bool GetDoNotAdjustFinalState() const { return DO_NOT_ADJUST_FINAL_STATE; };
    G4bool GetProduceFissionFragments() const { return PRODUCE_FISSION_FRAGMENTS; };
    G4bool GetUseWendtFissionModel() const { return USE_WENDT_FISSION_MODEL; };
    G4bool GetUseNRESP71Model() const { return USE_NRESP71_MODEL; };
    G4bool GetUseDBRC() const { return USE_DBRC; };
    G4bool GetCheckHPNames() const { return CHECK_HP_NAMES; };
    G4bool GetPHPCheck() const { return PHP_CHECK; };
    G4bool GetPHCUsePoisson() const { return PHP_USE_POISSON; };
    G4bool GetDEBUG() const { return DEBUG; };

    const G4String& GetNeutronHPPath() const { return fDataPath[0]; }; 
    const G4String& GetParticleHPPath(const G4ParticleDefinition*) const;
    G4int GetPHPIndex(const G4ParticleDefinition*) const;

    void SetUseOnlyPhotoEvaporation(G4bool val) { USE_ONLY_PHOTONEVAPORATION = val; };
    void SetSkipMissingIsotopes(G4bool val) { SKIP_MISSING_ISOTOPES = val; };
    void SetNeglectDoppler(G4bool val) { NEGLECT_DOPPLER = val; };
    void SetDoNotAdjustFinalState(G4bool val) { DO_NOT_ADJUST_FINAL_STATE = val; };
    void SetProduceFissionFragments(G4bool val)
    {
      // Make sure both fission fragment models are not active at same time
      PRODUCE_FISSION_FRAGMENTS = USE_WENDT_FISSION_MODEL ? false : val;
    };
    void SetUseWendtFissionModel(G4bool val)
    {
      USE_WENDT_FISSION_MODEL = val;
      // Make sure both fission fragment models are not active at same time
      if (USE_WENDT_FISSION_MODEL) PRODUCE_FISSION_FRAGMENTS = false;
    };
    void SetUseNRESP71Model(G4bool val) { USE_NRESP71_MODEL = val; };
    void SetUseDBRC(G4bool val) { USE_DBRC = val; };

    void DumpSetting();

    void RegisterElasticCrossSections(G4PhysicsTable* val) { theElasticCrossSections = val; };
    G4PhysicsTable* GetElasticCrossSections() const { return theElasticCrossSections; };
    void RegisterCaptureCrossSections(G4PhysicsTable* val) { theCaptureCrossSections = val; };
    G4PhysicsTable* GetCaptureCrossSections() const { return theCaptureCrossSections; };
    void RegisterInelasticCrossSections(const G4ParticleDefinition* part, G4PhysicsTable* ptr)
    {
      theInelasticCrossSections[GetPHPIndex(part)] = ptr;
    };
    G4PhysicsTable* GetInelasticCrossSections(const G4ParticleDefinition* part) const
    {
      return theInelasticCrossSections[GetPHPIndex(part)];
    };
    void RegisterFissionCrossSections(G4PhysicsTable* val) { theFissionCrossSections = val; };
    G4PhysicsTable* GetFissionCrossSections() const { return theFissionCrossSections; };

    std::vector<G4ParticleHPChannel*>* GetElasticFinalStates() const { return theElasticFSs; };
    void RegisterElasticFinalStates(std::vector<G4ParticleHPChannel*>* val)
    {
      theElasticFSs = val;
    };

    std::vector<G4ParticleHPChannelList*>*
    GetInelasticFinalStates(const G4ParticleDefinition* part) const
    {
      return theInelasticFSs[GetPHPIndex(part)];
    };
    void RegisterInelasticFinalStates(const G4ParticleDefinition* part,
                                      std::vector<G4ParticleHPChannelList*>* ptr)
    {
      theInelasticFSs[GetPHPIndex(part)] = ptr;
    };

    std::vector<G4ParticleHPChannel*>* GetCaptureFinalStates() const { return theCaptureFSs; };
    void RegisterCaptureFinalStates(std::vector<G4ParticleHPChannel*>* val)
    {
      theCaptureFSs = val;
    };
    std::vector<G4ParticleHPChannel*>* GetFissionFinalStates() const { return theFissionFSs; };
    void RegisterFissionFinalStates(std::vector<G4ParticleHPChannel*>* val)
    {
      theFissionFSs = val;
    };

    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>*
    GetThermalScatteringCoherentCrossSections() const
    {
      return theTSCoherentCrossSections;
    };
    void RegisterThermalScatteringCoherentCrossSections(
      std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* val)
    {
      theTSCoherentCrossSections = val;
    };
    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>*
    GetThermalScatteringIncoherentCrossSections() const
    {
      return theTSIncoherentCrossSections;
    };
    void RegisterThermalScatteringIncoherentCrossSections(
      std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* val)
    {
      theTSIncoherentCrossSections = val;
    };
    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>*
    GetThermalScatteringInelasticCrossSections() const
    {
      return theTSInelasticCrossSections;
    };
    void RegisterThermalScatteringInelasticCrossSections(
      std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* val)
    {
      theTSInelasticCrossSections = val;
    };

    std::map<G4int, std::map<G4double, std::vector<std::pair<G4double, G4double>*>*>*>*
    GetThermalScatteringCoherentFinalStates() const
    {
      return theTSCoherentFinalStates;
    };
    void RegisterThermalScatteringCoherentFinalStates(
      std::map<G4int, std::map<G4double, std::vector<std::pair<G4double, G4double>*>*>*>* val)
    {
      theTSCoherentFinalStates = val;
    };
    std::map<G4int, std::map<G4double, std::vector<E_isoAng*>*>*>*
    GetThermalScatteringIncoherentFinalStates() const
    {
      return theTSIncoherentFinalStates;
    };
    void RegisterThermalScatteringIncoherentFinalStates(
      std::map<G4int, std::map<G4double, std::vector<E_isoAng*>*>*>* val)
    {
      theTSIncoherentFinalStates = val;
    };
    std::map<G4int, std::map<G4double, std::vector<E_P_E_isoAng*>*>*>*
    GetThermalScatteringInelasticFinalStates() const
    {
      return theTSInelasticFinalStates;
    };
    void RegisterThermalScatteringInelasticFinalStates(
      std::map<G4int, std::map<G4double, std::vector<E_P_E_isoAng*>*>*>* val)
    {
      theTSInelasticFinalStates = val;
    };

    G4double GetMinADBRC() const { return theMinADBRC; };
    G4double GetMinEnergyDBRC() const { return theMinEnergyDBRC; };
    G4double GetMaxEnergyDBRC() const { return theMaxEnergyDBRC; };
    G4double GetMaxEnergyDoppler() const { return theMaxEnergyDoppler; };

    void SetMinADBRC(G4double val) { theMinADBRC = val; };
    void SetMinEnergyDBRC(G4double val) { theMinEnergyDBRC = val; };
    void SetMaxEnergyDBRC(G4double val) { theMaxEnergyDBRC = val; };
    void SetMaxEnergyDoppler(G4double val) { theMaxEnergyDoppler = val; };

    G4ParticleHPManager(G4ParticleHPManager &) = delete;
    G4ParticleHPManager & operator=(const G4ParticleHPManager &right) = delete;

  private:
    G4ParticleHPManager();
    void register_data_file(const G4String&, const G4String&);

    static G4ParticleHPManager* instance;

    std::map<G4String, G4String> mDataEvaluation;

    G4int verboseLevel{1};

    G4ParticleHPMessenger* messenger;
    G4bool USE_ONLY_PHOTONEVAPORATION{false};
    G4bool SKIP_MISSING_ISOTOPES{false};
    G4bool NEGLECT_DOPPLER{false};
    G4bool DO_NOT_ADJUST_FINAL_STATE{false};
    G4bool PRODUCE_FISSION_FRAGMENTS{false};
    G4bool USE_WENDT_FISSION_MODEL{false};
    G4bool USE_NRESP71_MODEL{false};
    G4bool USE_DBRC{false};
    G4bool CHECK_HP_NAMES{false};
    G4bool PHP_CHECK{true};
    G4bool PHP_USE_POISSON{false};
    G4bool DEBUG{false};
    G4bool isPrinted{false};

    G4PhysicsTable* theElasticCrossSections{nullptr};
    G4PhysicsTable* theCaptureCrossSections{nullptr};
    G4PhysicsTable* theInelasticCrossSections[6]{nullptr};
    G4PhysicsTable* theFissionCrossSections{nullptr};

    std::vector<G4ParticleHPChannel*>* theElasticFSs{nullptr};
    std::vector<G4ParticleHPChannelList*>* theInelasticFSs[6]{nullptr};
    std::vector<G4ParticleHPChannel*>* theCaptureFSs{nullptr};
    std::vector<G4ParticleHPChannel*>* theFissionFSs{nullptr};

    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* theTSCoherentCrossSections{nullptr};
    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* theTSIncoherentCrossSections{
      nullptr};
    std::map<G4int, std::map<G4double, G4ParticleHPVector*>*>* theTSInelasticCrossSections{nullptr};

    std::map<G4int, std::map<G4double, std::vector<std::pair<G4double, G4double>*>*>*>*
      theTSCoherentFinalStates{nullptr};
    std::map<G4int, std::map<G4double, std::vector<E_isoAng*>*>*>* theTSIncoherentFinalStates{
      nullptr};
    std::map<G4int, std::map<G4double, std::vector<E_P_E_isoAng*>*>*>* theTSInelasticFinalStates{
      nullptr};

    G4double theMinADBRC{200.};
    G4double theMinEnergyDBRC;
    G4double theMaxEnergyDBRC;
    G4double theMaxEnergyDoppler;

    G4String fDataPath[6]{""};
};
#endif
