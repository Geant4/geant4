#ifndef __PHYSICS_LIST_DEFINED__
#define __PHYSICS_LIST_DEFINED__

#include "G4VUserPhysicsList.hh"
#include "G4DosimPhysicsMessenger.hh"

enum NeutronPsocesses {LEinelastic,Precompound,Mars,Chiral,PrecompoundNew,PrecompoundAll,Undefined};

class G4HadronInelasticProcess;
class G4NeutronInelasticProcess;
class G4ProtonInelasticProcess;
class Mars01EminCut;
class G4ProcessManager;
class Hall;
class PhysicsList : public G4VUserPhysicsList
{
public:
  PhysicsList(Hall* pHall):pLowEnergyModel(NULL),
		pPrecompoundModel(NULL), 
		pMars5GeVModel(NULL),
		pChiralInvModel(NULL),
		pCutProcess(NULL),
		pPrecompoundNew(NULL),
		pPrecompoundAll(NULL),
		pProtonLowEnergyModel(NULL),
		pProtonPrecompoundModel(NULL),
		pProtonMars5GeVModel(NULL),
		pProtonChiralInvModel(NULL),
		pProtonCutProcess(NULL),
		pProtonPrecompNew(NULL),
		pProtonPrecompAll(NULL),
		pAlphaLowEnergy(NULL),
		pAlphaPrecompound(NULL),
		pAlphaMars(NULL),
		pAlphaCut(NULL),
		pAlphaChiral(NULL),
		pDeuteronLowEnergy(NULL),
		pDeuteronPrecompound(NULL),
		pDeuteronMars(NULL),
		pDeuteronCut(NULL),
		pDeuteronChiral(NULL),
		m_cRegisteredModel(Undefined),
			   m_pGeom(pHall){
    pMessenger = new G4DosimPhysicsMessenger(this);
  }
  virtual ~PhysicsList();
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  void ConstructParticleProcess();
  void ConstructGlobal();
  void ChangeModel(char Which);
  void ResetModel(G4ProcessManager* pProtMan,G4ProcessManager* pNeutMan,
		  G4ProcessManager* pAlphaMan,G4ProcessManager* pDeuteronMan);
private:
  G4NeutronInelasticProcess* pLowEnergyModel;
  G4NeutronInelasticProcess* pPrecompoundModel;
  G4NeutronInelasticProcess* pPrecompoundNew;
  G4NeutronInelasticProcess* pPrecompoundAll;
  G4NeutronInelasticProcess* pMars5GeVModel;
  G4NeutronInelasticProcess* pChiralInvModel;
  Mars01EminCut*             pCutProcess;
  G4ProtonInelasticProcess* pProtonLowEnergyModel;
  G4ProtonInelasticProcess* pProtonPrecompoundModel;
  G4ProtonInelasticProcess* pProtonPrecompNew;
  G4ProtonInelasticProcess* pProtonPrecompAll;
  G4ProtonInelasticProcess* pProtonMars5GeVModel;
  G4ProtonInelasticProcess* pProtonChiralInvModel;
  Mars01EminCut*            pProtonCutProcess;
  G4HadronInelasticProcess* pAlphaLowEnergy;
  G4HadronInelasticProcess* pAlphaPrecompound;
  G4HadronInelasticProcess* pAlphaMars;
  Mars01EminCut*            pAlphaCut;
  G4HadronInelasticProcess* pAlphaChiral;
  G4HadronInelasticProcess* pDeuteronLowEnergy;
  G4HadronInelasticProcess* pDeuteronPrecompound;
  G4HadronInelasticProcess* pDeuteronMars;
  Mars01EminCut*            pDeuteronCut;
  G4HadronInelasticProcess* pDeuteronChiral;
  char                  m_cRegisteredModel;
  G4DosimPhysicsMessenger *pMessenger;
  Hall*                     m_pGeom;
};

#endif
