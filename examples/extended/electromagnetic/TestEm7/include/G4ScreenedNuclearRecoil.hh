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
//
// $Id: G4ScreenedNuclearRecoil.hh,v 1.3 2007/12/07 17:51:10 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
//
// Class Description
// Process for screened electromagnetic nuclear elastic scattering; 
// Physics comes from:
// Marcus H. Mendenhall and Robert A. Weller, 
// "Algorithms  for  the rapid  computation  of  classical  cross  sections  
// for  screened  Coulomb  collisions  "
// Nuclear  Instruments  and  Methods  in  Physics  Research  B58  (1991)  11-17  
// The only input required is a screening function phi(r/a) which is the ratio
// of the actual interatomic potential for two atoms with atomic numbers Z1 and Z2,
// to the unscreened potential Z1*Z2*e^2/r where e^2 is elm_coupling in Geant4 units
// the actual screening tables are computed externally in a python module "screened_scattering.py"
// to allow very specific screening functions to be added if desired, without messing
// with the insides of this code.
//
// First version, April 2004, Marcus H. Mendenhall, Vanderbilt University
//
// Class Description - End


#ifndef G4ScreenedNuclearRecoil_h
#define G4ScreenedNuclearRecoil_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "c2_function.hh"

#include <map>
#include <vector>

class G4VParticleChange;

typedef struct G4ScreeningTables { 
	G4double z1, z2, m1, m2, au, emin;
	c2_function<G4double> *EMphiData; 
} G4ScreeningTables;

// A class for loading ScreenedCoulombCrossSections
class G4ScreenedCoulombCrossSectionInfo
{
public:
        G4ScreenedCoulombCrossSectionInfo() { }
        ~G4ScreenedCoulombCrossSectionInfo() { }

        const char *CVSHeaderVers() { return
                "";
        }

        const char *CVSFileVers() { return ""; }
};

// A class for loading ScreenedCoulombCrossSections
class G4ScreenedCoulombCrossSection : public G4ScreenedCoulombCrossSectionInfo
{
public:

	G4ScreenedCoulombCrossSection() : verbosity(1) { }
	G4ScreenedCoulombCrossSection(const G4ScreenedCoulombCrossSection &src) : 
		G4ScreenedCoulombCrossSectionInfo(),verbosity(src.verbosity) { }
	virtual ~G4ScreenedCoulombCrossSection();
	
	typedef std::map<G4int, G4ScreeningTables> ScreeningMap;
	
	// a local, fast-access mapping of a particle's Z to its full definition
	typedef std::map<G4int, class G4ParticleDefinition *> ParticleCache;
	
	// LoadData is called by G4ScreenedNuclearRecoil::GetMeanFreePath
	// It loads the data tables, builds the elemental cross-section tables.
	virtual void LoadData(G4String screeningKey, G4int z1, G4double m1, G4double recoilCutoff) = 0;
	
	// BuildMFPTables is called by G4ScreenedNuclearRecoil::GetMeanFreePath to build the MFP tables for each material
	void BuildMFPTables(void); // scan the MaterialsTable and construct MFP tables
	
	virtual G4ScreenedCoulombCrossSection *create() = 0; // a 'virtual constructor' which clones the class
	const G4ScreeningTables *GetScreening(G4int Z)  { return &(screeningData[Z]); }
	void SetVerbosity(G4int v) { verbosity=v; }
	
	// this process needs element selection weighted only by number density
	G4ParticleDefinition* SelectRandomUnweightedTarget(const G4MaterialCutsCouple* couple);
	
	enum { nMassMapElements=116 };
		
	G4double standardmass(G4int z1) { return z1 <= nMassMapElements ? massmap[z1] : 2.5*z1; }
	
	// get the mean-free-path table for the indexed material 
	c2_function<G4double> * operator [] (G4int materialIndex) { 
			return MFPTables.find(materialIndex)!=MFPTables.end() ? MFPTables[materialIndex] : (c2_function<G4double> *)0;
		}
	
protected:
	ScreeningMap screeningData; // screening tables for each element
	ParticleCache targetMap;
	G4int verbosity;
	std::map<G4int, c2_function<G4double> *> sigmaMap; // total cross section for each element
	std::map<G4int, c2_function<G4double> *> MFPTables; // MFP for each material
	
private:
	static const G4double massmap[nMassMapElements+1];

};

typedef struct G4CoulombKinematicsInfo {
	G4double impactParameter;
	G4ScreenedCoulombCrossSection *crossSection;
	G4double a1, a2, sinTheta, cosTheta, sinZeta, cosZeta, eRecoil;
	G4ParticleDefinition *recoilIon;	} G4CoulombKinematicsInfo;

class G4ScreenedCollisionStage {
public:
	virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
			const class G4Track& aTrack, const class G4Step& aStep)=0;
	virtual ~G4ScreenedCollisionStage() {}
};

class G4ScreenedCoulombClassicalKinematics: public G4ScreenedCoulombCrossSectionInfo, public G4ScreenedCollisionStage {

public:
	G4ScreenedCoulombClassicalKinematics() { }
	virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
		const class G4Track& aTrack, const class G4Step& aStep);
	
	G4bool DoScreeningComputation(class G4ScreenedNuclearRecoil *master,
		const G4ScreeningTables *screen, 
		G4double eps, G4double beta);
	virtual ~G4ScreenedCoulombClassicalKinematics() {}
};

class G4SingleScatter: public G4ScreenedCoulombCrossSectionInfo, public G4ScreenedCollisionStage {

public:
	G4SingleScatter() { }
	virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
		const class G4Track& aTrack, const class G4Step& aStep);
	virtual ~G4SingleScatter() {}
};

class G4ScreenedNuclearRecoil : public G4ScreenedCoulombCrossSectionInfo, public G4VDiscreteProcess
{
public:
	
	friend class G4ScreenedCollisionStage;

	G4ScreenedNuclearRecoil(const G4String& processName = "ScreenedElastic",
			const G4String &ScreeningKey="zbl", G4bool GenerateRecoils=1, 
			G4double RecoilCutoff=100.0*eV, G4double PhysicsCutoff=10.0*eV);
	
	virtual ~G4ScreenedNuclearRecoil();
	
	virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* );
	
	virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
	
	virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	
	virtual void BuildPhysicsTable(const G4ParticleDefinition&) { }
	
	virtual void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);
	
	virtual G4bool CheckNuclearCollision(G4double A, G4double A1, G4double apsis); // return true if hard collision

	virtual G4ScreenedCoulombCrossSection *GetNewCrossSectionHandler(void);
	
	G4double GetNIEL() const { return NIEL; } // Get non-ionizing energy loss for last step
	
	void ResetTables(); // clear all data tables to allow changing energy cutoff, materials, etc.
	
	std::string GetScreeningKey() const { return screeningKey; }
	void AllowEnergyDeposition(G4bool flag) { registerDepositedEnergy=flag; }
	G4bool GetAllowEnergyDeposition() const { return registerDepositedEnergy; }
	void EnableRecoils(G4bool flag) { generateRecoils=flag; }
	G4bool GetEnableRecoils() const { return generateRecoils; }
	void SetMFPScaling(G4double scale) { MFPScale=scale; }
	G4double GetMFPScaling() const { return MFPScale; }
	void AvoidNuclearReactions(G4bool flag) { avoidReactions=flag; }
	G4bool GetAvoidNuclearReactions() const { return avoidReactions; }
	void SetRecoilCutoff(G4double energy) { recoilCutoff=energy; }
	G4double GetRecoilCutoff() const { return recoilCutoff; }
	void SetPhysicsCutoff(G4double energy) { physicsCutoff=energy; ResetTables(); }
	G4double GetPhysicsCutoff() const { return physicsCutoff; }
	class G4ParticleChange &GetParticleChange() { return aParticleChange; }
	void AddToNIEL(G4double energy) { NIEL+=energy; }
	void SetCrossSectionHardening(G4double fraction, G4double HardeningFactor) {
		hardeningFraction=fraction;
		hardeningFactor=HardeningFactor;
	}
	G4double GetHardeningFraction() const { return hardeningFraction; }
	G4double GetHardeningFactor() const { return hardeningFactor; }
	G4double GetCurrentInteractionLength() const { return currentInteractionLength; }
	void SetExternalCrossSectionHandler(G4ScreenedCoulombCrossSection *cs) { 
		externalCrossSectionConstructor=cs; 
	}
	G4int GetVerboseLevel() const { return verboseLevel; }
	std::map<G4int, G4ScreenedCoulombCrossSection*> &GetCrossSectionHandlers() 
		{ return crossSectionHandlers; }
	void ClearStages(void); 
	void AddStage(G4ScreenedCollisionStage *stage) { collisionStages.push_back(stage); }
	G4CoulombKinematicsInfo &GetKinematics() { return kinematics; }
	void SetValidCollision(G4bool flag) { validCollision=flag; }
	G4bool GetValidCollision() const { return validCollision; }
	
protected:
	G4double highEnergyLimit, lowEnergyLimit;
	G4String screeningKey;
	G4bool generateRecoils, avoidReactions;
	G4double recoilCutoff, physicsCutoff;
	G4bool registerDepositedEnergy;
	G4double NIEL;
	G4double MFPScale;
	G4double hardeningFraction, hardeningFactor;
	
	G4ScreenedCoulombCrossSection *externalCrossSectionConstructor;
	std::vector<G4ScreenedCollisionStage *> collisionStages;
	
	std::map<G4int, G4ScreenedCoulombCrossSection*> crossSectionHandlers;		
	std::map<G4int, c2_function<G4double>*> meanFreePathTables;
	
	G4bool validCollision;
	G4CoulombKinematicsInfo kinematics;

};

// A customized G4CrossSectionHandler which gets its data from an external program
class G4NativeScreenedCoulombCrossSection: public G4ScreenedCoulombCrossSection 
{
public:
	G4NativeScreenedCoulombCrossSection();

	G4NativeScreenedCoulombCrossSection(const G4NativeScreenedCoulombCrossSection &src) 
		: G4ScreenedCoulombCrossSection(src), phiMap(src.phiMap) { }
	
	G4NativeScreenedCoulombCrossSection(const G4ScreenedCoulombCrossSection &src) : G4ScreenedCoulombCrossSection(src)  { }
	virtual ~G4NativeScreenedCoulombCrossSection();
	
	virtual void LoadData(G4String screeningKey, G4int z1, G4double m1, G4double recoilCutoff);
	virtual G4ScreenedCoulombCrossSection *create() 
		{ return new G4NativeScreenedCoulombCrossSection(*this); } 
	// get a list of available keys
	std::vector<G4String> GetScreeningKeys() const;
	
	typedef c2_function<G4double> &(*ScreeningFunc)(G4int z1, G4int z2, size_t nPoints, G4double rMax, G4double *au);
	
	void AddScreeningFunction(G4String name, ScreeningFunc fn) {
		phiMap[name]=fn;
	}
	
private:
		// this is a map used to look up screening function generators
		std::map<std::string, ScreeningFunc> phiMap;
};



#endif
