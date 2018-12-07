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
/// \file electromagnetic/TestEm7/include/G4ScreenedNuclearRecoil.hh
/// \brief Definition of the G4ScreenedNuclearRecoil class
//
//
//
// G4ScreenedNuclearRecoil.hh,v 1.24 2008/05/01 19:58:59 marcus Exp
// GEANT4 tag 
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
// of the actual interatomic potential for two atoms with atomic numbers Z1 and
// Z2,
// to the unscreened potential Z1*Z2*e^2/r where e^2 is elm_coupling in Geant4 
// units
// the actual screening tables are computed externally in a python module 
// "screened_scattering.py"
// to allow very specific screening functions to be added if desired, without 
// messing with the insides of this code.
//
// First version, April 2004, Marcus H. Mendenhall, Vanderbilt University
// May 1, 2008 -- Added code to allow process to have zero cross section above 
//  max energy, to coordinate with G4MSC.  -- mhm
//
// Class Description - End


#ifndef G4ScreenedNuclearRecoil_h
#define G4ScreenedNuclearRecoil_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "c2_function.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <map>
#include <vector>

class G4VNIELPartition;

typedef c2_const_ptr<G4double> G4_c2_const_ptr;
typedef c2_ptr<G4double> G4_c2_ptr;
typedef c2_function<G4double> G4_c2_function;

typedef struct G4ScreeningTables { 
        G4double z1, z2, m1, m2, au, emin;
        G4_c2_const_ptr EMphiData; 
} G4ScreeningTables;

// A class for loading ScreenedCoulombCrossSections
class G4ScreenedCoulombCrossSectionInfo
{
public:
        G4ScreenedCoulombCrossSectionInfo() { }
        ~G4ScreenedCoulombCrossSectionInfo() { }
        
        static const char* CVSHeaderVers() { return 
 "G4ScreenedNuclearRecoil.hh,v 1.24 2008/05/01 19:58:59 marcus Exp GEANT4 tag ";
        }
        static const char* CVSFileVers();
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
      virtual void LoadData(G4String screeningKey, G4int z1, G4double m1,
                            G4double recoilCutoff) = 0;
        
      // BuildMFPTables is called by G4ScreenedNuclearRecoil::GetMeanFreePath 
      //to build the MFP tables for each material
      void BuildMFPTables(void); // scan the MaterialsTable and construct MFP 
                                 //tables
        
      virtual G4ScreenedCoulombCrossSection *create() = 0; 
      // a 'virtual constructor' which clones the class
      const G4ScreeningTables *GetScreening(G4int Z)
                     { return &(screeningData[Z]); }
      void SetVerbosity(G4int v) { verbosity=v; }
        
      // this process needs element selection weighted only by number density
      G4ParticleDefinition* SelectRandomUnweightedTarget
                                    (const G4MaterialCutsCouple* couple);
        
      enum { nMassMapElements=116 };
                
      G4double standardmass(G4int z1) 
            { return z1 <= nMassMapElements ? massmap[z1] : 2.5*z1; }
        
      // get the mean-free-path table for the indexed material 
      const G4_c2_function * operator [] (G4int materialIndex) { 
              return MFPTables.find(materialIndex)!=MFPTables.end() ?
                     &(MFPTables[materialIndex].get()) : (G4_c2_function *)0;
      }
        
protected:
      ScreeningMap screeningData; // screening tables for each element
      ParticleCache targetMap;
      G4int verbosity;
      std::map<G4int, G4_c2_const_ptr > sigmaMap; 
      // total cross section for each element
      std::map<G4int, G4_c2_const_ptr > MFPTables; // MFP for each material
        
private:
      static const G4double massmap[nMassMapElements+1];

};

typedef struct G4CoulombKinematicsInfo {
        G4double impactParameter;
        G4ScreenedCoulombCrossSection *crossSection;
        G4double a1, a2, sinTheta, cosTheta, sinZeta, cosZeta, eRecoil;
        G4ParticleDefinition *recoilIon;        
        const G4Material *targetMaterial;
} G4CoulombKinematicsInfo;

class G4ScreenedCollisionStage {
public:
     virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
                    const class G4Track& aTrack, const class G4Step& aStep)=0;
     virtual ~G4ScreenedCollisionStage() {}
};

class G4ScreenedCoulombClassicalKinematics: 
  public G4ScreenedCoulombCrossSectionInfo, public G4ScreenedCollisionStage {

public:
     G4ScreenedCoulombClassicalKinematics();
     virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
             const class G4Track& aTrack, const class G4Step& aStep);
        
     G4bool DoScreeningComputation(class G4ScreenedNuclearRecoil *master,
             const G4ScreeningTables *screen,  G4double eps, G4double beta);
             
     virtual ~G4ScreenedCoulombClassicalKinematics() { }
     
protected:
     // the c2_functions we need to do the work.
     c2_const_plugin_function_p<G4double> &phifunc;
     c2_linear_p<G4double> &xovereps;
     G4_c2_ptr diff;        
};

class G4SingleScatter: public G4ScreenedCoulombCrossSectionInfo,
                       public G4ScreenedCollisionStage {

public:
     G4SingleScatter() { }
     virtual void DoCollisionStep(class G4ScreenedNuclearRecoil *master,
             const class G4Track& aTrack, const class G4Step& aStep);
     virtual ~G4SingleScatter() {}
};

/**
     \brief A process which handles screened Coulomb collisions between nuclei
 
*/

class G4ScreenedNuclearRecoil : public G4ScreenedCoulombCrossSectionInfo,
                                public G4VDiscreteProcess
{
public:
        
     friend class G4ScreenedCollisionStage;

     /// \brief Construct the process and set some physics parameters for it.
     /// \param processName the name to assign the process
     /// \param ScreeningKey the name of a screening function to use.  
     /// The default functions are "zbl" (recommended for soft scattering),
     /// "lj" (recommended for backscattering) and "mol" (Moliere potential)
     /// \param GenerateRecoils if frue, ions struck by primary are converted 
     /// into new moving particles.
     /// If false, energy is deposited, but no new moving ions are created.
     /// \param RecoilCutoff energy below which no new moving particles will be
     /// created, even if a GenerateRecoils is true.
     /// Also, a moving primary particle will be stopped if its energy falls
     /// below this limit.
     /// \param PhysicsCutoff the energy transfer to which screening tables are
     /// calucalted.  
     /// There is no really
     /// compelling reason to change it from the 10.0 eV default.  
     /// However, see the paper on running this
     /// in thin targets for further discussion, and its interaction 
     /// with SetMFPScaling()
     
     G4ScreenedNuclearRecoil(const G4String& processName = "ScreenedElastic",
                const G4String &ScreeningKey="zbl", G4bool GenerateRecoils=1, 
                G4double RecoilCutoff=100.0*eV, G4double PhysicsCutoff=10.0*eV);
                
     /// \brief destructor
     virtual ~G4ScreenedNuclearRecoil();
     
     /// \brief used internally by Geant4 machinery
     virtual G4double GetMeanFreePath(const G4Track&, G4double,
                                      G4ForceCondition* );
                                      
     /// \brief used internally by Geant4 machinery
     virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
                                             const G4Step& aStep);
     /// \brief test if a prticle of type \a aParticleType can use this process
     /// \param aParticleType the particle to test
     virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
     /// \brief Build physics tables in advance.  Not Implemented.
     /// \param aParticleType the type of particle to build tables for 
     virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
     /// \brief Export physics tables for persistency.  Not Implemented.
     /// \param aParticleType the type of particle to build tables for 
     virtual void DumpPhysicsTable(const G4ParticleDefinition& aParticleType);
     /// \brief deterine if the moving particle is within  the strong force
     /// range of the selected nucleus
     /// \param A the nucleon number of the beam
     /// \param A1 the nucleon number of the target
     /// \param apsis the distance of closest approach
     
     virtual G4bool CheckNuclearCollision(G4double A, G4double A1,
                            G4double apsis); // return true if hard collision

     virtual G4ScreenedCoulombCrossSection *GetNewCrossSectionHandler(void);
        
     /// \brief Get non-ionizing energy loss for last step
     G4double GetNIEL() const { return NIEL; } 
        
     /// \brief clear precomputed screening tables
     void ResetTables(); 
     // clear all data tables to allow changing energy cutoff, materials, etc.

     /// \brief set the upper energy beyond which this process has no 
     /// cross section
     ///
     /// This funciton is used to coordinate this process with G4MSC. 
     /// Typically, G4MSC should 
     /// not be allowed to operate in a range which overlaps that of this 
     /// process.  The criterion which is most reasonable
     /// is that the transition should be somewhere in the modestly 
     /// relativistic regime (500 MeV/u for example).
     /// param energy energy per nucleon for the cutoff
     
     void SetMaxEnergyForScattering(G4double energy) {processMaxEnergy=energy;}
     
     /// \brief find out what screening function we are using
     std::string GetScreeningKey() const { return screeningKey; }
     
     /// \brief enable or disable all energy deposition by this process
     /// \param flag if true, enable deposition of energy (the default). 
     /// If false, disable deposition.
     
     void AllowEnergyDeposition(G4bool flag) {registerDepositedEnergy=flag;}
     
     /// \brief get flag indicating whether deposition is enabled
     G4bool GetAllowEnergyDeposition() const {return registerDepositedEnergy;}
     
     /// \brief enable or disable the generation of recoils.  
     /// If recoils are disabled, the energy they would have received is just
     /// deposited.
     /// param flag if true, create recoil ions in cases in which the energy 
     /// is above the recoilCutoff.  
     /// If false, just deposit the energy.
     
     void EnableRecoils(G4bool flag) { generateRecoils=flag; }
     
     /// \brief find out if generation of recoils is enabled.
     G4bool GetEnableRecoils() const { return generateRecoils; }
     
     /// \brief set the mean free path scaling as specified
     /// \param scale the factor by which the default MFP will be scaled.  
     /// Set to less than 1 for very thin films, typically, 
     /// to sample multiple scattering,
     /// or to greater than 1 for quick simulations with a very long flight path
     
     void SetMFPScaling(G4double scale) { MFPScale=scale; }
     
     /// \brief get the MFPScaling parameter
     G4double GetMFPScaling() const { return MFPScale; }
     
     /// \brief enable or disable whether this process will skip collisions 
     /// which are close enough they need hadronic phsyics.
     /// Default is true (skip close collisions).
     /// Disabling this results in excess nuclear stopping power.
     /// \param flag true results in hard collisions being skipped.  
     /// false allows hard collisions.
     
     void AvoidNuclearReactions(G4bool flag) { avoidReactions=flag; }
     
     /// \brief get the flag indicating whether hadronic collisions are ignored.
     G4bool GetAvoidNuclearReactions() const { return avoidReactions; }
     
     /// \brief set the minimum energy (per nucleon) at which recoils can 
     /// be generated, 
     /// and the energy (per nucleon) below which all ions are stopped. 
     /// \param energy energy per nucleon 
     
     void SetRecoilCutoff(G4double energy) { recoilCutoff=energy; }
     
     /// \brief get the recoil cutoff
     G4double GetRecoilCutoff() const { return recoilCutoff; }
     
     /// \brief set the energy to which screening tables are computed.  
     /// Typically, this is 10 eV or so, and not often changed.
     /// \param energy the cutoff energy 
     
     void SetPhysicsCutoff(G4double energy) { physicsCutoff=energy; 
                                              ResetTables(); }
                                              
     /// \brief get the physics cutoff energy.
     G4double GetPhysicsCutoff() const { return physicsCutoff; }
     
     /// \brief set the pointer to a class for paritioning energy into NIEL
     /// \brief part the pointer to the class.
     
     void SetNIELPartitionFunction(const G4VNIELPartition *part);
     
     /// \brief set the cross section boost to provide faster computation of
     /// backscattering
     /// \param fraction the fraction of particles to have their cross section 
     /// boosted.
     /// \param HardeningFactor the factor by which to boost the scattering 
     /// cross section.
     
     void SetCrossSectionHardening(G4double fraction, G4double HardeningFactor){
                hardeningFraction=fraction;
                hardeningFactor=HardeningFactor;
     }
     
     /// \brief get the fraction of particles which will have boosted scattering
     G4double GetHardeningFraction() const { return hardeningFraction; }
     
     /// \brief get the boost factor in use.
     G4double GetHardeningFactor() const { return hardeningFactor; }
     
     /// \brief the the interaciton length used in the last scattering.
     G4double GetCurrentInteractionLength() const {
                                    return currentInteractionLength; }
                                    
     /// \brief set a function to compute screening tables, 
     /// if the user needs non-standard behavior.
     /// \param cs a class which constructs the screening tables.
     
     void SetExternalCrossSectionHandler(G4ScreenedCoulombCrossSection *cs) { 
                            externalCrossSectionConstructor=cs; 
     }
     
     /// \brief get the verbosity.
     G4int GetVerboseLevel() const { return verboseLevel; }

     std::map<G4int, G4ScreenedCoulombCrossSection*> &GetCrossSectionHandlers() 
                { return crossSectionHandlers; }
                
     void ClearStages(void);
      
     void AddStage(G4ScreenedCollisionStage *stage) { 
                                     collisionStages.push_back(stage); }
                                     
     G4CoulombKinematicsInfo &GetKinematics() { return kinematics; }
     
     void SetValidCollision(G4bool flag) { validCollision=flag; }
     
     G4bool GetValidCollision() const { return validCollision; }

     /// \brief get the pointer to our ParticleChange object. 
     /// for internal use, primarily.
     class G4ParticleChange &GetParticleChange() 
              { return static_cast<G4ParticleChange &>(*pParticleChange); }
              
     /// \brief take the given energy, and use the material information 
     /// to partition it into NIEL and ionizing energy.
     
     void DepositEnergy(G4int z1, G4double a1, const G4Material *material,
                        G4double energy);
        
protected:
     /// \brief the energy per nucleon above which the MFP is constant 
     G4double highEnergyLimit;
     
     /// \brief the energy per nucleon below which the MFP is zero
     G4double lowEnergyLimit;
     
     /// \brief the energy per nucleon beyond which the cross section is zero,
     /// to cross over to G4MSC
     G4double processMaxEnergy;
     G4String screeningKey;
     G4bool generateRecoils, avoidReactions;
     G4double recoilCutoff, physicsCutoff;
     G4bool registerDepositedEnergy;
     G4double IonizingLoss, NIEL;
     G4double MFPScale;
     G4double hardeningFraction, hardeningFactor;
        
     G4ScreenedCoulombCrossSection *externalCrossSectionConstructor;
     std::vector<G4ScreenedCollisionStage *> collisionStages;
        
     std::map<G4int, G4ScreenedCoulombCrossSection*> crossSectionHandlers;
        
     G4bool validCollision;
     G4CoulombKinematicsInfo kinematics;
     const G4VNIELPartition *NIELPartitionFunction;
};

// A customized G4CrossSectionHandler which gets its data from 
// an external program

class G4NativeScreenedCoulombCrossSection: public G4ScreenedCoulombCrossSection 
{
public:
     G4NativeScreenedCoulombCrossSection();

     G4NativeScreenedCoulombCrossSection(
                              const G4NativeScreenedCoulombCrossSection &src) 
          : G4ScreenedCoulombCrossSection(src), phiMap(src.phiMap) { }
        
     G4NativeScreenedCoulombCrossSection(
                                     const G4ScreenedCoulombCrossSection &src) 
          : G4ScreenedCoulombCrossSection(src)  { }
          
     virtual ~G4NativeScreenedCoulombCrossSection();
        
     virtual void LoadData(G4String screeningKey, G4int z1, G4double m1,
                           G4double recoilCutoff);
                           
     virtual G4ScreenedCoulombCrossSection *create() 
                { return new G4NativeScreenedCoulombCrossSection(*this); } 
                
     // get a list of available keys
     std::vector<G4String> GetScreeningKeys() const;
        
     typedef G4_c2_function &(*ScreeningFunc)(G4int z1, G4int z2, 
                                 size_t nPoints, G4double rMax, G4double *au);
        
     void AddScreeningFunction(G4String name, ScreeningFunc fn) {
                           phiMap[name]=fn;
     }
        
private:
     // this is a map used to look up screening function generators
     std::map<std::string, ScreeningFunc> phiMap;
};

G4_c2_function &ZBLScreening(G4int z1, G4int z2, size_t npoints, 
                             G4double rMax, G4double *auval);
                             
G4_c2_function &MoliereScreening(G4int z1, G4int z2, size_t npoints, 
                                 G4double rMax, G4double *auval);
                                 
G4_c2_function &LJScreening(G4int z1, G4int z2, size_t npoints,
                            G4double rMax, G4double *auval);
                                                        
G4_c2_function &LJZBLScreening(G4int z1, G4int z2, size_t npoints,
                               G4double rMax, G4double *auval);

#endif
