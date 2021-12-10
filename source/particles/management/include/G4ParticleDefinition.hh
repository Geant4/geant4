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
// G4ParticleDefinition
//
// Class description:
//
// This class contains all the static data of a particle.
// It uses the process manager in order to collect all the processes
// this kind of particle can undertake.

// Authors: G.Cosmo, 2 December 1995 - Design, based on object model
//          M.Asai, 29 January 1996 - First implementation
// History:
// - 1996-2003, H.Kurashige - Revisions
// - 11.03.2003, H.Kurashige - Restructuring for Cuts per Region
// - 25.01.2013, G.Cosmo, A.Dotti - Introduced thread-safety for MT
// - 15.06.2017, K.L.Genser - Added support for MuonicAtom
// --------------------------------------------------------------------
#ifndef G4ParticleDefinition_hh
#define G4ParticleDefinition_hh 1

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PDefManager.hh"

class G4ProcessManager;
class G4DecayTable;
class G4ParticleTable;
class G4ParticlePropertyTable;
class G4VTrackingManager;

using G4ParticleDefinitionSubInstanceManager = G4PDefManager;

class G4ParticleDefinition 
{
  friend class G4ParticlePropertyTable;

  public:

    // Only one type of constructor can be used for G4ParticleDefinition.
    // If you want to create new particle, you must set name of the particle 
    // at construction. Most of members seen as arguments of the constructor 
    // (except last 3 arguments concerning with decay ) are  "constant" 
    // and can not be changed later. (No "SET" methods are available)
    // Each type of particle must be constructed as a unique object
    // of special class derived from G4ParticleDefinition.
    // See G4ParticleTypes for detail
 
    G4ParticleDefinition(const G4String&  aName,  
                         G4double         mass,     
                         G4double         width,
                         G4double         charge,   
                         G4int            iSpin,
                         G4int            iParity,
                         G4int            iConjugation,
                         G4int            iIsospin,   
                         G4int            iIsospinZ, 
                         G4int            gParity,
                         const G4String&  pType,
                         G4int            lepton,
                         G4int            baryon,
                         G4int            encoding,
                         G4bool           stable,
                         G4double         lifetime,
                         G4DecayTable*    decaytable,
                         G4bool           shortlived = false,
                         const G4String&  subType = "",
                         G4int            anti_encoding = 0,
                         G4double         magneticMoment = 0.0);

    virtual ~G4ParticleDefinition();
      
    G4ParticleDefinition(const G4ParticleDefinition&) = delete;
    G4ParticleDefinition& operator=(const G4ParticleDefinition &) = delete;
      // Can not use "copy constructor", equality nor "default constructor"!

    G4bool operator==(const G4ParticleDefinition& right) const;
    G4bool operator!=(const G4ParticleDefinition& right) const;

    // With the following Getxxxx methods, one can get values 
    // for members which can not be changed
 
    const G4String& GetParticleName() const { return theParticleName; }

    G4double GetPDGMass() const { return thePDGMass; }
    G4double GetPDGWidth() const { return thePDGWidth; } 
    G4double GetPDGCharge() const { return thePDGCharge; }

    G4double GetPDGSpin() const { return thePDGSpin; }
    G4int    GetPDGiSpin() const { return thePDGiSpin; }
    G4int    GetPDGiParity() const { return thePDGiParity; }
    G4int    GetPDGiConjugation() const { return thePDGiConjugation; }
    G4double GetPDGIsospin() const { return thePDGIsospin; }
    G4double GetPDGIsospin3() const { return thePDGIsospin3; }
    G4int    GetPDGiIsospin() const { return thePDGiIsospin; }
    G4int    GetPDGiIsospin3() const { return thePDGiIsospin3; }
    G4int    GetPDGiGParity() const { return thePDGiGParity; }
 
    G4double GetPDGMagneticMoment() const { return thePDGMagneticMoment; }
    inline void SetPDGMagneticMoment(G4double mageticMoment); 
    G4double CalculateAnomaly()  const;
      // Gives the anomaly of magnetic moment for spin 1/2 particles 

    const G4String& GetParticleType() const { return theParticleType; }
    const G4String& GetParticleSubType() const { return theParticleSubType; }
    G4int GetLeptonNumber() const { return theLeptonNumber; }
    G4int GetBaryonNumber() const { return theBaryonNumber; }

    G4int GetPDGEncoding() const { return thePDGEncoding; }
    G4int GetAntiPDGEncoding() const { return theAntiPDGEncoding; }
    inline void SetAntiPDGEncoding(G4int aEncoding);

    inline G4int GetQuarkContent(G4int flavor) const;
    inline G4int GetAntiQuarkContent(G4int flavor) const;
      // Returns the number of quark with flavor contained in this particle. 
      // The value of flavor is assigned as follows 
      // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t
 
    G4bool IsShortLived() const { return fShortLivedFlag; }

    inline G4bool GetPDGStable() const; 
    void SetPDGStable(const G4bool aFlag) { thePDGStable=aFlag; }

    inline G4double GetPDGLifeTime() const;
    void SetPDGLifeTime(G4double aLifeTime) { thePDGLifeTime=aLifeTime; }

    inline G4double GetIonLifeTime() const;
      // Get life time of a generic ion through G4NuclideTable.

    inline G4DecayTable* GetDecayTable() const;
    inline void SetDecayTable(G4DecayTable* aDecayTable); 
      // Set/Get Decay Table
      //   !! Decay Table can be modified !!  

    G4ProcessManager* GetProcessManager() const; 
    void SetProcessManager(G4ProcessManager* aProcessManager); 
      // Set/Get Process Manager
      //   !! Process Manager can be modified !!  

    G4VTrackingManager* GetTrackingManager() const;
    void SetTrackingManager(G4VTrackingManager* aTrackingManager);
      // Set/Get Tracking Manager; nullptr means the default
      //   !! Tracking Manager can be modified !!

    inline G4ParticleTable* GetParticleTable() const;
      // Get pointer to the particle table

    inline G4int GetAtomicNumber() const;
    inline G4int GetAtomicMass() const;
      // Get AtomicNumber and AtomicMass
      // These properties are defined for nucleus

    void DumpTable() const;
      // Prints information of data members.

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // Control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

    void SetApplyCutsFlag(G4bool);
    inline G4bool GetApplyCutsFlag() const;

    inline G4bool IsGeneralIon() const;
      // True only if the particle is G4Ions
      // (it means that theProcessManager is same as one for G4GenricIon)

    inline G4bool IsMuonicAtom() const;
      // True only if the particle is a G4MuonicAtom
      // (it means that theProcessManager is same as the one for G4MuonicAtom)

    inline G4ProcessManager* GetMasterProcessManager() const;
      // Returns the process manager master pointer.
    inline void SetMasterProcessManager(G4ProcessManager* aNewPM);
      // Sets the shadow master pointer (not to be used by user)

    inline G4int GetInstanceID() const;
      // Returns the instance ID

    static const G4PDefManager& GetSubInstanceManager();
      // Returns the private data instance manager

    static void Clean();
      // Clear memory allocated by sub-instance manager

    void SetParticleDefinitionID(G4int id=-1);
    inline G4int GetParticleDefinitionID() const;

    inline G4bool IsHypernucleus() const;
    inline G4int GetNumberOfLambdasInHypernucleus() const;
    inline G4bool IsAntiHypernucleus() const;
    inline G4int GetNumberOfAntiLambdasInAntiHypernucleus() const;
      // The first two methods return "false" and 0, respectively,
      // if the particle is not an hypernucleus; else, they return
      // "true" and the number of Lambdas bound in the nucleus.
      // Similarly, the last two methods return "false" and 0,
      // respectively, if the particle is not an anti-hypernucleus;
      // else, they return "true" and the number of anti-Lambdas
      // bound in the anti-nucleus.
      // Notice that, for the time being, we are assuming that
      // (anti-)Lambda is the only type of (anti-)hyperon present
      // in all (anti-)hypernuclei.

  protected:

    G4ParticleDefinition();
      // Cannot be used

    G4int FillQuarkContents();
      // Calculates quark and anti-quark contents
      // return value is the PDG encoding for this particle.
      // It means error if the return value is different from
      // this->thePDGEncoding.

    inline void SetParticleSubType(const G4String& subtype);

    inline void SetAtomicNumber(G4int);
    inline void SetAtomicMass(G4int);

    enum { NumberOfQuarkFlavor = 6 };

    G4int  theQuarkContent[NumberOfQuarkFlavor];
    G4int  theAntiQuarkContent[NumberOfQuarkFlavor];
      //  the number of quark (minus Sign means anti-quark) contents
      //  The value of flavor is assigned as follows 
      //    0:d, 1:u, 2:s, 3:c, 4:b, 5:t

  protected:

    G4bool isGeneralIon = false;
    G4bool isMuonicAtom = false;

  private:

    // --- Shadow of master pointers

    G4ProcessManager* theProcessManagerShadow = nullptr;
      // Each worker thread can access this field from the master thread
      // through this pointer.

    G4int g4particleDefinitionInstanceID = 0;
      // This field is used as instance ID.

    G4PART_DLL static G4PDefManager subInstanceManager;
      // This field helps to use the class G4PDefManager introduced above.

    //  --- Following values can not be changed
    //  --- i.e. No Setxxxx Methods for them 

    G4String theParticleName = "";
      // The name of the particle.
      // Each object must have its specific name!!

    //  --- Following member values must be defined with Units

    G4double thePDGMass = 0.0;
      // The mass of the particle, in units of equivalent energy.

    G4double thePDGWidth = 0.0;
      // The decay width of the particle, usually the width of a
      // Breit-Wigner function, assuming that you are near the
      // mass center anyway. (in units of equivalent energy)

    G4double thePDGCharge = 0.0;
      // The charge of the particle.(in units of Coulomb)

    //   --- Following members are quantum number
    //       i.e. discrete numbers can be allowed
    //       So, you can define them only by using integer in constructor 

    G4int thePDGiSpin = 0;
      // The total spin of the particle, also often denoted as
      // capital J, in units of 1/2.
    G4double thePDGSpin = 0.0;
      // The total spin of the particle, in units of 1.

    G4int thePDGiParity = 0;
      // The parity quantum number, in units of 1. If the parity
      // is not defined for this particle, we will set this to 0.

    G4int thePDGiConjugation = 0;
      // This charge conjugation quantum number in units of 1.

    G4int thePDGiGParity = 0;
      // The value of the G-parity quantum number.

    G4int thePDGiIsospin = 0;
    G4int thePDGiIsospin3 = 0;
      // The isospin and its 3rd-component in units of 1/2.
    G4double thePDGIsospin = 0.0;
    G4double thePDGIsospin3 = 0.0;
      // The isospin quantum number in units of 1.
 
    G4double thePDGMagneticMoment = 0.0;
      // The magnetic moment.

    G4int theLeptonNumber = 0;
      // The lepton quantum number.

    G4int theBaryonNumber = 0;
      // The baryon quantum number.

    G4String theParticleType = "";
      // More general textual type description of the particle.

    G4String theParticleSubType = "";
      // Textual type description of the particle
      // eg. pion, lamda etc.

    G4int thePDGEncoding = 0;
      // The Particle Data Group integer identifier of this particle
 
    G4int theAntiPDGEncoding = 0;
      // The Particle Data Group integer identifier of the anti-particle

    // --- Following members can be changed after construction

    G4bool fShortLivedFlag = false;
      // Particles which have true value of this flag 
      // will not be tracked by TrackingManager 

    G4bool thePDGStable = false;
      // Is an indicator that this particle is stable. It must
      // not decay. If the user tries to assign a kind of decay
      // object to it, it will refuse to take it.

    G4double thePDGLifeTime = 0.0;
      // Is related to the decay width of the particle. The mean
      // life time is given in seconds.

    G4DecayTable* theDecayTable = nullptr;
      // Points DecayTable 

    G4ParticleTable* theParticleTable = nullptr;

    G4int theAtomicNumber = 0;
    G4int theAtomicMass = 0;
 
    G4int verboseLevel = 1;
    G4bool fApplyCutsFlag = false;
};

#include "G4ParticleDefinition.icc"

#endif
