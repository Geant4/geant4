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
// $Id: G4ParticleDefinition.hh 103108 2017-03-16 13:00:35Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ---------------- G4ParticleDefinition ----------------
// first implementation by Makoto Asai - 29 January 1996
// revised - G.Cosmo - 29 February 1996
// revised - H.Kurashige - 19 April 1996
// revised -  H.Kurashige - 4 July 1996
// added GetEnergyCuts() and GetLengthCuts() - G.Cosmo - 11 July 1996
// added Set/GetVerboseLevel() - H.Kurashige - 11 November 1997
// added SetCuts() and ResetCuts - H.Kurashige - 15 November 1996
// change SetProcessManager as public - H.Kurashige - 06 June 1998
// added  GetEnergyThreshold - H.Kurashige - 08 June 1998
// added  ShortLived flag and ApplyCuts flag - H.Kurashige - 27 June 1998
// fixed  some improper codings - H.Kurashige - 08 April 1999
// added  sub-type - H.Kurashige - 15 February 2000
// added  RestoreCuts - H.Kurashige - 09 March 2001
// restructuring for Cuts per Region - H.Kurashige - 11 March 2003 
// added  MagneticMoment - H.Kurashige - March 2007
// modified for thread-safety for MT - G.Cosmo, A.Dotti - January 2013
// ------------------------------------------------------------

#ifndef G4ParticleDefinition_h
#define G4ParticleDefinition_h 1

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PDefManager.hh"

class G4ProcessManager;
class G4DecayTable;
class G4ParticleTable;
class G4ParticlePropertyTable;

//typedef G4PDefSplitter<G4PDefData>  G4PDefManager;
typedef G4PDefManager G4ParticleDefinitionSubInstanceManager;

class G4ParticleDefinition 
{
  // Class Description:
  //
  // This class containes all the static data of a particle.
  // It also has uses a process manager in order to collect
  // all the processes this kind of particle can undertake.

  friend class  G4ParticlePropertyTable;

 public: // With Description

  // Only one type of constructor can be used for G4ParticleDefinition.
  // If you want to create new particle, you must set name of the particle 
  // at construction. Most of members seen as arguments of the constructor 
  // (except last 3 arguments concerning with decay ) are  "constant" 
  // and can not be changed later. (No "SET" methods are available)
  // Each type of particle must be constructed as a unique object
  // of special class derived from G4ParticleDefinition.
  // see G4ParticleTypes for detail 
 
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
			   G4DecayTable     *decaytable,
			   G4bool           shortlived = false,
                           const G4String&  subType ="",
                           G4int            anti_encoding =0,
			   G4double         magneticMoment = 0.0);

      virtual ~G4ParticleDefinition();
      
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
      void     SetPDGMagneticMoment(G4double mageticMoment); 
      G4double CalculateAnomaly()  const;
        // Gives the anomaly of magnetic moment for spin 1/2 particles 

      const G4String& GetParticleType() const { return theParticleType; }
      const G4String& GetParticleSubType() const { return theParticleSubType; }
      G4int    GetLeptonNumber() const { return theLeptonNumber; }
      G4int    GetBaryonNumber() const { return theBaryonNumber; }

      G4int    GetPDGEncoding() const { return thePDGEncoding; }
      G4int    GetAntiPDGEncoding() const { return theAntiPDGEncoding; }
      void     SetAntiPDGEncoding(G4int aEncoding);

 
      G4int    GetQuarkContent(G4int flavor) const;
      G4int    GetAntiQuarkContent(G4int flavor) const;
        // Returns the number of quark with flavor contained in this particle. 
        // The value of flavor is assigned as follows 
        // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t
 
      G4bool   IsShortLived() const { return fShortLivedFlag; }

      G4bool   GetPDGStable() const; 
      void     SetPDGStable(const G4bool aFlag) { thePDGStable=aFlag; }

      G4double GetPDGLifeTime() const;
      void     SetPDGLifeTime(G4double aLifeTime) { thePDGLifeTime=aLifeTime; }

      G4double GetIonLifeTime() const;
        // Get life time of a generic ion through G4NuclideTable.

      G4DecayTable* GetDecayTable() const;
      void          SetDecayTable(G4DecayTable* aDecayTable); 
        // Set/Get Decay Table
        //   !! Decay Table can be modified !!  

      G4ProcessManager* GetProcessManager() const; 
      void SetProcessManager(G4ProcessManager* aProcessManager); 
        // Set/Get Process Manager
        //   !! Process Manager can be modified !!  

      G4ParticleTable* GetParticleTable() const;
        // Get pointer to the particle table

      G4int GetAtomicNumber() const;
      G4int GetAtomicMass() const;
        // Get AtomicNumber and AtomicMass
        // These properties are defined for nucleus

      void DumpTable() const;
        //  Prints information of data members.

      void  SetVerboseLevel(G4int value);
      G4int GetVerboseLevel() const;
        // controle flag for output message
        //  0: Silent
        //  1: Warning message
        //  2: More

      void   SetApplyCutsFlag(G4bool);
      G4bool GetApplyCutsFlag() const;

      G4bool IsGeneralIon() const;
      // true only if the particle is G4Ions
      // (it means that theProcessManager is same as one for G4GenricIon)

      G4int operator==(const G4ParticleDefinition &right) const;
      G4int operator!=(const G4ParticleDefinition &right) const;

  public :  // without description

      inline G4ProcessManager* GetMasterProcessManager() const;
      // Returns the process manager master pointer.
      inline void SetMasterProcessManager(G4ProcessManager* aNewPM);
      // Sets the shadow master pointer (not to be used by user)

      inline G4int GetInstanceID() const;
      // Returns the instance ID.

      static const G4PDefManager& GetSubInstanceManager();
      // Returns the private data instance manager.

      static void Clean();
      // Clear memory allocated by sub-instance manager.

 private:
      // --- Shadow of master pointers.

      G4ProcessManager *theProcessManagerShadow;
      // Each worker thread can access this field from the master thread
      // through this pointer.

      G4int g4particleDefinitionInstanceID;
      // This field is used as instance ID.

      G4PART_DLL static G4PDefManager subInstanceManager;
      // This field helps to use the class G4PDefManager introduced above.

  protected:

      G4int FillQuarkContents();
        // Calculates quark and anti-quark contents
        // return value is PDG encoding for this particle.
        // It means error if the return value is deffernt from
        // this->thePDGEncoding.

      void SetParticleSubType(const G4String& subtype);

      void SetAtomicNumber(G4int );
      void SetAtomicMass(G4int );

      //  !!!  can not use "copy constructor" nor "default constructor" !!!!
      //
      G4ParticleDefinition(const G4ParticleDefinition &right);
      G4ParticleDefinition();

  private:

      //  !!!  Assignment operation is forbidden !!!
      //
      const G4ParticleDefinition & operator=(const G4ParticleDefinition &r);

  protected:

      enum {NumberOfQuarkFlavor = 6};
      G4int  theQuarkContent[NumberOfQuarkFlavor];
      G4int  theAntiQuarkContent[NumberOfQuarkFlavor];
      //  the number of quark (minus Sign means anti-quark) contents
      //  The value of flavor is assigned as follows 
      //    0:d, 1:u, 2:s, 3:c, 4:b, 5:t

  private:

      //  --- Following values can not be changed
      //  --- i.e. No Setxxxx Methods for them 

      G4String theParticleName;
      //  The name of the particle.
      //  Each object must have its specific name!!

      //  --- Following member values must be defined with Units

      G4double thePDGMass;
      //  The mass of the particle, in units of equivalent energy.

      G4double thePDGWidth;
      //  The decay width of the particle, usually the width of a
      //  Breit-Wigner function, assuming that you are near the
      //  mass center anyway. (in units of equivalent energy)

      G4double thePDGCharge;
      //  The charge of the particle.(in units of Coulomb)

      //   --- Following members are quantum number
      //       i.e. discrete numbers can be allowded
      //       So, you can defined only by using integer in constructor 

      G4int thePDGiSpin;
      //  The total spin of the particle, also often denoted as
      //  capital J, in units of 1/2.
      G4double thePDGSpin;
      //  The total spin of the particle, in units of 1.

      G4int thePDGiParity;
      //  The parity quantum number, in units of 1. If the parity
      //  is not defined for this particle, we will set this to 0.

      G4int thePDGiConjugation;
      //  This charge conjugation quantum number in units of 1.

      G4int thePDGiGParity;
      //  The value of the G-parity quantum number.

      G4int thePDGiIsospin;
      G4int thePDGiIsospin3;
      //  The isospin and its 3rd-component in units of 1/2.
      G4double thePDGIsospin;
      G4double thePDGIsospin3;
      //  The isospin quantum number in units of 1.
 
      G4double thePDGMagneticMoment;
      //  The magnetic moment.

      G4int theLeptonNumber;
      //  The lepton quantum number.

      G4int theBaryonNumber;
      //  The baryon quantum number.

      G4String theParticleType;
      //  More general textual type description of the particle.

      G4String theParticleSubType;
      // Textual type description of the particle
      // eg. pion, lamda etc.

      G4int thePDGEncoding;
      //  The Particle Data Group integer identifier of this particle
 
      G4int theAntiPDGEncoding;
      //  The Particle Data Group integer identifier of the anti-particle

      // --- Following members can be changed after construction

      G4bool fShortLivedFlag;
      //  Particles which have true value of this flag 
      //  will not be tracked by TrackingManager 

      G4bool thePDGStable;
      //  Is an indicator that this particle is stable. It must
      //  not decay. If the user tries to assign a kind of decay
      //  object to it, it will refuse to take it.

      G4double thePDGLifeTime;
      //  Is related to the decay width of the particle. The mean
      //  life time is given in seconds.

      G4DecayTable *theDecayTable;
      //  Points DecayTable 

   private:

      G4ParticleTable* theParticleTable;

      G4int theAtomicNumber;
      G4int theAtomicMass;
 
      G4int verboseLevel;
      G4bool fApplyCutsFlag;

   protected:
      G4bool isGeneralIon;

   public:
      void SetParticleDefinitionID(G4int id=-1);
      G4int GetParticleDefinitionID() const;
};

#include "G4ParticleDefinition.icc"

#endif
