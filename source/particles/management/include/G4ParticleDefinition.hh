// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleDefinition.hh,v 1.6 1999-12-15 14:51:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ---------------- G4ParticleDefinition ----------------
// first implementation by Makoto Asai, 29 January 1996
// revised by G.Cosmo, 29 February 1996
// revised by H.Kurashige, 19 April 1996
// revised by H.Kurashige, 4 July 1996
// added GetEnergyCuts() and GetLengthCuts() by G.Cosmo, 11 July 1996
// added Set/GetVerboseLevel()    H.Kurashige 11 Nov. 1997
// added SetCuts() and ResetCuts  H.Kurashige 15 Nov.1996
// change SetProcessManager as public H.Kurashige 06 June 1998
// added  GetEnergyThreshold  H.Kurashige 08 June 1998
// added  ShortLived flag and ApplyCuts flag  H.Kurashige 27  June 1998
// fixed  some improper codings   H.Kurashige 08 Apr. 1999
// ------------------------------------------------------------

#ifndef G4ParticleDefinition_h
#define G4ParticleDefinition_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ProcessManager;
class G4DecayTable;
class G4ParticleTable;
class G4Material;

class G4ParticleDefinition 
{
  // Class Description
  //  This class containes all the static data of a particle.
  //  It also has uses a process manager in order to collect
  //  all the processes this kind of particle can undertake.
  //

 public: // With Description
  //  Only one type of constructor can be used for G4ParticleDefinition.
  //  If you want to create new particle, you must set name of the particle 
  //  at construction. Most of members seen as arguments of the constructor 
  //  (except last 3 arguments concerning with decay ) are  "constant" 
  //  and can not be changed later. (No "SET" methods are available)
  //  Each type of particle must be constructed as a unique static object
  //  of special class derived from G4ParticleDefinition.
  //  see G4ParticleTypes for detail 
 
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
			   G4bool           shortlived = false);

       virtual ~G4ParticleDefinition();
      
  public: // With Description
      // These methods concerning cut values are provided 
      // to invoke corresponding methods for each particle type.
      // Actual implementation can be seen in the class 
      // G4ParticleWithCuts  
      virtual void          ResetCuts();
      virtual void          SetCuts(G4double );
      virtual void          ReCalcCuts();
      virtual G4double      	GetLengthCuts() const;
      virtual G4double*	        GetEnergyCuts() const;
      virtual G4double      	GetEnergyThreshold(const G4Material* ) const;
 
      // applyCuts flag
      G4bool                GetApplyCutsFlag() const;
      void                  SetApplyCutsFlag(G4bool flag);
      
  public: // With Description
  // By these following Getxxxx methods, you can get values 
  // for members which can not be changed
  //  
      const G4String& GetParticleName() const { return theParticleName; }

      G4double GetPDGMass() const { return thePDGMass; }
      G4double GetPDGWidth() const { return thePDGWidth; } 
      G4double GetPDGCharge() const { return thePDGCharge; }

      G4double GetPDGSpin() const { return thePDGSpin; }
      G4int    GetPDGiSpin() const { return thePDGiSpin; }
      G4int    GetPDGiParity() const { return thePDGiParity; }
      G4int    GetPDGiConjugation() const { return thePDGiConjugation; }
      G4double GetPDGIsospin() const { return thePDGiIsospin; }
      G4double GetPDGIsospin3() const { return thePDGiIsospin3; }
      G4int    GetPDGiIsospin() const { return thePDGiIsospin; }
      G4int    GetPDGiIsospin3() const { return thePDGiIsospin3; }
      G4int    GetPDGiGParity() const { return thePDGiGParity; }

      G4String GetParticleType() const { return theParticleType; }
      G4int    GetLeptonNumber() const { return theLeptonNumber; }
      G4int    GetBaryonNumber() const { return theBaryonNumber; }

      G4int    GetPDGEncoding() const { return thePDGEncoding; }
      G4int    GetAntiPDGEncoding() const { return theAntiPDGEncoding; }
      void     SetAntiPDGEncoding(G4int aEncoding);

 
      G4int    GetQuarkContent(G4int flavor) const;
      G4int    GetAntiQuarkContent(G4int flavor) const;
      //  return the number of quark with flavor contained in this particle. 
      //  The value of flavor is assigned as follows 
      // 1:d, 2:u, 3:s, 4:c, 5:b, 6:t, 7:l(down type quark) 8:h(up type quark) 
 
  public: // With Description
      // ShortLived flag
      G4bool   IsShortLived() const { return fShortLivedFlag; }

      G4bool   GetPDGStable() const { return thePDGStable; }
      void     SetPDGStable(const G4bool aFlag) { thePDGStable=aFlag; }

      G4double GetPDGLifeTime() const { return thePDGLifeTime; }
      void     SetPDGLifeTime(G4double aLifeTime) { thePDGLifeTime = aLifeTime; }

  public:// With Description
      G4DecayTable* GetDecayTable();
      void          SetDecayTable(G4DecayTable* aDecayTable); 
      // Set/Get Decay Table
      //   !! Decay Table can be modified !!  

  public: // With Description
      G4ProcessManager* GetProcessManager() const; 
      void SetProcessManager(G4ProcessManager* aProcessManager); 
      // Set/Get Process Manager
      //   !! Process Manager can be modified !!  

      G4ParticleTable* GetParticleTable();
      // get pointer to the particle table

      void DumpTable() const;
      //  Prints information of data members.

  protected:
      G4int   FillQuarkContents();
      //  calculate quark and anti-quark contents
      //  return value is PDG encoding for this particle.
      //  It means error if the return value is deffernt from
      //  this->thePDGEncoding.

 public:
      void  SetVerboseLevel(G4int value);
      G4int GetVerboseLevel() const;
      // controle flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  private:
  //  !!!  can not use "copy constructor" nor "default constructor" !!!!
       G4ParticleDefinition(const G4ParticleDefinition &right);
       G4ParticleDefinition();

  private:  
  //  !!!  Assignment operation is forbidden !!!
      const G4ParticleDefinition & operator=(const G4ParticleDefinition &right);

  public:
      G4int operator==(const G4ParticleDefinition &right) const;
      G4int operator!=(const G4ParticleDefinition &right) const;

  private:
  //  Values following can not be changed
  //  i.e. No Setxxxx Methods for them 

      G4String theParticleName;
      //  The name of the particle.
      //  Each object must have its specific name!!

    //    --- following member values must be defined with Units
      G4double thePDGMass;
      //  The mass of the particle, in units of equivalent energy.

      G4double thePDGWidth;
      //  The decay width of the particle, usually the width of a
      //  Breit-Wigner function, assuming that you are near the
      //  mass center anyway. (in units of equivalent energy)

      G4double thePDGCharge;
      //  The charge of the particle.(in units of Coulomb)

    //   ---- following members are quantum number
    //         i.e. discrete numbers can be allowded
    //        So, you can defined only by using integer in constructor 
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

      G4int thePDGiIsospin;
      G4int thePDGiIsospin3;
      //  The isospin and its 3rd-component in units of 1/2.
      G4double thePDGIsospin;
      G4double thePDGIsospin3;
      //  The isospin quantum number in units of 1.

      G4int thePDGiGParity;
      //  The value of the G-parity quantum number.

      G4int theLeptonNumber;
      //  The lepton quantum number.

      G4int theBaryonNumber;
      //  The baryon quantum number.

      G4String theParticleType;
      //  More general textual type description of the particle.

      G4int thePDGEncoding;
      //  The Particle Data Group integer identifier of this particle
 
      G4int theAntiPDGEncoding;
      //  The Particle Data Group integer identifier of the anti-particle

  protected:
      enum {NumberOfQuarkFlavor = 8};
      G4int  theQuarkContent[NumberOfQuarkFlavor];
      G4int  theAntiQuarkContent[NumberOfQuarkFlavor];
      //  the number of quark (minus Sign means anti-quark) contents
      //  The value of flavor is assigned as follows 
      //    0:d, 1:u, 2:s, 3:c, 
      //    4:b, 5:t, 6:l(down type quark) 7:h(up type quark) 


  private:
    // Following members can be changed after construction

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

      class G4DecayTable *theDecayTable;
      //  Points DecayTable 
 
   private:
      class G4ProcessManager *theProcessManager;
      //  Points to G4ProcessManager

      G4ParticleTable* theParticleTable;
 
      // ApplyCuts flag 
      G4bool                fApplyCutsFlag;
     
  
 private:
   G4int verboseLevel;
};

inline 
G4ParticleTable* G4ParticleDefinition::GetParticleTable()
{ 
  return theParticleTable; 
}

inline
G4DecayTable* G4ParticleDefinition::GetDecayTable()
{ 
  return theDecayTable; 
}

inline
void          G4ParticleDefinition::SetDecayTable(G4DecayTable* aDecayTable) 
{ 
  theDecayTable = aDecayTable; 
}

inline 
void G4ParticleDefinition::SetApplyCutsFlag(G4bool flag)
{
  fApplyCutsFlag = flag;
}

inline 
G4bool G4ParticleDefinition::GetApplyCutsFlag() const
{
  return fApplyCutsFlag;
}

inline 
void G4ParticleDefinition::SetVerboseLevel(G4int value)
{
   verboseLevel = value;
}

inline 
G4int G4ParticleDefinition::GetVerboseLevel() const
{
   return verboseLevel;
}

inline
G4ProcessManager* G4ParticleDefinition::GetProcessManager() const
{
    return theProcessManager;
}

inline
void G4ParticleDefinition::SetProcessManager(G4ProcessManager *aProcessManager)
{
      theProcessManager = aProcessManager;
}

inline
G4int G4ParticleDefinition::GetQuarkContent(G4int flavor) const 
{ 
  if ((flavor>0) && (flavor<NumberOfQuarkFlavor)){
    return theQuarkContent[flavor-1];
  }else {
    if (verboseLevel >0) {
      G4cout << "Invalid Quark Flavor for G4ParticleDefinition::GetQuarkContent";
      G4cout << ": flavor=" << flavor <<G4endl;
    }
    return 0;
  }  
}

inline
G4int G4ParticleDefinition::GetAntiQuarkContent(G4int flavor) const 
{ 
  if ((flavor>0) && (flavor<NumberOfQuarkFlavor)){
    return theAntiQuarkContent[flavor-1];
  }else {
    if (verboseLevel >0) {
      G4cout <<"Invalid Quark Flavor for G4ParticleDefinition::GetAntiQuarkContent";
      G4cout << ": flavor=" << flavor <<G4endl;
    }
    return 0;
  }  
}

inline void             G4ParticleDefinition::ResetCuts(){}

inline void             G4ParticleDefinition::SetCuts(G4double ){}

inline void             G4ParticleDefinition::ReCalcCuts(){}

inline G4double         G4ParticleDefinition::GetLengthCuts() const {return -1.0;}

inline G4double*	G4ParticleDefinition::GetEnergyCuts() const {return 0;}

inline G4double      	G4ParticleDefinition::GetEnergyThreshold(const G4Material* ) const     
{      
  return -1.0 * eV;   
}
 
inline void             G4ParticleDefinition::SetAntiPDGEncoding(G4int aEncoding)
{ 
  theAntiPDGEncoding = aEncoding; 
}

#endif





















