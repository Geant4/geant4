// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDecayChannel.hh,v 1.3 1999-08-30 08:27:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May  1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4VDecayChannel_h
#define G4VDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"

class    G4ParticleDefinition;
class    G4DecayProducts;
class    G4ParticleTable;

class G4VDecayChannel
{
  // This class is a abstract class
  public:
    //Constructors 
      G4VDecayChannel(const G4String &aName, G4int Verbose = 1);
      G4VDecayChannel(const G4String  &aName, 
		     const G4String& theParentName,
		     G4double        theBR,
		     G4int           theNumberOfDaughters,
		     const G4String& theDaughterName1,
		     const G4String& theDaughterName2 = "",
		     const G4String& theDaughterName3 = "",
		     const G4String& theDaughterName4 = "" );

    //  Destructor
      virtual ~G4VDecayChannel();

  private:
    //  copy constructor and assignment operatotr
     G4VDecayChannel(const G4VDecayChannel &);
     G4VDecayChannel & operator=(const G4VDecayChannel &);

  public:
    // equality operators
    G4int operator==(const G4VDecayChannel &right) const {return (this == &right);};
    G4int operator!=(const G4VDecayChannel &right) const {return (this != &right);};

    // less-than operator is defined for G4DecayTable
    G4int operator<(const G4VDecayChannel &right) const;

  public:
   virtual G4DecayProducts* DecayIt(G4double parentMass = -1.0) = 0;

  public:    
     //get kinematics name
     G4String  GetKinematicsName() const;
     //get branching ratio
     G4double   GetBR() const;
     //get number of daughter particles
     G4int      GetNumberOfDaughters() const;     

     //get the pointer to the parent particle
     G4ParticleDefinition * GetParent();
     //get the pointer to a daughter particle 
     G4ParticleDefinition * GetDaughter(G4int anIndex);

     //get the name of the parent particle
     const G4String& GetParentName() const;
     //get the name of a daughter particle
     const G4String& GetDaughterName(G4int anIndex) const;

     // get mass of parent
     G4double GetParentMass() const; 
     G4double GetDaughterMass(G4int anIndex) const; 

     //set the parent particle (by name or by pointer) 
     void SetParent(const G4ParticleDefinition * particle_type);
     void SetParent(const G4String &particle_name);
     //set branching ratio
     void SetBR(G4double value); 
     //set number of daughter particles
     void SetNumberOfDaughters(G4int value);     
     //set a daughter particle (by name or by pointer) 
     void SetDaughter(G4int anIndex, 
		      const G4ParticleDefinition * particle_type);
     void SetDaughter(G4int anIndex, 
		      const G4String &particle_name);

  protected:
    // kinematics name
    G4String   kinematics_name;
    // branching ratio  [0.0 - 1.0]
    G4double   rbranch;
    // number of daughters
    G4int      numberOfDaughters;
    // parent particle
    G4String*  parent_name;
    //daughter particles
    G4String** daughters_name;

  protected:
    // celar daughters array
    void ClearDaughtersName();

  protected:
    // pointer to particle table
    G4ParticleTable*       particletable;

    // temporary buffers of pointers to G4ParticleDefinition
    G4ParticleDefinition*  parent;
    G4ParticleDefinition** daughters;

    // parent mass
    G4double               parent_mass;
    G4double*              daughters_mass;
    

    // fill daughters array
    void FillDaughters();
    // fill parent
    void FillParent();

  public:
    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel()  const;
    void  DumpInfo();

  private:  
    // controle flag for output message
    G4int verboseLevel;
    //  0: Silent
    //  1: Warning message
    //  2: More

    static const G4String   noName;
};

inline
 G4int G4VDecayChannel::operator<(const G4VDecayChannel &right) const
{
  return (this->rbranch < right.rbranch);
}

inline 
  G4ParticleDefinition* G4VDecayChannel::GetDaughter(G4int anIndex)
 { 
  //pointers to daughter particles are filled, if they are not set yet 
  if (daughters == 0) FillDaughters();

  //get the pointer to a daughter particle
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return daughters[anIndex];
  } else {
    if (verboseLevel>0)
      G4cout << "G4VDecayChannel::GetDaughter  index out of range "<<anIndex<<endl;
    return 0;
  }
}

inline
 const G4String& G4VDecayChannel::GetDaughterName(G4int anIndex) const
{
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return *daughters_name[anIndex];
  } else {
    if (verboseLevel>0){
      G4cout << "G4VDecayChannel::GetDaughterName ";
      G4cout << "index out of range " << anIndex << endl;
    }
    return noName;
  }
}

inline
 G4double G4VDecayChannel::GetDaughterMass(G4int anIndex) const
{
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return daughters_mass[anIndex];
  } else {
    if (verboseLevel>0){
      G4cout << "G4VDecayChannel::GetDaughterMass ";
      G4cout << "index out of range " << anIndex << endl;
    }
    return 0.0;
  }
}

inline 
  G4ParticleDefinition* G4VDecayChannel::GetParent()
{ 
  //the pointer to the parent particle is filled, if it is not set yet 
   if (parent == 0) FillParent();
  //get the pointer to the parent particle
  return parent;
}

inline
 const G4String& G4VDecayChannel::GetParentName() const
{
  return *parent_name;
}

inline
 G4double G4VDecayChannel::GetParentMass() const
{
  return parent_mass;
}


inline
  void G4VDecayChannel::SetParent(const G4String &particle_name)
{
  if (parent_name != 0) delete parent_name;
  parent_name = new G4String(particle_name);
  parent = 0;
}

inline
 G4int G4VDecayChannel::GetNumberOfDaughters() const 
{ 
  return  numberOfDaughters;
}

inline
 G4String G4VDecayChannel::GetKinematicsName() const { return kinematics_name; }

inline
 void  G4VDecayChannel::SetBR(G4double value){ rbranch = value; }

inline
 G4double G4VDecayChannel::GetBR() const { return rbranch; }

inline
 void  G4VDecayChannel::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4VDecayChannel::GetVerboseLevel() const { return verboseLevel; }



#endif






