#ifndef G4QGSModel_h
#define G4QGSModel_h 1

// Class Description
// Model for hadron (p,n,pi,K) nuclear reactions in geant4. IT implements
// A. Kaydalov's quark gluon string model.
// To be used in your physics list, in case you need this kind of physics.
// Class Description - End

#include "G4ExcitedStringVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4PomeronCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4VPartonStringModel.hh"
#include "G4QGSParticipants.hh"
#include "G4DiffractiveStringBuilder.hh"
#include "G4SoftStringBuilder.hh"
#include "G4PartonPair.hh"

//*********************************************************************************************** 


//*****************************************************************************************

class G4QGSModel : public G4VPartonStringModel
    {
// Constructors   
public:
    G4QGSModel();
    G4QGSModel(const G4QGSModel &right);
    virtual ~G4QGSModel();

// Method
public:
    virtual G4V3DNucleus* GetWoundedNucleus() const;
 
public:
    virtual void Init(const G4Nucleus& Nucleus, const G4DynamicParticle& Projectile);
    virtual G4ExcitedStringVector * GetStrings();
 
private:
   G4QGSParticipants theParticipants;
   G4DiffractiveStringBuilder theDiffractiveStringBuilder;
   G4SoftStringBuilder theSoftStringBuilder;

private:
   // cash theCurrentVelocity for lorentztrafo HPW 
   G4ThreeVector theCurrentVelocity;

   };

//-------------------------------------------------------------------------------------------


//*****************************************************************************************
    
#endif


