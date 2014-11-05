// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFPrimaryGeneratorAction.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Definition of the FFPrimaryGeneratorAction class
//!
//  ================ End Documentation Comments ================
//
//  Modified:
//
//  23-06-14                                              BWendt
//  Added function "GetNeutronSourceCenter()" and supporting class-level fields
//
// -------------------------------------------------------------

#ifndef FFPRIMARYGENERATORACTION
#define FFPRIMARYGENERATORACTION

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include "G4VUserPrimaryGeneratorAction.hh"


class FFPrimaryGeneratorAction
:   public G4VUserPrimaryGeneratorAction
{
public:
// Constructor
    FFPrimaryGeneratorAction();

// Functions
    virtual void GeneratePrimaries(G4Event* event);

// Destructor
    virtual ~FFPrimaryGeneratorAction();
    
private:
// Fields
#ifndef NDEBUG
    G4long fEventNumber;
#endif // NDEBUG
    G4VPhysicalVolume* fH2OPhysical;
    G4VPhysicalVolume* fNeutronPhysical;
    G4Tubs* fNeutronSolid;
    G4ParticleGun* const fParticleGun;
    G4VPhysicalVolume* fTankPhysical;
    
// Functions
    G4ThreeVector GetNeutronSourceCenter(void);
};

#endif //FFPRIMARYGENERATORACTION


