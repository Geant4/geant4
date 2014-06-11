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
// -------------------------------------------------------------

#ifndef FFPRIMARYGENERATORACTION
#define FFPRIMARYGENERATORACTION

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Tubs.hh"

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
    G4Tubs* neutronSource;
    G4ParticleGun* const particleGun;
};

#endif //FFPRIMARYGENERATORACTION


