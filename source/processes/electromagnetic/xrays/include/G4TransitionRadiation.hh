// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransitionRadiation.hh,v 1.2 1999-04-13 09:26:51 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4TransitionRadiation  -- header file
//
// Class for description of  transition radiation generated
// by  charged particle crossed interface between material 1
// and material 2 (1 -> 2). Transition radiation could be of kind:
// - optical back
// - optical forward
// - X-ray   forward (for relativistic case Tkin/mass >= 10^2)
//
// GEANT 4 class header file --- Copyright CERN 1995
// CERB Geneva Switzerland
//
// for information related to this code, please, contact
// CERN, CN Division, ASD Group
// History:
// 18.12.97, V. Grichine (Vladimir.Grichine@cern.ch)


#ifndef G4TransitionRadiation_h
#define G4TransitionRadiation_h


#include "G4VDiscreteProcess.hh"
#include "G4Material.hh"
// #include "G4OpBoundaryProcess.hh"

class G4TransitionRadiation : public   G4VDiscreteProcess    
{
public:

// Constructors


  G4TransitionRadiation( const G4String& processName = "TR") ;


//  G4TransitionRadiation(const G4TransitionRadiation& right) ;

// Destructor

   ~G4TransitionRadiation() ;

// Operators
// G4TransitionRadiation& operator=(const G4TransitionRadiation& right) ;
// G4int operator==(const G4TransitionRadiation& right)const ;
// G4int operator!=(const G4TransitionRadiation& right)const ;

// Methods

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType)
        {
          return ( aParticleType.GetPDGCharge() != 0.0 );
        }

	G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double previousStepSize,
				 G4ForceCondition* condition)
        {
          *condition = Forced;
	  return DBL_MAX;      // so TR doesn't limit mean free path
        }

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				       const G4Step&  aStep)
        { 
          ClearNumberOfInteractionLengthLeft();
          return &aParticleChange;
        }




virtual
G4double SpectralAngleTRdensity( G4double energy,
                                 G4double varAngle ) const = 0 ;

G4double IntegralOverEnergy( G4double energy1,
                             G4double energy2,
                             G4double varAngle     ) const ;

G4double IntegralOverAngle( G4double energy,
                            G4double varAngle1,
                            G4double varAngle2     ) const ;

G4double AngleIntegralDistribution( G4double varAngle1,
                                    G4double varAngle2     ) const ;

G4double EnergyIntegralDistribution( G4double energy1,
                                     G4double energy2     )  const   ;



// Access functions


protected :

G4int fMatIndex1 ;                   // index of the 1st material
G4int fMatIndex2 ;                   // index of the 2nd material

private :

G4double fGamma ;

// Local constants
static const G4int fSympsonNumber ; // Accuracy of Sympson integration 10
static const G4int fGammaNumber   ; // = 15
static const G4int fPointNumber   ; // = 100

G4double fMinEnergy ;                //  min TR energy
G4double fMaxEnergy ;                //  max TR energy
G4double fMaxTheta  ;                //  max theta of TR quanta

G4double fSigma1 ;                   // plasma energy Sq of matter1
G4double fSigma2 ;                   // plasma energy Sq of matter2


} ;    

#endif   // G4TransitionRadiation_h
