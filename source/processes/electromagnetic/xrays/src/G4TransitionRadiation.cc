// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TransitionRadiation.cc,v 1.1 1999-01-07 16:11:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4TransitionRadiation class -- implementation file

// GEANT 4 class implementation file --- Copyright CERN 1995
// CERN Geneva Switzerland

// For information related to this code, please, contact
// CERN, CN Division, ASD Group
// History:
// 1st version 11.09.97 V. Grichine (Vladimir.Grichine@cern.ch )
// 2nd version 16.12.97 V. Grichine


#include <math.h>
// #include "G4ios.hh"
// #include <fstream.h>
// #include <stdlib.h>

#include "G4TransitionRadiation.hh"
#include "G4Material.hh"

// Init gamma array


// Local constants

const G4int   G4TransitionRadiation::fSympsonNumber = 100 ;
const G4int   G4TransitionRadiation::fGammaNumber = 15 ;
const G4int   G4TransitionRadiation::fPointNumber = 100 ;


///////////////////////////////////////////////////////////////////////
//
// Constructor for selected couple of materials
//

G4TransitionRadiation::
G4TransitionRadiation( const G4String& processName )
  : G4VDiscreteProcess(processName)
{
  //  fMatIndex1 = pMat1->GetIndex() ;
  //  fMatIndex2 = pMat2->GetIndex() ;
}

//////////////////////////////////////////////////////////////////////
//
// Destructor
//

G4TransitionRadiation::~G4TransitionRadiation()
{
	;
}


///////////////////////////////////////////////////////////////////
//
// Sympson integral of TR spectral-angle density over energy between
// the limits energy 1 and energy2 at fixed varAngle = 1 - cos(Theta)

G4double
G4TransitionRadiation::IntegralOverEnergy( G4double energy1,
                                           G4double energy2,
                                           G4double varAngle     )  const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(energy2 - energy1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
    sumEven += SpectralAngleTRdensity(energy1 + 2*i*h,varAngle)  ;
    sumOdd  += SpectralAngleTRdensity(energy1 + (2*i - 1)*h,varAngle) ;
  }
  sumOdd += SpectralAngleTRdensity(energy1 + (2*fSympsonNumber - 1)*h,varAngle) ;
  return h*(  SpectralAngleTRdensity(energy1,varAngle)
            + SpectralAngleTRdensity(energy2,varAngle)
            + 4.0*sumOdd + 2.0*sumEven    )/3.0 ;
}



///////////////////////////////////////////////////////////////////
//
// Sympson integral of TR spectral-angle density over energy between
// the limits varAngle1 and varAngle2 at fixed energy

G4double
G4TransitionRadiation::IntegralOverAngle( G4double energy,
                                          G4double varAngle1,
                                          G4double varAngle2     ) const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(varAngle2 - varAngle1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
    sumEven += SpectralAngleTRdensity(energy,varAngle1 + 2*i*h)  ;
    sumOdd  += SpectralAngleTRdensity(energy,varAngle1 + (2*i - 1)*h) ;
  }
  sumOdd += SpectralAngleTRdensity(energy,varAngle1 + (2*fSympsonNumber - 1)*h) ;

  return h*(  SpectralAngleTRdensity(energy,varAngle1)
            + SpectralAngleTRdensity(energy,varAngle2)
            + 4.0*sumOdd + 2.0*sumEven    )/3.0 ;
}

///////////////////////////////////////////////////////////////////
//
// The number of transition radiation photons generated in the
// angle interval between varAngle1 and varAngle2
//

G4double G4TransitionRadiation::
AngleIntegralDistribution( G4double varAngle1,
                           G4double varAngle2     )   const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(varAngle2 - varAngle1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
   sumEven += IntegralOverEnergy(fMinEnergy,
                                 fMinEnergy +0.3*(fMaxEnergy-fMinEnergy),
                                 varAngle1 + 2*i*h)
            + IntegralOverEnergy(fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                                 fMaxEnergy,
                                 varAngle1 + 2*i*h);
   sumOdd  += IntegralOverEnergy(fMinEnergy,
                                 fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                                 varAngle1 + (2*i - 1)*h)
            + IntegralOverEnergy(fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                                 fMaxEnergy,
                                 varAngle1 + (2*i - 1)*h) ;
  }
  sumOdd += IntegralOverEnergy(fMinEnergy,
                               fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               varAngle1 + (2*fSympsonNumber - 1)*h)
          + IntegralOverEnergy(fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               fMaxEnergy,
                               varAngle1 + (2*fSympsonNumber - 1)*h) ;

  return h*(IntegralOverEnergy(fMinEnergy,
                               fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               varAngle1)
          + IntegralOverEnergy(fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               fMaxEnergy,
                               varAngle1)
          + IntegralOverEnergy(fMinEnergy,
                               fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               varAngle2)
          + IntegralOverEnergy(fMinEnergy + 0.3*(fMaxEnergy - fMinEnergy),
                               fMaxEnergy,
                               varAngle2)
            + 4.0*sumOdd + 2.0*sumEven    )/3.0 ;
}

///////////////////////////////////////////////////////////////////
//
// The number of transition radiation photons, generated in the
// energy interval between energy1 and energy2
//

G4double G4TransitionRadiation::
EnergyIntegralDistribution( G4double energy1,
                            G4double energy2     )  const
{
  G4int i ;
  G4double h , sumEven = 0.0 , sumOdd = 0.0 ;
  h = 0.5*(energy2 - energy1)/fSympsonNumber ;
  for(i=1;i<fSympsonNumber;i++)
  {
   sumEven += IntegralOverAngle(energy1 + 2*i*h,0.0,0.01*fMaxTheta )
            + IntegralOverAngle(energy1 + 2*i*h,0.01*fMaxTheta,fMaxTheta);
   sumOdd  += IntegralOverAngle(energy1 + (2*i - 1)*h,0.0,0.01*fMaxTheta)
            + IntegralOverAngle(energy1 + (2*i - 1)*h,0.01*fMaxTheta,fMaxTheta) ;
  }
  sumOdd += IntegralOverAngle(energy1 + (2*fSympsonNumber - 1)*h,
                              0.0,0.01*fMaxTheta)
          + IntegralOverAngle(energy1 + (2*fSympsonNumber - 1)*h,
                              0.01*fMaxTheta,fMaxTheta) ;

  return h*(IntegralOverAngle(energy1,0.0,0.01*fMaxTheta)
          + IntegralOverAngle(energy1,0.01*fMaxTheta,fMaxTheta)
          + IntegralOverAngle(energy2,0.0,0.01*fMaxTheta)
          + IntegralOverAngle(energy2,0.01*fMaxTheta,fMaxTheta)
            + 4.0*sumOdd + 2.0*sumEven    )/3.0 ;
}




// end of G4TransitionRadiation implementation file --------------------------
