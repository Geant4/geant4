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
// $Id: G4FieldTrack.cc,v 1.15 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "G4FieldTrack.hh"

std::ostream& operator<<( std::ostream& os, const G4FieldTrack& SixVec)
{
     const G4double *SixV = SixVec.SixVector;
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " "
                  << SixV[2] << " ";  // Position
     os << " V= " << SixV[3] << " " << SixV[4] << " "
                  << SixV[5] << " ";  // Momentum
     os << " v2= "
        << G4ThreeVector(SixV[3], SixV[4], SixV[5]).mag(); // mom magnitude
     os << " mdm= " << SixVec.fMomentumDir.mag(); 
     os << " l= " << SixVec.GetCurveLength();
     os << " ) ";
     return os;
}

G4FieldTrack::G4FieldTrack( const G4ThreeVector& pPosition, 
			          G4double       LaboratoryTimeOfFlight,
			    const G4ThreeVector& pMomentumDirection,
			          G4double       kineticEnergy,
			          G4double       restMass_c2,
		                  G4double       charge, 
			    const G4ThreeVector& Spin,
			          G4double       magnetic_dipole_moment,
			          G4double       curve_length )
 : fKineticEnergy(kineticEnergy),
   fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(LaboratoryTimeOfFlight), 
   fProperTimeOfFlight(0.),
   // fMomentumDir(pMomentumDirection),
   fChargeState(  charge, magnetic_dipole_moment ) 
{
  G4double momentum  = std::sqrt(kineticEnergy*(kineticEnergy+2.0*restMass_c2));
  G4ThreeVector pMomentum= momentum * pMomentumDirection; 
  SetCurvePnt( pPosition, pMomentum, curve_length );
    // Sets momentum direction as well.

  fMomentumDir=pMomentumDirection; 
    // Set the momentum direction again - keeping value from argument exactly

  InitialiseSpin( Spin ); 
}

G4FieldTrack::G4FieldTrack( const G4ThreeVector& pPosition, 
                            const G4ThreeVector& pMomentumDirection,    
                                  G4double       curve_length, 
                                  G4double       kineticEnergy,
                            const G4double       restMass_c2,
                                  G4double,   // velocity
                                  G4double       pLaboratoryTimeOfFlight,
                                  G4double       pProperTimeOfFlight,
                            const G4ThreeVector* pSpin)
 : fKineticEnergy(kineticEnergy),
   fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(pLaboratoryTimeOfFlight), 
   fProperTimeOfFlight(pProperTimeOfFlight),
   // fMomentumDir(pMomentumDirection), 
   fChargeState( DBL_MAX ) //  charge not set 
{
  G4double momentum  = std::sqrt(kineticEnergy*(kineticEnergy+2.0*restMass_c2));
  G4ThreeVector pMomentum= momentum * pMomentumDirection; 

  SetCurvePnt( pPosition, pMomentum, curve_length );
  // Sets momentum direction as well.

  // Set the momentum direction again
  //   -- to avoid numerical issues from multiplying by momentum and dividing again
  fMomentumDir=pMomentumDirection; 

  G4ThreeVector Spin(0.0, 0.0, 0.0); 
  if( !pSpin ) Spin= G4ThreeVector(0.,0.,0.); 
  else         Spin= *pSpin;
  InitialiseSpin( Spin ); 
}

G4FieldTrack::G4FieldTrack( char )                  //  Nothing is set !!
  : fKineticEnergy(0.), fRestMass_c2(0.), fLabTimeOfFlight(0.),
    fProperTimeOfFlight(0.), fChargeState( DBL_MAX )
{
  G4ThreeVector Zero(0.0, 0.0, 0.0);
  SetCurvePnt( Zero, Zero, 0.0 );
  InitialiseSpin( Zero ); 
}

void G4FieldTrack::
     SetChargeAndMoments(G4double charge, 
			 G4double magnetic_dipole_moment, // default= DBL_MAX - do not change
			 G4double electric_dipole_moment, //   ditto
			 G4double magnetic_charge )       //   ditto
{
  fChargeState.SetChargeAndMoments( charge,  magnetic_dipole_moment, 
		      electric_dipole_moment,  magnetic_charge ); 

  // fpChargeState->SetChargeAndMoments( charge,  magnetic_dipole_moment, 
  //	      electric_dipole_moment,  magnetic_charge ); 

  // TO-DO: Improve the implementation using handles
  //   -- and handle to the old one (which can be shared by other copies) and
  //      must not be left to hang loose 
  // 
  // fpChargeState= new G4ChargeState(  charge, magnetic_dipole_moment, 
  //			     electric_dipole_moment, magnetic_charge  ); 
}
