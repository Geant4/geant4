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
// G4FieldTrack implementation
//
// Author: John Apostolakis, CERN - First version, 14.10.1996
// -------------------------------------------------------------------

#include "G4FieldTrack.hh"

std::ostream& operator<<( std::ostream& os, const G4FieldTrack& SixVec)
{
     const G4double* SixV = SixVec.SixVector;
     const G4int precPos= 9;   // For position
     const G4int precEp=  9;   // For Energy / momentum
     const G4int precLen= 12;  // For Length along track
     const G4int precSpin= 9;  // For polarisation
     const G4int precTime= 6;  // For time of flight
     const G4long oldpr= os.precision(precPos);
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " "
                  << SixV[2] << " ";  // Position
     os.precision(precEp);     
     os << " P= " << SixV[3] << " " << SixV[4] << " "
                  << SixV[5] << " ";  // Momentum
     os << " Pmag= "
        << G4ThreeVector(SixV[3], SixV[4], SixV[5]).mag(); // mom magnitude
     os << " Ekin= " << SixVec.fKineticEnergy ;
     os.precision(precLen);
     os << " l= " << SixVec.GetCurveLength();
     os.precision(6);
     os << " m0= " <<   SixVec.fRestMass_c2;
     os << " (Pdir-1)= " <<  SixVec.fMomentumDir.mag()-1.0;
     if( SixVec.fLabTimeOfFlight > 0.0 )
     {
       os.precision(precTime);
     }
     else
     {
       os.precision(3);
     }
     os << " t_lab= "    << SixVec.fLabTimeOfFlight; 
     os << " t_proper= " << SixVec.fProperTimeOfFlight ;
     G4ThreeVector pol= SixVec.GetPolarization();
     if( pol.mag2() > 0.0 )
     {
        os.precision(precSpin);
        os << " PolV= " << pol; // SixVec.GetPolarization();
     }
     else
     {
        os << " PolV= (0,0,0) "; 
     }
     os << " ) ";
     os.precision(oldpr);
     return os;
}

G4FieldTrack::G4FieldTrack( const G4ThreeVector& pPosition, 
			          G4double       LaboratoryTimeOfFlight,
			    const G4ThreeVector& pMomentumDirection,
			          G4double       kineticEnergy,
			          G4double       restMass_c2,
		                  G4double       charge, 
			    const G4ThreeVector& vecPolarization,
			          G4double       magnetic_dipole_moment,
                                  G4double       curve_length,
                                  G4double       pdgSpin )
:  fDistanceAlongCurve(curve_length),
   fKineticEnergy(kineticEnergy),
   fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(LaboratoryTimeOfFlight), 
   fProperTimeOfFlight(0.),
   // fMomentumDir(pMomentumDirection),
   fChargeState(  charge, magnetic_dipole_moment, pdgSpin ) 
   // fChargeState(  charge, magnetic_dipole_moment ) , 
   // fPDGSpin( pdgSpin )
{
  UpdateFourMomentum( kineticEnergy, pMomentumDirection ); 
    // Sets momentum direction as well.

  SetPosition( pPosition );
  SetPolarization( vecPolarization ); 
}

G4FieldTrack::G4FieldTrack( const G4ThreeVector& pPosition, 
                            const G4ThreeVector& pMomentumDirection,    
                                  G4double       curve_length, 
                                  G4double       kineticEnergy,
                            const G4double       restMass_c2,
                                  G4double,   // velocity
                                  G4double       pLaboratoryTimeOfFlight,
                                  G4double       pProperTimeOfFlight,
                            const G4ThreeVector* pPolarization,
                                  G4double       pdgSpin )
 : fDistanceAlongCurve(curve_length),
   fKineticEnergy(kineticEnergy),
   fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(pLaboratoryTimeOfFlight), 
   fProperTimeOfFlight(pProperTimeOfFlight),
   fChargeState( DBL_MAX, DBL_MAX, -1.0 ) //  charge not set 
{
  UpdateFourMomentum( kineticEnergy, pMomentumDirection ); 
    // Sets momentum direction as well.
    
  SetPosition( pPosition );    
  fChargeState.SetPDGSpin( pdgSpin );   

  G4ThreeVector PolarVec(0.0, 0.0, 0.0); 
  if( pPolarization )  { PolarVec= *pPolarization; }
  SetPolarization( PolarVec );
}

G4FieldTrack::G4FieldTrack( char )                  //  Nothing is set !!
  : fKineticEnergy(0.), fRestMass_c2(0.), fLabTimeOfFlight(0.),
    fProperTimeOfFlight(0.), fChargeState( DBL_MAX , DBL_MAX, -1 )
{
  G4ThreeVector Zero(0.0, 0.0, 0.0);
  SetCurvePnt( Zero, Zero, 0.0 );
  SetPolarization( Zero ); 
  // fInitialMomentumMag = 0.00; // Invalid
  // fLastMomentumMag = 0.0; 
}

void G4FieldTrack::
     SetChargeAndMoments(G4double charge, 
			 G4double magnetic_dipole_moment, // default = DBL_MAX
			 G4double electric_dipole_moment, // ditto
			 G4double magnetic_charge )       // ditto
{
  fChargeState.SetChargesAndMoments( charge,  
                                     magnetic_dipole_moment, 
                                     electric_dipole_moment,  
                                     magnetic_charge ); 

  // NOTE: Leaves Spin unchanged !
  // 
  // G4double pdgSpin= fChargeState.GetSpin();
  // New Property of ChargeState (not well documented! )

  // IDEA: Improve the implementation using handles
  //   -- and handle to the old one (which can be shared by other copies) and
  //      must not be left to hang loose 
  // 
  // fpChargeState= new G4ChargeState(  charge, magnetic_dipole_moment, 
  //			     electric_dipole_moment, magnetic_charge  ); 
}

// Load values from array
//  
// Note that momentum direction must-be/is normalised
//
void G4FieldTrack::LoadFromArray(const G4double valArrIn[ncompSVEC],
                                       G4int noVarsIntegrated)
{
  // Fill the variables not integrated with zero -- so it's clear !!
  //
  G4double valArr[ncompSVEC];
  for(G4int i=0; i<noVarsIntegrated; ++i)
  {
     valArr[i] = valArrIn[i];
  }
  for(G4int i=noVarsIntegrated; i<ncompSVEC; ++i)
  {
     valArr[i] = 0.0; 
  }

  SixVector[0] = valArr[0];
  SixVector[1] = valArr[1];
  SixVector[2] = valArr[2];
  SixVector[3] = valArr[3];
  SixVector[4] = valArr[4];
  SixVector[5] = valArr[5];

  G4ThreeVector Momentum(valArr[3],valArr[4],valArr[5]);

  G4double momentum_square= Momentum.mag2();
  fMomentumDir= Momentum.unit();

  fKineticEnergy = momentum_square
                 / (std::sqrt(momentum_square+fRestMass_c2*fRestMass_c2)
                   + fRestMass_c2 ); 
    // The above equation is stable for small and large momenta

  // The following components may or may not be
  // integrated over -- integration is optional
  // fKineticEnergy = valArr[6];

  fLabTimeOfFlight = valArr[7];
  fProperTimeOfFlight = valArr[8];
  G4ThreeVector vecPolarization= G4ThreeVector(valArr[9],valArr[10],valArr[11]);
  SetPolarization( vecPolarization ); 

  // fMomentumDir=G4ThreeVector(valArr[13],valArr[14],valArr[15]);
  // fDistanceAlongCurve= valArr[]; 
}
