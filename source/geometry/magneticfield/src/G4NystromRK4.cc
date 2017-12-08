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
// $Id: G4NystromRK4.cc 107821 2017-12-05 14:14:47Z gunter $
//
// History:
// - Created:      I.Gavrilenko    15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#include <iostream>
#include "G4NystromRK4.hh"

//////////////////////////////////////////////////////////////////
// Constructor - with optional distance ( has default value)
//////////////////////////////////////////////////////////////////

G4NystromRK4::G4NystromRK4(G4Mag_EqRhs* magEqRhs, G4double distanceConstField)
  : G4MagIntegratorStepper(magEqRhs, 6),            // number of variables
    m_fEq( magEqRhs ),
    m_magdistance( distanceConstField ),
    m_cof( 0.0 ),
    m_mom( 0.0 ),
    m_imom( 0.0 ),
    m_cachedMom( false )
{
  m_fldPosition[0]  = m_iPoint[0] = m_fPoint[0] = m_mPoint[0] = 9.9999999e+99 ;
  m_fldPosition[1]  = m_iPoint[1] = m_fPoint[1] = m_mPoint[1] = 9.9999999e+99 ;
  m_fldPosition[2]  = m_iPoint[2] = m_fPoint[2] = m_mPoint[2] = 9.9999999e+99 ;
  m_fldPosition[3]  = -9.9999999e+99;
  m_lastField[0] = m_lastField[1] = m_lastField[2] = 0.0;

  m_magdistance2 = distanceConstField*distanceConstField;
}

////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////

G4NystromRK4::~G4NystromRK4()
{
}

/////////////////////////////////////////////////////////////////////////////////
// Integration in one  step 
/////////////////////////////////////////////////////////////////////////////////

void 
G4NystromRK4::Stepper
(const G4double P[],const G4double dPdS[],G4double Step,G4double Po[],G4double Err[])
{
  const G4double perMillion = 1.0e-6;
  G4double R[4] = {   P[0],   P[1] ,    P[2],  P[7] };   // x, y, z, t
  G4double A[3] = {dPdS[0], dPdS[1], dPdS[2]};

  m_iPoint[0]=R[0]; m_iPoint[1]=R[1]; m_iPoint[2]=R[2];

  constexpr G4double one_sixth= 1./6.;
  const G4double S  =     Step   ;
  const G4double S5 =  .5*Step   ;
  const G4double S4 = .25*Step   ;
  const G4double S6 =     Step * one_sixth;   // Step / 6.;
  
  // Ensure that the location and cached field value are correct
  getField( R );

  // Ensure that the momentum is set correctly.

  // - Quick check momentum magnitude (squared) against previous value
  G4double newmom2 = (P[3]*P[3]+P[4]*P[4]+P[5]*P[5]); 
  G4double oldmom2 = m_mom * m_mom;
  if( std::fabs(newmom2 - oldmom2) > perMillion * oldmom2 ) {
     m_mom   = std::sqrt(newmom2) ;
     m_imom  = 1./m_mom;
     m_cof   = m_fEq->FCof()*m_imom;
  }

#ifdef  G4DEBUG_FIELD
  CheckCachedMomemtum( P, m_mom );
  CheckFieldPosition( P, m_fldPosition );
#endif
  
  // Point 1
  //
  G4double K1[3] = { m_imom*dPdS[3], m_imom*dPdS[4], m_imom*dPdS[5] };
  
  // Point2
  //
  G4double p[4] = {R[0]+S5*(A[0]+S4*K1[0]),
		   R[1]+S5*(A[1]+S4*K1[1]),
		   R[2]+S5*(A[2]+S4*K1[2]),
		   P[7]                   }; 
  getField(p);

  G4double A2[3] = {A[0]+S5*K1[0],A[1]+S5*K1[1],A[2]+S5*K1[2]};
  G4double K2[3] = {(A2[1]*m_lastField[2]-A2[2]*m_lastField[1])*m_cof,
		    (A2[2]*m_lastField[0]-A2[0]*m_lastField[2])*m_cof,
		    (A2[0]*m_lastField[1]-A2[1]*m_lastField[0])*m_cof};
 
  m_mPoint[0]=p[0]; m_mPoint[1]=p[1]; m_mPoint[2]=p[2];

  // Point 3 with the same magnetic field
  //
  G4double A3[3] = {A[0]+S5*K2[0],A[1]+S5*K2[1],A[2]+S5*K2[2]};
  G4double K3[3] = {(A3[1]*m_lastField[2]-A3[2]*m_lastField[1])*m_cof,
		    (A3[2]*m_lastField[0]-A3[0]*m_lastField[2])*m_cof,
		    (A3[0]*m_lastField[1]-A3[1]*m_lastField[0])*m_cof};
  
  // Point 4
  //
  p[0] = R[0]+S*(A[0]+S5*K3[0]);
  p[1] = R[1]+S*(A[1]+S5*K3[1]);
  p[2] = R[2]+S*(A[2]+S5*K3[2]);             

  getField(p);
  
  G4double A4[3] = {A[0]+S*K3[0],A[1]+S*K3[1],A[2]+S*K3[2]};
  G4double K4[3] = {(A4[1]*m_lastField[2]-A4[2]*m_lastField[1])*m_cof,
		    (A4[2]*m_lastField[0]-A4[0]*m_lastField[2])*m_cof,
		    (A4[0]*m_lastField[1]-A4[1]*m_lastField[0])*m_cof};
  
  // New position
  //
  Po[0] = P[0]+S*(A[0]+S6*(K1[0]+K2[0]+K3[0]));
  Po[1] = P[1]+S*(A[1]+S6*(K1[1]+K2[1]+K3[1]));
  Po[2] = P[2]+S*(A[2]+S6*(K1[2]+K2[2]+K3[2]));

  m_fPoint[0]=Po[0]; m_fPoint[1]=Po[1]; m_fPoint[2]=Po[2];

  // New direction
  //
  Po[3] = A[0]+S6*(K1[0]+K4[0]+2.*(K2[0]+K3[0]));
  Po[4] = A[1]+S6*(K1[1]+K4[1]+2.*(K2[1]+K3[1]));
  Po[5] = A[2]+S6*(K1[2]+K4[2]+2.*(K2[2]+K3[2]));

  // Errors
  //
  Err[3] = S*std::fabs(K1[0]-K2[0]-K3[0]+K4[0]);
  Err[4] = S*std::fabs(K1[1]-K2[1]-K3[1]+K4[1]);
  Err[5] = S*std::fabs(K1[2]-K2[2]-K3[2]+K4[2]);
  Err[0] = S*Err[3]                       ;
  Err[1] = S*Err[4]                       ;
  Err[2] = S*Err[5]                       ;
  Err[3]*= m_mom                          ;
  Err[4]*= m_mom                          ;
  Err[5]*= m_mom                          ;

  // Normalize momentum
  //
  G4double normF = m_mom/std::sqrt(Po[3]*Po[3]+Po[4]*Po[4]+Po[5]*Po[5]);
  Po [3]*=normF; Po[4]*=normF; Po[5]*=normF; 

  // Pass Energy, time unchanged -- time is not integrated !!
  Po[6]=P[6]; Po[7]=P[7];
}


/////////////////////////////////////////////////////////////////////////////////
// Estimate the maximum distance from the curve to the chord
/////////////////////////////////////////////////////////////////////////////////

G4double 
G4NystromRK4::DistChord() const 
{
  G4double ax = m_fPoint[0]-m_iPoint[0];  
  G4double ay = m_fPoint[1]-m_iPoint[1];  
  G4double az = m_fPoint[2]-m_iPoint[2];
  G4double dx = m_mPoint[0]-m_iPoint[0]; 
  G4double dy = m_mPoint[1]-m_iPoint[1]; 
  G4double dz = m_mPoint[2]-m_iPoint[2];
  G4double d2 = (ax*ax+ay*ay+az*az)    ; 

  if(d2!=0.) {
    G4double ds = (ax*dx+ay*dy+az*dz)/d2;
    dx         -= (ds*ax)               ;
    dy         -= (ds*ay)               ;
    dz         -= (ds*az)               ;
  }
  return std::sqrt(dx*dx+dy*dy+dz*dz);
}

/////////////////////////////////////////////////////////////////////////////////
// Derivatives calculation - caching the momentum value
/////////////////////////////////////////////////////////////////////////////////

void 
G4NystromRK4::ComputeRightHandSide(const G4double P[],G4double dPdS[])
{
  G4double P4vec[4]= { P[0], P[1], P[2], P[7] }; // Time is P[7]
  getField(P4vec);
  m_mom   = std::sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5])     ; 
  m_imom  = 1./m_mom                                ;
  m_cof   = m_fEq->FCof()*m_imom                    ;
  m_cachedMom = true                                ; // Caching the value
  dPdS[0] = P[3]*m_imom                             ; // dx /ds
  dPdS[1] = P[4]*m_imom                             ; // dy /ds
  dPdS[2] = P[5]*m_imom                             ; // dz /ds
  dPdS[3] = m_cof*(P[4]*m_lastField[2]-P[5]*m_lastField[1]) ; // dPx/ds
  dPdS[4] = m_cof*(P[5]*m_lastField[0]-P[3]*m_lastField[2]) ; // dPy/ds
  dPdS[5] = m_cof*(P[3]*m_lastField[1]-P[4]*m_lastField[0]) ; // dPz/ds
}

////////////////////////////////////////////////////////////////////////////
// Check that the location is (almost) unmoved from 'last' field evaluation
////////////////////////////////////////////////////////////////////////////

G4bool
G4NystromRK4::CheckFieldPosition( const G4double Position[3],
                                  const G4double lastPosition[3] )
{
  G4bool ok= true;
  G4double dx = Position[0] - lastPosition[0];
  G4double dy = Position[1] - lastPosition[1];
  G4double dz = Position[2] - lastPosition[2];
  G4double distMag2 = dx*dx+dy*dy+dz*dz;
  if( distMag2 > m_magdistance2) {
     const G4double allowedDist = std::sqrt( m_magdistance2 );
     G4double dist= std::sqrt( distMag2 );
     G4cerr << " NystromRK4::Stepper> ERROR> Moved from correct field position by "
               << dist <<  "( larger than allowed = " << allowedDist << " ) "
               << G4endl;
     ok= false;
  }
  return ok;
}

////////////////////////////////////////////////////
// Check magnitude of momentum against saved value
////////////////////////////////////////////////////

G4bool G4NystromRK4::CheckCachedMomemtum( const G4double PosMom[6],
                                                G4double savedMom )
{
  constexpr G4double perThousand = 1.0e-3;
  G4bool ok= true;
  G4double new_mom2= (PosMom[3]*PosMom[3]+PosMom[4]*PosMom[4]+PosMom[5]*PosMom[5]);
  G4double new_mom=  std::sqrt(new_mom2); 
  if( std::fabs(new_mom - savedMom ) > perThousand * savedMom ) {
     G4cerr << " Nystrom::Stepper WARNING: momentum magnitude is invalid / has changed "
            << G4endl
            << " new value    (p-mag) = "   << new_mom << G4endl
            << " cached value (p-mag) = "  << savedMom   << G4endl;
     if( savedMom > 0.0 ) {
        G4cerr << " ratio  (new/old) = " << new_mom / savedMom << G4endl;
     }
     ok= false;
  }
  return ok;
}
