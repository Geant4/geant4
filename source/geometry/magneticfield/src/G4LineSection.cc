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
// $Id: G4LineSection.cc 85845 2014-11-05 15:43:58Z gcosmo $
//
// --------------------------------------------------------------------

#include "G4LineSection.hh" 

G4double G4LineSection::Dist( G4ThreeVector OtherPnt ) const
{
  G4double       dist_sq;  
  G4ThreeVector  VecAZ;
  G4double sq_VecAZ, inner_prod, unit_projection ; 

  VecAZ= OtherPnt - EndpointA;
  sq_VecAZ = VecAZ.mag2();

  inner_prod= VecAtoB.dot( VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 
  //
  if( fABdistanceSq != 0.0 )
  {
    //  unit_projection= inner_prod * InvsqDistAB();
    unit_projection = inner_prod/fABdistanceSq;

    if( (0. <= unit_projection ) && (unit_projection <= 1.0 ) )
    {
      dist_sq= sq_VecAZ -  unit_projection * inner_prod ;
    }
    else
    {
     //  The perpendicular from the point to the line AB meets the line
     //   in a point outside the line segment!
     
      if( unit_projection < 0. ) // A is the closest point
      {
        dist_sq= sq_VecAZ;  
      }
      else                       // B is the closest point
      {
	G4ThreeVector   EndpointB = EndpointA + VecAtoB;
        G4ThreeVector   VecBZ =     OtherPnt - EndpointB;
        dist_sq =  VecBZ.mag2();
      }
    }
  }
  else
  {
     dist_sq = (OtherPnt - EndpointA).mag2() ;   
  }  
  if( dist_sq < 0.0 ) dist_sq = 0.0 ;

  return std::sqrt(dist_sq) ;  
}
