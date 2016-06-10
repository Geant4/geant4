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
// $Id: G4HarmonicPolMagField.cc 68055 2013-03-13 14:43:28Z gcosmo $
//
// -------------------------------------------------------------------

#include "G4HarmonicPolMagField.hh"
#include "globals.hh"

G4HarmonicPolMagField::G4HarmonicPolMagField()
{
}

G4HarmonicPolMagField* G4HarmonicPolMagField::Clone() const
{
    return new G4HarmonicPolMagField;
}
/////////////////////////////////////////////////////////////////////////  

G4HarmonicPolMagField::~G4HarmonicPolMagField()
{
}

/////////////////////////////////////////////////////////////////////////


void G4HarmonicPolMagField::GetFieldValue(const G4double yTrack[7],
                                                G4double B[3]     ) const  
{
   G4int i ;
   G4double a = 1.00 ;   // mm -> m 
   G4double x = a*yTrack[0], y = a*yTrack[1], z = a*yTrack[2] ;
   G4double x2 = x*x, y2 = y*y, z2 = z*z ;
   G4double x3 = x2*x, y3 = y2*y, z3 = z2*z ;
   G4double xy = x*y, xz = x*z, yz = y*z, xyz = x*y*z ;
   static G4ThreadLocal G4double
   c[24] = {
             .010, .010, .010,                                       // 3(0)
	     .0001, .0001, .0001, .0001, .0001,                      // 5(1)
	     .00001, .00001, .00001, .00001, .00001, .00001, .00001, // 7(2)
	     .000001, .000001, .000001, .000001, .000001, .000001,
	     .0000001, .0000001, .0000001                            // 9(3)
            } ;                                            // total :   24
   
   //   for(i=0;i<24;i++)
   //   {
   //      c[i] = 1.0*c[i] ;
   //   }
   B[0] =  c[1]
          -2*c[3]*x + c[4]*z +c[6]*y -2*c[7]*x
          -6*c[8]*xz + c[9]*(z2-x2) -2*c[10]*xy + c[11]*yz - 2*c[12]*xz
          +c[13]*(y2-x2) - 6*c[14]*xy
          -4*c[15]*(3*x*z2-x3) +c[16]*(z3-3*x2*z) - 6*c[17]*xyz +c[18]*y*(z2-x2)
          -2*c[19]*(x*z2+x*y2-2*x3/3) + c[20]*z*(y2-x2) - 6*c[21]*xyz
          +c[22]*(y3-3*x2*y) - 4*c[23]*(3*x*y2-x3) ;
   
   B[1] =  c[2]
          +c[5]*z + c[6]*x + 2*c[7]*y
          +c[10]*(z2-x2) + c[11]*xz +2*c[12]*yz +2*c[13]*xy + 3*c[14]*(y2-x2)
	  +c[17]*(z3-3*x2*z) + c[18]*(x*z2-x3/3) +2*c[19]*y*(z2-x2)
          +2*c[20]*xyz
          +3*c[21]*z*(y2-x2) + c[22]*(3*x*y2-x3) + 4*c[23]*(y3-3*x2*y) ;
   
   B[2] =  c[0]
          +c[3]*z + c[4]*x + c[5]*y
          +3*c[8]*(z2-x2) + 2*c[9]*xz + 2*c[10]*yz + c[11]*xy + c[12]*(y2-x2)
          +4*c[15]*(z3-3*x2*z) + c[16]*(3*x*z2-x3) + 3*c[17]*(y*z2-x2*y)
          +2*c[18]*xyz
          +2*c[19]*z*(y2-x2) + c[20]*(x*y2-x3/3) + c[21]*(y3-3*x2*y) ;
   for(i=0;i<3;i++)
   {
      B[i] = 0.1*B[i] ;
   }
}

// -----------------------------------------------------------------
