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
//
// -------------------------------------------------------------------

#include "G4QuadrupoleMagField.hh"
#include "G4RotationMatrix.hh"

static G4RotationMatrix IdentityMatrix; 

G4QuadrupoleMagField::G4QuadrupoleMagField(G4double pGradient)
{
   fGradient = pGradient ;
   fOrigin   = G4ThreeVector( 0.0, 0.0, 0.0) ;
   fpMatrix  = &IdentityMatrix;
}


/////////////////////////////////////////////////////////////////////////

G4QuadrupoleMagField::G4QuadrupoleMagField(G4double pGradient,
                                           G4ThreeVector pOrigin,
                                           G4RotationMatrix* pMatrix)
{
   fGradient    = pGradient ;
   fOrigin      = pOrigin ;
   fpMatrix     = pMatrix ;
}

G4Field* G4QuadrupoleMagField::Clone() const
{
    return new G4QuadrupoleMagField(fGradient, fOrigin, fpMatrix);
}

/////////////////////////////////////////////////////////////////////////

G4QuadrupoleMagField::~G4QuadrupoleMagField()
{
}

void G4QuadrupoleMagField::GetFieldValue( const G4double y[7],
                                                G4double B[3]  ) const
//  with displaced origin and rotation
{
    G4ThreeVector r_global = G4ThreeVector(
                                           y[0] - fOrigin.x(),
                                           y[1] - fOrigin.y(),
                                           y[2] - fOrigin.z());
  
  const G4ThreeVector r_local = (*fpMatrix) * r_global;
  const G4ThreeVector B_local( fGradient * r_local.y(),fGradient * r_local.x(),0);
  const G4ThreeVector B_global = fpMatrix->inverse() * B_local;

   B[0] = B_global.x() ;
   B[1] = B_global.y() ;
   B[2] = B_global.z() ;
}
