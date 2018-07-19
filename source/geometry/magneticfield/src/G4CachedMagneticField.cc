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
// $Id: G4CachedMagneticField.cc 96751 2016-05-04 09:39:38Z gcosmo $
//
// --------------------------------------------------------------------

#include "G4CachedMagneticField.hh"

G4CachedMagneticField::G4CachedMagneticField(G4MagneticField *pMagField, 
					     G4double         distance)
  : G4MagneticField(), 
    fLastLocation(DBL_MAX,DBL_MAX,DBL_MAX),
    fLastValue(DBL_MAX,DBL_MAX,DBL_MAX),
    fCountCalls(0),  fCountEvaluations(0)
{
  fpMagneticField= pMagField;
  fDistanceConst= distance;

  // G4cout << " Cached-B-Field constructor> Distance = " << distance << G4endl;
  this->ClearCounts(); 
}

G4Field* G4CachedMagneticField::Clone() const
{
  // Cannot use copy constructor: I need to clone the associated magnetic field
  G4MagneticField* aF = static_cast<G4MagneticField*>(fpMagneticField->Clone());
  G4CachedMagneticField* cloned = new G4CachedMagneticField(aF, fDistanceConst);

  cloned->fLastLocation = fLastLocation;
  cloned->fLastValue = fLastValue;
  return cloned;
}

G4CachedMagneticField::~G4CachedMagneticField()
{
}

void
G4CachedMagneticField::ReportStatistics()
{
  G4cout << " Cached field: " << G4endl
	 << "   Number of calls:        " << fCountCalls << G4endl
	 << "   Number of evaluations : " << fCountEvaluations << G4endl;                     
}

G4CachedMagneticField::G4CachedMagneticField(const G4CachedMagneticField &rightCMF)
  : G4MagneticField(rightCMF)
{
  fpMagneticField= rightCMF.fpMagneticField;  // NOTE: sharing pointer here!
  fDistanceConst = rightCMF.fDistanceConst;
  fLastLocation  = rightCMF.fLastLocation;
  fLastValue     = rightCMF.fLastValue;
  this->ClearCounts(); 
}

G4CachedMagneticField& G4CachedMagneticField::operator = (const G4CachedMagneticField &p)
{
  if (&p == this) return *this;
  G4MagneticField::operator=(p);
  fpMagneticField= p.fpMagneticField;  // NOTE: sharing pointer here!
  fDistanceConst = p.fDistanceConst;
  fLastLocation  = p.fLastLocation;
  fLastValue     = p.fLastValue;
  this->ClearCounts(); 
  return *this;
}

void
G4CachedMagneticField::GetFieldValue( const G4double Point[4],
                                            G4double *Bfield ) const
{
  G4ThreeVector newLocation( Point[0], Point[1], Point[2] );

  // G4cout << "Cache-B-field called at " << newLocation << G4endl;

  G4double      distSq= (newLocation-fLastLocation).mag2();
  fCountCalls++;
  if( distSq < fDistanceConst*fDistanceConst ) { 
     Bfield[0] = fLastValue.x();
     Bfield[1] = fLastValue.y();
     Bfield[2] = fLastValue.z();
  }else{
     // G4CachedMagneticField* thisNonC= const_cast<G4CachedMagneticField*>(this);
     fpMagneticField->GetFieldValue( Point, Bfield );
     // G4cout << " Evaluating. " << G4endl;
     fCountEvaluations++;
     // thisNonC->
     fLastLocation= G4ThreeVector( Point[0],  Point[1],  Point[2] );
     // thisNonC->
     fLastValue=    G4ThreeVector( Bfield[0], Bfield[1], Bfield[2] );
  }
}
