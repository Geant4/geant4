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
#ifndef NTSTField_hh
#define NTSTField_hh

#include "G4UniformMagField.hh"

class NTSTField : public G4UniformMagField {
	public:
	NTSTField( const G4ThreeVector& FieldVector )
		: G4UniformMagField(FieldVector),
		  count(0) {;}
	NTSTField( G4double vField,
                   G4double vTheta,
                   G4double vPhi     )
		: G4UniformMagField( vField, vTheta, vPhi ),
		  count(0) {;}
	virtual ~NTSTField() {;}
	
	
	void GetFieldValue( const G4double yTrack[3] ,
                                  G4double *MagField ) const 
	{
		count++;
		G4UniformMagField::GetFieldValue( yTrack, MagField );
	}
	
	G4int GetCount() const { return count; }
	void ClearCount() { count = 0; }
	
	protected:
	mutable G4int count;
};

#endif
