//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
