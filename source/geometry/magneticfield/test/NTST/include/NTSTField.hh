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
