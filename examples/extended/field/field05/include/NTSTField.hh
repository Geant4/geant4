#ifndef NTSTField_hh
#define NTSTField_hh

#include "G4DELPHIMagField.hh"

class NTSTField : public G4DELPHIMagField { 
     public:
	NTSTField(): G4DELPHIMagField(), count(0) {;}
        virtual ~NTSTField() {;}
	
	void GetFieldValue( const G4double yTrack[4] ,
                                  G4double *MagField ) const 
	{
		count++;
		G4DELPHIMagField::GetFieldValue( yTrack, MagField );
	}
	
	G4int GetCount() const { return count; }
	void ClearCount() { count = 0; }
	
     private:
	mutable G4int count;
};

#endif
