// Class describing magnetic field parametrised by harmonic polynom up to
// 3rd order. The function MagneticField(yTrack,B) calculates the magnetic field
// induction vector B for the trajectory point yTrack according to formula given
// in: M.Metcalf, Analysis of the SFM Field, OM Development Note AP-10 (revised)
// 09.01.1974
//
// 3.2.97 V.Grichine

#ifndef G4HARMONICPOLMAGFIELD_HH
#define G4HARMONICPOLMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4HarmonicPolMagField : public G4MagneticField
{
public:
		       
         G4HarmonicPolMagField() ;	       
           
        ~G4HarmonicPolMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                    G4double B[]              ) const  ;
			   

protected:

private:
              
} ;

#endif
