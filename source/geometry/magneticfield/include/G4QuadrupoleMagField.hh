// Class for creation of quadrupole magnetic field
// fGradient - is gardient value for  quadrupole magnetic lense. Then the magnetic
// field components are:
// Bx = B[0] = fGradient*X ,
// By = B[1] = fGradient*Y ,
// Bz = B[2] = 0 .
// Here X,Y,Z are the coordinates of a space point of interest.
//
// 3.2.97 V.Grichine

#ifndef G4QUADRUPOLEMAGFIELD_HH
#define G4QUADRUPOLEMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4QuadrupoleMagField : public G4MagneticField
{
public:
		       
	 G4QuadrupoleMagField(G4double pGradient) ;	       
           
        ~G4QuadrupoleMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                          G4double B[]        ) const ;
			   

protected:

private:
          G4double fGradient ;
              
} ;

#endif
