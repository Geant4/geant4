
// Class describing line current magnetic field.
// fFieldConstant determines the coefficient in the field law.
// The line current is directed along Z axis and crosses the XY plane in the point
// (0,0) .
//
// 3.2.97 V. Grichine

#ifndef G4LINECURRENTMAGFIELD_HH
#define G4LINECURRENTMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4LineCurrentMagField : public G4MagneticField
{
public:
		       
         G4LineCurrentMagField(G4double pFieldConstant) ;	       
           
        ~G4LineCurrentMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                          G4double B[]       ) const ;
			   

protected:

private:
          G4double fFieldConstant ;
              
} ;

#endif
