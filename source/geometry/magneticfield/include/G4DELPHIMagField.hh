// Class describing the DELPHI magnetic field. This axial symmetry
// field mainly directed along Z axis. The function MagneticField(yTrack,B)
// calculates the magnetic induction vector B in point corresponding to
// yTrack according to parametrization given in:
// P.Billoir, Precise tracking in a quasi-honogeneous magnetic field, DELPHI 
// 87-6 PROG 65, 3 February 1987
//
// 3.2.97 V. Grichine 

#ifndef G4DELPHIMAGFIELD_HH
#define G4DELPHIMAGFIELD_HH

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4DELPHIMagField : public G4MagneticField
{
public:
		       
         G4DELPHIMagField() ;	       
           
        ~G4DELPHIMagField() ;
	     
         void GetFieldValue(const G4double yTrack[] ,
	                    G4double B[]              ) const ;
			   

protected:

private:
              
} ;

#endif
