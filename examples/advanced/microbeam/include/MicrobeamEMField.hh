// -------------------------------------------------------------------
// $Id: MicrobeamEMField.hh,v 1.3 2006-06-01 22:25:19 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamEMField_h
#define MicrobeamEMField_h 1

#include "globals.hh"
#include "G4ElectroMagneticField.hh"

class MicrobeamEMField
#ifndef STANDALONE
 : public G4ElectroMagneticField
#endif

{
  
public:
  MicrobeamEMField();
  void  GetFieldValue( const  double Point[4], double *Bfield ) const;
		       
  G4bool DoesFieldChangeEnergy() const {return true;}

};

#endif
