#include "G4LineCurrentMagField.hh"
#include "globals.hh"
#include "geomdefs.hh"



G4LineCurrentMagField::G4LineCurrentMagField(G4double pFieldConstant)
{
   fFieldConstant = pFieldConstant ;
}
////////////////////////////////////////////////////////////////////////

G4LineCurrentMagField::~G4LineCurrentMagField()
{
   ;
}

///////////////////////////////////////////////////////////////////////////


void
   G4LineCurrentMagField::GetFieldValue( const G4double yTrack[7],
                                               G4double B[3]      ) const  
{
   //   G4double fFieldConstant = 100 ;
   G4double a = 1.00*mm ;   // mm -> m 
   G4double x = a*yTrack[0], y = a*yTrack[1], z = a*yTrack[2] ;
   G4double x2 = x*x, y2 = y*y, r2 = x2 + y2 ;
   G4double r = sqrt(r2+a*a) ;
   G4double Br = fFieldConstant/r;
   B[0] = -Br*y/r ;
   B[1] = Br*x/r ;
   B[2] = 0 ;
   return ;
}

// -----------------------------------------------------------------

