#include "G4QuadrupoleMagField.hh"
#include "globals.hh"
#include "geomdefs.hh"



G4QuadrupoleMagField::G4QuadrupoleMagField(G4double pGradient)
{
   fGradient = pGradient ;
}

/////////////////////////////////////////////////////////////////////////

G4QuadrupoleMagField::~G4QuadrupoleMagField()
{
   ;
}

////////////////////////////////////////////////////////////////////////


void
   G4QuadrupoleMagField::GetFieldValue( const G4double y[7],
                                              G4double B[3]  ) const  
{
   //   G4double fGradient = 0.001 ;   // Tesla/mm
   B[0] = fGradient*y[1] ;
   B[1] = fGradient*y[0] ;
   B[2] = 0 ;
   return ;
}



