#ifndef G4tgrPlaceDivRep_h
#define G4tgrPlaceDivRep_h
/*---------------------------------------------------------------------------   
ClassName:   G4tgrPlaceDivRep
Author:      P. Arce
Changes:     12/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to descripe the position of a G4tgrVolume inside another G4tgrVolume as a replica
(filling the whole parent volume)*/ 
//----------------------------------------------------------------------------  
#include "globals.hh"
#include "geomdefs.hh"
#include "G4tgrPlace.hh"
#include "CLHEP/Vector/ThreeVector.h"

enum DivType {DivByNdiv,DivByWidth,DivByNdivAndWidth};

class G4tgrPlaceDivRep : public G4tgrPlace
{

 public:
  G4tgrPlaceDivRep();
  ~G4tgrPlaceDivRep(){ };

  // creates an object passing the only data that is fixed (ndiv, width, offset may be have to be recalculated)
  G4tgrPlaceDivRep( const std::vector<G4String>& wl );

  EAxis BuildAxis( const G4String& axisName );

  //! access functions
  EAxis GetAxis() const {return theAxis;}
  int GetNDiv() const {return theNDiv;}
  double GetWidth() const {return theWidth;}
  double GetOffset() const {return theOffset;}
  DivType GetDivType() const {return theDivType;}

  void SetParentName( G4String& parentName ) {theParentName = parentName;}
  void SetNDiv( int ndiv ) { theNDiv = ndiv; }
  void SetWidth( double width ) { theWidth = width; }
  void SetAxis( EAxis axis ) { theAxis = axis; }
  void SetOffset( double offset ) { theOffset = offset; }
  void SetDivType( DivType typ ) { theDivType = typ; }

 private:
  int theNDiv;
  double theWidth;
  EAxis theAxis;
  double theOffset;
  DivType theDivType;
  
};
#endif
