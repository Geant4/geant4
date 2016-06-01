//  Contents ---------------------------------------------------------
//
//	G4Assembly
//
//  Description:
//   
//	C++ header file for ...
//	Uses the xxxxx classes.
//	A G4Assembly is ...
//  End --------------------------------------------------------------

//  Interface Dependencies -------------------------------------------

#ifndef G4ASSEMBLY_HH
#define G4ASSEMBLY_HH

#include "G4PlacedSolid.hh"
#include "G4OrderedTable.hh"   
#include "G4BREPSolid.hh"

typedef RWTPtrOrderedVector<G4PlacedSolid> G4PlacedVector;  
    
//  End Interface Dependencies ---------------------------------------


//  Class  //
class G4Assembly
{

public:
  G4Assembly();
  ~G4Assembly();

  void SetPlacedVector(G4PlacedVector&);

  G4PlacedSolid* GetPlacedSolid(G4int solidNumber)
  {
     return placedVec[solidNumber];
  }

  G4int GetNumberOfSolids()
  {
     return numberOfSolids;
  }


private:  
  G4int           numberOfSolids;
  G4PlacedVector  placedVec;

};

#endif












