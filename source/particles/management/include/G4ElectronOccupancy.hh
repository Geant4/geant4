// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectronOccupancy.hh,v 1.2 1999-08-30 08:27:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 17 Aug 1999
// ----------------------------------------------------------------
//     This class has information of occupation of electrons 
//     in atomic orbits
//  
//     GetOccupancy(N) gives the number of electron
//     in N-th orbit
//       For example : Carbon atom should be 
//          GetOccupancy(0)      --> 2
//          GetOccupancy(1)      --> 4
//          GetOccupancy(2..7)   --> 0
//
//     GetTotalOccupancy() gives the total number of electrons
// --------------------------------------------------------------- 



#ifndef G4ElectronOccupancy_h
#define G4ElectronOccupancy_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ios.hh"

class G4ElectronOccupancy 
{
 public:
   enum { MaxSizeOfOrbit = 7};

   G4ElectronOccupancy( G4int sizeOrbit = MaxSizeOfOrbit   );
   G4ElectronOccupancy( const G4ElectronOccupancy& right );

   virtual    	       ~G4ElectronOccupancy();

  //  new/delete operators are oberloded to use G4Allocator
     inline void *operator new(size_t);
     inline void operator delete(void *aElectronOccupancy);

 
  //- operators
     G4ElectronOccupancy & operator=(const G4ElectronOccupancy &right);
     G4int operator==(const G4ElectronOccupancy &right) const;
     G4int operator!=(const G4ElectronOccupancy &right) const;
   
   // The following methods returns
   //     0:  if the orbit(atom) is vacant 
   //    >0:  number of electrons in orbit
   G4int  GetTotalOccupancy() const;
   G4int  GetOccupancy(G4int orbit) const;
 
   //
   G4int  AddElectron(G4int orbit, G4int number = 1);
   G4int  RemoveElectron(G4int orbit, G4int number = 1);
   
   G4int  GetSizeOfOrbit() const;
   void   DumpInfo() const;

 private:
   G4int  theSizeOfOrbit;
   G4int  theTotalOccupancy;
   G4int* theOccupancies;

};

extern G4Allocator<G4ElectronOccupancy> aElectronOccupancyAllocator;

// ------------------------
// Inlined operators
// ------------------------

inline void * G4ElectronOccupancy::operator new(size_t)
{
  void * aElectronOccupancy;
  aElectronOccupancy = (void *) aElectronOccupancyAllocator.MallocSingle();
  return aElectronOccupancy;
}

inline void G4ElectronOccupancy::operator delete(void * aElectronOccupancy)
{
  aElectronOccupancyAllocator.FreeSingle((G4ElectronOccupancy *) aElectronOccupancy);
}

inline
 G4int  G4ElectronOccupancy::GetSizeOfOrbit() const
{
  return  theSizeOfOrbit;
}

inline
 G4int G4ElectronOccupancy::GetTotalOccupancy() const
{
  return  theTotalOccupancy;
}

inline
 G4int  G4ElectronOccupancy::GetOccupancy(G4int orbit) const
{
  G4int value = 0;
  if ((orbit >=0)&&(orbit<theSizeOfOrbit)){
    value = theOccupancies[orbit];
  }
  return value;  
}

inline 
 G4int  G4ElectronOccupancy::AddElectron(G4int orbit, G4int number)
{
  G4int value =0;
  if ((orbit >=0)&&(orbit<theSizeOfOrbit)){
    theOccupancies[orbit] += number;
    theTotalOccupancy += number; 
    value = number;   
  }
  return value;
}

inline 
 G4int  G4ElectronOccupancy::RemoveElectron(G4int orbit, G4int number)
{
  G4int value =0;
  if ((orbit >=0)&&(orbit<theSizeOfOrbit) ){
    if ( theOccupancies[orbit] < number ) number = theOccupancies[orbit];
    theOccupancies[orbit] -= number;
    theTotalOccupancy -= number;    
    value = number;
  }
  return value;
}
#endif






