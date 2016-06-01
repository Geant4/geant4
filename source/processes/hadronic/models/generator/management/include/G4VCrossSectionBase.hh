#ifndef G4VCrossSectionBase_h
#define G4VCrossSectionBase_h 1

#include "globals.hh"

const G4int maxpsig = 13;
const G4int maxreac = 14;

//****************************************************************************************

class G4VCrossSectionBase 
   {
public:
   //Constructors
   G4VCrossSectionBase();
  ~G4VCrossSectionBase();

public:   
   virtual G4double GetTotCrossSection(G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;
   virtual G4double GetElCrossSection (G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;
   virtual G4double GetInCrossSection (G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;   
   virtual G4double GetAnnihCrossSection(G4int ProjectileEncoding, G4int TargetEncoding, G4double
   Energy) = 0;
   };

//****************************************************************************************
#endif



