// Rich advanced example for Geant4
// AerogelRefData.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef AerogelRefData_h
#define AerogelRefData_h 1
#include "G4ios.hh" 
#include "globals.hh"
#include <vector>
#include "AerogelTypeSpec.hh"

class AerogelRefData {

public: 
 
  AerogelRefData(G4String AerogelRefInputFileName);
  virtual ~AerogelRefData();

  void ReadStdAerogelRefIndex();
  vector<G4double> GetCurAerogelRefIndValueVect(G4int AerogelTypenum );
  G4double GetCurAerogelRefIndValue( G4int rbinw , G4int AerogelTypenum );
  G4double GetRefnominal(G4int);
  AerogelType GetAerogelType(G4int);
  G4int GetNumberOfRefIndBins() {return NumberOfRefIndBins; }
  G4double GetAerogelRefphotE(G4int BinNumw ) 
       {return StdAerogelRefphotE[BinNumw]; }
  G4double GetStdAerogelRefIndValue(G4int BinNumv ) 
       {return StdAerogelRefIndexValue[BinNumv];} 
  vector<G4double> GetAerogelRefphotEVect() 
         { return StdAerogelRefphotE ; }
  vector<G4double> GetStdAerogelRefIndValueVect() 
         {return StdAerogelRefIndexValue ; }

  G4double GetAerogelRefIndShift(G4int iAgtype) 
              {return  AerogelRefIndShift[iAgtype];}
  G4double GetAerogelWavelengthRef(G4int jAgtype) 
              {return AerogelWavelengthRef[jAgtype]; }
  G4double  GetStdAgelRefIndexE(G4double RefEnergy );
  G4String GetAerogelRefIndexFileName() {return AerogelRefIndexFileName; }
private:

  G4String AerogelRefIndexFileName;
  G4int NumberOfRefIndBins;
  vector<G4double>StdAerogelRefphotE ;
  vector<G4double>StdAerogelRefIndexValue ;
  G4double StdAerogelNominalRefractiveIndex;
  vector<G4double> AerogelRefIndShift;
  vector<G4double> AerogelWavelengthRef;
};

#endif 

