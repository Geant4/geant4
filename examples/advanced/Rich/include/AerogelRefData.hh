//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
  std::vector<G4double> GetCurAerogelRefIndValueVect(G4int AerogelTypenum );
  G4double GetCurAerogelRefIndValue( G4int rbinw , G4int AerogelTypenum );
  G4double GetRefnominal(G4int);
  AerogelType GetAerogelType(G4int);
  G4int GetNumberOfRefIndBins() {return NumberOfRefIndBins; }
  G4double GetAerogelRefphotE(G4int BinNumw ) 
       {return StdAerogelRefphotE[BinNumw]; }
  G4double GetStdAerogelRefIndValue(G4int BinNumv ) 
       {return StdAerogelRefIndexValue[BinNumv];} 
  std::vector<G4double> GetAerogelRefphotEVect() 
         { return StdAerogelRefphotE ; }
  std::vector<G4double> GetStdAerogelRefIndValueVect() 
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
  std::vector<G4double>StdAerogelRefphotE ;
  std::vector<G4double>StdAerogelRefIndexValue ;
  G4double StdAerogelNominalRefractiveIndex;
  std::vector<G4double> AerogelRefIndShift;
  std::vector<G4double> AerogelWavelengthRef;
};

#endif 

