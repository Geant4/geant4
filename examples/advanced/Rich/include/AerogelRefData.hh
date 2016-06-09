//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

