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
//
// $Id: G4EvaporationLevelDensityParameter.hh,v 1.4.2.1 2009/03/04 14:56:06 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4EvaporationLevelDensityParameter_h
#define G4EvaporationLevelDensityParameter_h 1


#include "G4VLevelDensityParameter.hh"
#include "G4CameronTruranHilfShellCorrections.hh"

class G4EvaporationLevelDensityParameter : public G4VLevelDensityParameter
{
public:

  G4EvaporationLevelDensityParameter();

  virtual ~G4EvaporationLevelDensityParameter();

private:  
	
  G4EvaporationLevelDensityParameter(const G4EvaporationLevelDensityParameter &right);

  const G4EvaporationLevelDensityParameter & operator=(const G4EvaporationLevelDensityParameter &right);
  G4bool operator==(const G4EvaporationLevelDensityParameter &right) const;
  G4bool operator!=(const G4EvaporationLevelDensityParameter &right) const;
  
public:
  G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const;

private:

  G4double ShellCorrection(const G4int Z, const G4int N) const
  { 
    G4CameronTruranHilfShellCorrections* SPtr = G4CameronTruranHilfShellCorrections::GetInstance();
    return SPtr->GetShellZ(Z) + SPtr->GetShellN(N);
  }
		  		
private:

  static const G4double ConstEvapLevelDensityParameter;
  static const G4double alpha;
  static const G4double beta;
  static const G4double gamma;
  static const G4double Bs;
};


#endif
