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
#ifndef DicomIntersectVolume__HH
#define DicomIntersectVolume__HH

#include "G4UImessenger.hh"

#include <vector>
#include <fstream>
#include "G4PhantomParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PhantomParameterisation;
class G4UIcmdWithAString;
class G4VSolid;

class DicomIntersectVolume : public G4UImessenger
{
public:
  DicomIntersectVolume();
  ~DicomIntersectVolume();

  void SetNewValue(G4UIcommand * command,
		   G4String newValues);

private:

  void BuildUserSolid( std::vector<G4String> params );
  void BuildG4Solid( std::vector<G4String> params );
  G4PhantomParameterisation* GetPhantomParam(G4bool bMustExist);
  G4bool IsPhantomVolume( G4VPhysicalVolume* pv );
  std::vector<G4VPhysicalVolume*> GetPhysicalVolumes( const G4String& name, bool exists, G4int nVols );
  std::vector<G4LogicalVolume*> GetLogicalVolumes( const G4String& name, bool exists, G4int nVols );
  std::vector<G4String> GetWordsInString( const G4String& stemp);

private:
  G4UIcmdWithAString* theUserVolumeCmd;
  G4UIcmdWithAString* theG4VolumeCmd;

  G4VSolid* theSolid;

  G4VPhysicalVolume* thePhantomVolume;

  std::ofstream fout;

  G4PhantomParameterisation* theRegularParam;

  G4ThreeVector thePhantomMinusCorner;

  G4bool* theVoxelIsInside;
};

#endif
