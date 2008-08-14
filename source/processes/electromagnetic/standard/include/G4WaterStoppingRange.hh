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
// $Id: G4WaterStoppingRange.hh,v 1.2 2008-08-14 16:05:53 antoni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4WaterStoppingRange_h
#define G4WaterStoppingRange_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4WaterStoppingRange
//
// Description: Data on stopping powers for light ions in compounds
//
// Author:      A.Ivantchenko 10.07.2008
//
// Modifications:
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the ICRU73 report
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "globals.hh"
#include "G4LPhysicsFreeVector.hh"
#include <vector>


class G4WaterStoppingRange
{
public:

  G4WaterStoppingRange(G4bool splineFlag = true);

  ~G4WaterStoppingRange();

  G4double GetRange(G4int ionZ, G4double kinEnergy);

  G4PhysicsVector* GetPhysicsVector(G4int ionZ);

private:

  void Initialise();

  // hide assignment operator
  G4WaterStoppingRange & operator=(const G4WaterStoppingRange &right);
  G4WaterStoppingRange(const G4WaterStoppingRange&);

  G4bool   spline;
  G4int    Z[16];
  G4double A[16];

  std::vector<G4LPhysicsFreeVector*>  R;
};

#endif
