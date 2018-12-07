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
/// \file runAndEvent/RE02/include/RE02RunAction.hh
/// \brief Definition of the RE02RunAction class
//
//
// 
//  
//

#ifndef RE02RunAction_h
#define RE02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

//=======================================================================
// RE02RunAction
//   
//  Generate Run object and Dumping Run summary.
//
//  T.Aso Created. 2007.Nov.
//
//=======================================================================
//
/// User run action class
///
/// - G4Run* GenerateRun()
///     instanciates a run of RE02Run with a sensitive detector name
///
/// - void BeginOfRunAction(const G4Run*)
///     shows the run number 
///
/// - void EndOfRunAction(const G4Run*)
///     shows accumulated information of primitive scorers and
///     outputs total energy deposit into "totED.txt" file
///
/// - G4int CopyNo(G4int ix, G4int iy, G4int iz)
///     returns a copy number of the segment number of the water phantom
///     indicated by three dimentional indexes
//
class RE02RunAction : public G4UserRunAction
{
public:
  // constructor and destructor
  RE02RunAction();
  virtual ~RE02RunAction();

public:
  // virtual method from G4UserRunAction.
  virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

public:
  // Utility method for converting segment number of
  // water phantom to copyNo of HitsMap.
  G4int CopyNo(G4int ix, G4int iy, G4int iz)
  {  return (iy*(fNx*fNz)+ix*fNz+iz); }

private:
  // Data member 
  // - vector of MultiFunctionalDetecor names.
  std::vector<G4String> fSDName;  

  // for conversion of sengment number to copyNo.
  G4int fNx, fNy, fNz;

};

//

#endif





