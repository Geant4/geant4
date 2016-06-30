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
// $Id: G4StrawTubeXTRadiator.hh 97385 2016-06-02 09:59:53Z gcosmo $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Process describing a straw tube radiator of X-ray transition radiation.  
// Thicknesses of plates and gas gaps are gamma distributed.
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
// 22.04.05 V. Grichine, first version 
// 28.09.07, V.Ivanchenko general cleanup without change of algorithms
//

#ifndef G4StrawTubeXTRadiator_h
#define G4StrawTubeXTRadiator_h 1

#include <complex>
#include "G4VXTRenergyLoss.hh"

class G4SandiaTable;

class G4StrawTubeXTRadiator : public G4VXTRenergyLoss
{
public:

  explicit G4StrawTubeXTRadiator (G4LogicalVolume* anEnvelope, G4Material*, 
             G4Material*, G4double,G4double,G4Material*,G4bool unishut = false,
             const G4String & processName = "StrawTubeXTRadiator");
  virtual ~G4StrawTubeXTRadiator ();

  // Auxiliary functions for plate/gas material parameters

  G4double  GetMediumFormationZone(G4double,G4double,G4double) ;
  void      ComputeMediumPhotoAbsCof() ;
  G4double  GetMediumLinearPhotoAbs(G4double) ;
  G4complex GetMediumComplexFZ(G4double,G4double,G4double) ;

  // Pure virtual function from base class
  G4double GetStackFactor(G4double energy, G4double gamma, 
                          G4double varAngle) override;

protected:

  G4int      fMatIndex3;
  G4double   fSigma3;

  G4SandiaTable* fMediumPhotoAbsCof;
};

#endif
