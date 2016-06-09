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
//
// $Id: G4StrawTubeXTRadiator.hh,v 1.1 2005/04/22 09:44:18 grichine Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
//


#ifndef G4StrawTubeXTRadiator_h
#define G4StrawTubeXTRadiator_h 1

#include <complex>
#include "G4VXTRenergyLoss.hh"

class G4StrawTubeXTRadiator : public G4VXTRenergyLoss
{
public:

  G4StrawTubeXTRadiator (G4LogicalVolume* anEnvelope, G4Material*, G4Material*,
                        G4double,G4double,G4Material*,G4bool unishut = false,
                        const G4String & processName = "StrawTubeXTRadiator");
  ~G4StrawTubeXTRadiator ();

// Auxiliary functions for plate/gas material parameters

  G4double  GetMediumFormationZone(G4double,G4double,G4double) ;
  void      ComputeMediumPhotoAbsCof() ;
  G4double  GetMediumLinearPhotoAbs(G4double) ;
  G4complex GetMediumComplexFZ(G4double,G4double,G4double) ;





  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);


protected:

  G4int      fMatIndex3;
  G4double   fSigma3;
  G4double** fMediumPhotoAbsCof ;
  G4int      fMediumIntervalNumber ;
  
};

#endif
