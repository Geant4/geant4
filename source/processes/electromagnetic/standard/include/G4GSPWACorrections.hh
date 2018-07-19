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
// $Id: $
//
// ----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4GSPWACorrections
//
// Author:        Mihaly Novak
//
// Creation date: 17.10.2017
//
// Modifications:
//
// Class description: class to describe and store correction factors to the
//   integrated quantities of G4GoudsmitSaundersonMscModel (screening parameter,
//   first and second moments) derived by using accurate Dirac-PWA based
//   integrated quantities.
//
// ----------------------------------------------------------------------------

#ifndef G4GSPWACorrections_h
#define G4GSPWACorrections_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"

#include <vector>
#include <string>
#include <sstream>

class G4Material;
class G4Element;


class G4GSPWACorrections {
public:
  G4GSPWACorrections(G4bool iselectron=true);

 ~G4GSPWACorrections();

  void     Initialise();

  void     GetPWACorrectionFactors(G4double logekin, G4double beta2, G4int matindx,
                                   G4double &corToScr, G4double &corToQ1, G4double &corToG2PerG1);
private:
  void     InitDataPerElement();

  void     InitDataPerMaterials();

  void     LoadDataElement(const G4Element*);

  void     InitDataMaterial(const G4Material*);

  void     ClearDataPerElement();

  void     ClearDataPerMaterial();

  // either per material or per Z
  struct DataPerMaterial {
    std::vector<G4double>   fCorScreening;    // correction factor to Moliere screening parameter
    std::vector<G4double>   fCorFirstMoment;  // correction factor to first moment
    std::vector<G4double>   fCorSecondMoment; // correction factor to second
  };


// data members
private:
  G4bool   fIsElectron;
  static constexpr G4int     gMaxZet    = 98;                 // max. Z for which correction data were computed (98)
  static constexpr G4int     gNumEkin   = 31;                 // number of kinetic energy grid points for Mott correction
  static constexpr G4int     gNumBeta2  = 16;                 // \beta^2 values between [fMinBeta2-fMaxBeta2]
  static constexpr G4double  gMinEkin   =   1.*CLHEP::keV;    // minimum kinetic energy value
  static constexpr G4double  gMidEkin   = 100.*CLHEP::keV;    // kinetic energy at the border of the E_{kin}-\beta^2 grids
  static constexpr G4double  gMaxBeta2  =   0.9999;           // maximum \beta^2 value
  //
  G4double                   fMaxEkin;        // from max fMaxBeta2 = 0.9999 (~50.5889 [MeV])
  G4double                   fLogMinEkin;     // \ln[fMinEkin]
  G4double                   fInvLogDelEkin;  // 1/[\ln(fMidEkin/fMinEkin)/(fNumEkin-fNumBeta2)]
  G4double                   fMinBeta2;       // <= E_{kin}=100 [keV] (~0.300546)
  G4double                   fInvDelBeta2;    // 1/[(fMaxBeta2-fMinBeta2)/(fNumBeta2-1)]
  //
  static const std::string   gElemSymbols[];
  //
  std::vector<DataPerMaterial*>  fDataPerElement;   // size will be gMaxZet+1; won't be null only at used Z indices
  std::vector<DataPerMaterial*>  fDataPerMaterial;  // size will #materials; won't be null only at used mat. indices

};

#endif
