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
// File name:     G4GSMottCorrection
//
// Author:        Mihaly Novak
//
// Creation date: 23.08.2017
//
// Modifications:
//
// Class description:
//   An object of this calss is used in the G4GoudsmitSaundersonTable when Mott-correction
//   was required by the user in the G4GoudsmitSaundersonMscModel.
//   The class is responsible to handle pre-computed Mott correction (rejection) functions
//   obtained as a ratio of GS angular distributions computed based on the Screened-Rutherford
//   DCS to GS angular distributions computed based on a more accurate corrected DCS_{cor}.
//   The DCS used to compute the accurate Goudsmit-Saunderson angular distributions is [1]:
//   DCS_{cor} = DCS_{SR}x[ DCS_{R}/DCS_{Mott}] where :
//    # DCS_{SR} is the relativistic Screened-Rutherford DCS (first Born approximate
//      solution of the Klein-Gordon i.e. relativistic Schrodinger equation =>
//      scattering of spinless e- on exponentially screened Coulomb potential)
//      note: the default (without using Mott-correction) GS angular distributions
//      are based on this DCS_{SR} with Moliere's screening parameter!
//    # DCS_{R} is the Rutherford DCS which is the same as above but without
//      screening
//    # DCS_{Mott} is the Mott DCS i.e. solution of the Dirac equation with a bare
//      Coulomb potential i.e. scattering of particles with spin (e- or e+) on a
//      point-like unscreened Coulomb potential [2]
//    # moreover, the screening parameter of the DCS_{cor} was determined such that
//  the DCS_{cor} with this corrected screening parameter reproduce the first
//  transport cross sections obtained from the corresponding most accurate DCS [3].
//  Unlike the default GS, the Mott-corrected angular distributions are particle type
//  (different for e- and e+ <= the DCS_{Mott} and the screening correction) and target
//  (Z and material) dependent.
//
// References:
//   [2] I.Kawrakow, E.Mainegra-Hing, D.W.O.Rogers, F.Tessier,B.R.B.Walters, NRCC
//       Report PIRS-701 (2013)
//   [2]  N.F. Mott, Proc. Roy. Soc. (London) A 124 (1929) 425.
//   [3] F.Salvat, A.Jablonski, C.J. Powell, CPC 165(2005) 157-190
//
// -----------------------------------------------------------------------------

#ifndef G4GSMottCorrection_h
#define G4GSMottCorrection_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"

#include <vector>
#include <string>
#include <sstream>

class G4Material;
class G4Element;


class G4GSMottCorrection {
public:
  G4GSMottCorrection(G4bool iselectron=true);

 ~G4GSMottCorrection();

  void     Initialise();

  void     GetMottCorrectionFactors(G4double logekin, G4double beta2, G4int matindx,
                                    G4double &mcToScr, G4double &mcToQ1, G4double &mcToG2PerG1);

  G4double GetMottRejectionValue(G4double logekin, G4double G4beta2, G4double q1, G4double cost,
                                 G4int matindx, G4int &ekindx, G4int &deltindx);

  static G4int GetMaxZet() { return gMaxZet; }

private:
  void InitMCDataPerElement();

  void InitMCDataPerMaterials();

  void LoadMCDataElement(const G4Element*);

  void ReadCompressedFile(std::string fname, std::istringstream &iss);

  void InitMCDataMaterial(const G4Material*);
  //
  // dat structures
  struct DataPerDelta {
    G4double         fSA;             // a,b,c,d spline interpolation parameters for the last \sin(0.5\theta) bin
    G4double         fSB;
    G4double         fSC;
    G4double         fSD;
    G4double        *fRejFuntion;     // rejection func. for a given E_{kin}, \delta, e^-/e^+ over the \sin(0.5\theta) grid
  };

  struct DataPerEkin {
    G4double         fMCScreening;    // correction factor to Moliere screening parameter
    G4double         fMCFirstMoment;  // correction factor to first moment
    G4double         fMCSecondMoment; // correction factor to second
    DataPerDelta   **fDataPerDelta;   // per delta value data structure for each delta values
  };

  // either per material or per Z
  struct DataPerMaterial {
    DataPerEkin  **fDataPerEkin;    // per kinetic energy data structure for each kinetic energy value
  };
  //
  void AllocateDataPerMaterial(DataPerMaterial*);
  void DeAllocateDataPerMaterial(DataPerMaterial*);
  void ClearMCDataPerElement();
  void ClearMCDataPerMaterial();
  //
  // data members:
  // - Mott correction data are computed over a :
  //  I.  Kinetic energy grid [both rejection functions and correction factors]:
  //      1. kinetic energy grid from 1[keV] - 100[keV] with log-spacing 16 points:
  //                 # linear interpolation on \ln[E_{kin}] will be used
  //      2. \beta^2 grid from E_{kin} = 100[keV](~0.300546) - \beta^2=0.9999(~50.5889MeV]) with linear spacing 16 points:
  //                 # linear interpolation on \beta^2 will be used
  //      3. the overall kinetic energy grid is from E_{kin}=1[keV] - E_{kin}<=\beta^2=0.9999(~50.5889MeV]) with 31 points
  //  II. Delta value grid [rejection functions at a given kinetic energy(also depends on \theta;Z,e-/e+)]:
  //      1. \delta=2 Q_{1SR} (\eta_{MCcor})/ [1-2 Q_{1SR} (\eta_{MCcor})] where Q_{1SR} is the first moment i.e.
  //         Q_{1SR}(\eta_{MCcor}) =s/\lambda_{el}G_{1SR}(\eta_{MCcor}) where s/\lambda_{el} is the mean number of elastic
  //         scattering along the path s and G_{1SR}(\eta_{MCcor}) is the first, Screened-Rutherford transport coefficient
  //         but computed by using the Mott-corrected Moliere screening parameter
  //      2. the delta value grid is from [0(1e-3) - 0.9] with linear spacing of 28 points:
  //                 # linear interpolation will be used on \delta
  // III. \sin(0.5\theta) grid[rejection function at a given kinetic energy - delta value pair (also depends on Z,e-/e+)]:
  //      1. 32 \sin(0.5\theta) pints between [0,1] with linear spacing: # linear interpolation on \sin(0.5\theta) will
  //         be used exept the last bin where spline is used (the corresponding 4 spline parameters are also stored)
private:
  G4bool                     fIsElectron;
  static constexpr G4int     gNumEkin   = 31;                 // number of kinetic energy grid points for Mott correction
  static constexpr G4int     gNumBeta2  = 16;                 // \beta^2 values between [fMinBeta2-fMaxBeta2]
  static constexpr G4int     gNumDelta  = 28;                 // \delta values between [0(1.e-3)-0.9]
  static constexpr G4int     gNumAngle  = 32;                 //
  static constexpr G4int     gMaxZet    = 98;                 // max. Z for which Mott-correction data were computed (98)
  static constexpr G4double  gMinEkin   =   1.*CLHEP::keV;   // minimum kinetic energy value
  static constexpr G4double  gMidEkin   = 100.*CLHEP::keV;   // kinetic energy at the border of the E_{kin}-\beta^2 grids
  static constexpr G4double  gMaxBeta2  =   0.9999;           // maximum \beta^2 value
  static constexpr G4double  gMaxDelta  =   0.9;              // maximum \delta value (the minimum is 0(1.e-3))
  //
  G4double                   fMaxEkin;        // from max fMaxBeta2 = 0.9999 (~50.5889 [MeV])
  G4double                   fLogMinEkin;     // \ln[fMinEkin]
  G4double                   fInvLogDelEkin;  // 1/[\ln(fMidEkin/fMinEkin)/(fNumEkin-fNumBeta2)]
  G4double                   fMinBeta2;       // <= E_{kin}=100 [keV] (~0.300546)
  G4double                   fInvDelBeta2;    // 1/[(fMaxBeta2-fMinBeta2)/(fNumBeta2-1)]
  G4double                   fInvDelDelta;    // 1/[0.9/(fNumDelta-1)]
  G4double                   fInvDelAngle;    // 1/[(1-0)/fNumAngle-1]
  //
  static const std::string   gElemSymbols[];
  //
  std::vector<DataPerMaterial*>  fMCDataPerElement;   // size will be gMaxZet+1; won't be null only at used Z indices
  std::vector<DataPerMaterial*>  fMCDataPerMaterial;  // size will #materials; won't be null only at used mat. indices
};

#endif // G4GSMottCorrection_h
