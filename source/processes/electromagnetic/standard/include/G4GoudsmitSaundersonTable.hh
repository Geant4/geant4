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
// $Id: G4GoudsmitSaundersonTable.hh 93663 2015-10-28 09:50:49Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4GoudsmitSaundersonTable
//
// Author:        Mihaly Novak / (Omrane Kadri)
//
// Creation date: 20.02.2009
//
// Class description:
//   Class to handle multiple scattering angular distributions precomputed by
//   using Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened
//   Rutherford DCS for elastic scattering of electrons/positrons [1,2]. This
//   class is used by G4GoudsmitSaundersonMscModel to sample the angular
//   deflection of electrons/positrons after travelling a given path.
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 18.05.2015 M. Novak This class has been completely replaced (only the original
//            class name was kept; class description was also inserted):
//            A new version of Kawrakow-Bielajew Goudsmit-Saunderson MSC model
//            based on the screened Rutherford DCS for elastic scattering of
//            electrons/positrons has been introduced[1,2]. The corresponding MSC
//            angular distributions over a 2D parameter grid have been recomputed
//            and the CDFs are now stored in a variable transformed (smooth) form
//            together with the corresponding rational interpolation parameters.
//            The new version is several times faster, more robust and accurate
//            compared to the earlier version (G4GoudsmitSaundersonMscModel class
//            that use these data has been also completely replaced)
//
// References:
//   [1] A.F.Bielajew, NIMB, 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//
// -----------------------------------------------------------------------------

#ifndef G4GoudsmitSaundersonTable_h
#define G4GoudsmitSaundersonTable_h 1

#include <vector>

#include "G4Types.hh"

class G4GoudsmitSaundersonTable
{
public:

  G4GoudsmitSaundersonTable(){};
  ~G4GoudsmitSaundersonTable();

  // initialie:
  //   - loads the precomputed MSC angular CDFs into memory
  //   - init. material dependent MSC parameters (Moliere's screening)
  //  (- only Master thread and only once)
  void Initialise();

  // samples cos(theta) i.e. angular deflection from the precomputed angular
  // distributions in the real multiple scattering case
  G4double SampleCosTheta(G4double, G4double, G4double, G4double, G4double, G4double);
  G4double SampleCosThetaII(G4double, G4double, G4double, G4double, G4double, G4double);

  // returns with the screening parameter value that results with the first
  // transport coefficient (G1) received as input parameter according to the
  // screened Rutherford DCS. Used only when fgIsUsePWATotalXsecData is TRUE
  // in G4GoudsmitSaundersonMscModel i.e. when PWA screeing is used instead of
  // Moliere's one.
  G4double GetScreeningParam(G4double);

  // samples angular deflection cos(theta) and sin(theta) for electrons/positrons
  // involving sampling of no scattering, single scattering, "few" scattering and
  // real multiple scattering
  void Sampling(G4double, G4double, G4double,  G4double&, G4double&);

  // material dependent MSC parameters (computed at initialisation) regarding
  // Moliere's screening parameter
  G4double GetMoliereBc(G4int matindx){return (*fgMoliereBc)[matindx];}
  G4double GetMoliereXc2(G4int matindx){return (*fgMoliereXc2)[matindx];}

private:

  //  hide assignment operator and cpy ctr.
  G4GoudsmitSaundersonTable & operator=(const  G4GoudsmitSaundersonTable &right);
  G4GoudsmitSaundersonTable(const  G4GoudsmitSaundersonTable&);

  // load precomputed CDFs of MSC angular distributions over a 2D parameter grid
  // CDFs are stored in a variable transformed, equally probable intervall form
  // together with the corresponding rational interpolation paraneters
  void LoadMSCData();
  void LoadMSCDataII();

  // initialisation of material dependent Moliere's MSC parameters
  void InitMoliereMSCParams();

private:
   //@{
   /** size of grids of some parameters */
   static const G4int fgNumLambdas =  76;  /** number of \f$ s/\lambda_{e} $\f-values      */
   static const G4int fgNumLamG1   =  21;  /** number of \f$ s/\lambda_{e}G_{1} $\f-values */
   static const G4int fgNumLamG1II   =  22;  /** number of \f$ s/\lambda_{e}G_{1} $\f-values */
   static const G4int fgNumUvalues = 101;  /** number of u-vaues                           */
   static const G4int fgNumScreeningParams = 160; /** number of A-vaues                    */
   //@}

   //@{
   /** girds of fixed parameter values */
   /** the grid \f$ s/\lambda_{e} $\f-values; size = fgNumLambdas = 76        */
   static const G4double fgLambdaValues[];
   /** the grid of \f$ s/\lambda_{e}G_{1} $\f-values; size = fgNumLamG1 = 11 */
   static const G4double fgLamG1Values[];
   static const G4double fgLamG1ValuesII[];

   /** the grid of u-values; size = fgNumUvalues = 101 */
   static const G4double fgUValues[];
   //@}

   // precomputed G1(A) function as a table -> run time interpolation to determine
   // the screening parameter value A that gives back the given first transport
   // coefficient G1
   static const G4double fgG1Values[];
   static const G4double fgScreeningParam[];
   static const G4double fgSrcAValues[];
   static const G4double fgSrcBValues[];
   //@{
   /** Precomputed equaly probable inverse CDF-s over the 3D parameter grid plus
    *  precomputed parameters necessary for proper rational interpolation of the
    *  inverse CDF.
    */
   static G4double fgInverseQ2CDFs[fgNumLambdas*fgNumLamG1*fgNumUvalues];
   static G4double fgInterParamsA2[fgNumLambdas*fgNumLamG1*fgNumUvalues];
   static G4double fgInterParamsB2[fgNumLambdas*fgNumLamG1*fgNumUvalues];
   static G4double fgInverseQ2CDFsII[fgNumLambdas*fgNumLamG1II*fgNumUvalues];
   static G4double fgInterParamsA2II[fgNumLambdas*fgNumLamG1II*fgNumUvalues];
   static G4double fgInterParamsB2II[fgNumLambdas*fgNumLamG1II*fgNumUvalues];

   //@}


   //@{
   /** Precomputed \f$ b_lambda_{c} $\f and \f$ \chi_c^{2} $\f material dependent
   *   Moliere parameters that can be used to compute the screening parameter,
   *   the elastic scattering cross section (or \f$ \lambda_{e} $\f) under the
   *   screened Rutherford cross section approximation. (These are used in
   *   G4GoudsmitSaundersonMscModel if fgIsUsePWATotalXsecData is FALSE.)
   */
   static std::vector<G4double> *fgMoliereBc;
   static std::vector<G4double> *fgMoliereXc2;
   //@}

   // flag to check if data are alredy in memory
   static G4bool fgIsInitialised;

};

#endif
