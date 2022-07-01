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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eDPWAElasticDCS
//
// Author:        Mihaly Novak
//
// Creation date: 02.07.2020
//
// Modifications:
//
// Class Description:
//
// Contains numerical Differential Cross Sections (DCS) for e-/e+ Coulomb
// scattering computed by Dirac Partial Wave Analysis (DPWA) [1]:
// - electrostatic interaction, with a local exchange correction in the case of
//   electrons (using Dirac-Fock e- densities; finite nuclear size with Fermi
//   charge distribution; exchange potential with Furness and McCarthy for e-)[2]
// - correlation-polarization (projectiles cause the polarization of the charge
//   cloud of the target atom and the induced dipole moment acts back on the
//   projectile) was accounted by using Local-Density Approximation (LDA) [2]
// - absorption: not included since it's an inelastic channel [2] (the cor-
//   responding excitations needs to be modelled by a separate, independent,
//   inelastic model).
// Using the above mentioned DPWA computation with a free atom approximation
// might lead to questionable results below few hundred [eV] where possible
// solid state or bounding effects might start to affect the potential.
// Nevertheless, the lower energy was set to 10 eV in order to provide(at least)
// some model even at low energies (with this caution). The highest projectile
// kinetic energy is 100 [MeV].
//
// The class provides interface methods for elastic, first-, second-transport
// cross section computations as well as for sampling cosine of polar angular
// deflections. These interface methods are also available for resricted cross
// section computations and angular deflection sampling.
//
// References:
//
// [1] Salvat, F., Jablonski, A. and Powell, C.J., 2005. ELSEPA—Dirac partial-
//     wave calculation of elastic scattering of electrons and positrons by
//     atoms, positive ions and molecules. Computer physics communications,
//     165(2), pp.157-190.
// [2] Salvat, F., 2003. Optical-model potential for electron and positron
//     elastic scattering by atoms. Physical Review A, 68(1), p.012708.
// [3] Benedito, E., Fernández-Varea, J.M. and Salvat, F.,2001. Mixed simulation
//     of the multiple elastic scattering of electrons and positrons using
//     partial-wave differential cross-sections. Nuclear Instruments and Methods
//     in Physics Research Section B: Beam Interactions with Materials and Atoms,
//     174(1-2), pp.91-110.
//
// -------------------------------------------------------------------

#ifndef G4eDPWAElasticDCS_h
#define G4eDPWAElasticDCS_h 1


#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "globals.hh"
#include "G4String.hh"
#include "G4Physics2DVector.hh"


#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"


#include "G4Log.hh"
#include "G4Exp.hh"


class G4eDPWAElasticDCS {

public:

  // CTR:
  // - iselectron   : data for e- (for e+ otherwise)
  // - isrestricted : sampling of angular deflection on restricted interavl is
  //                  required (i.e. in case of mixed-simulation models)
  G4eDPWAElasticDCS(G4bool iselectron=true, G4bool isrestricted=false);

   // DTR
 ~G4eDPWAElasticDCS();

  // initialise for a given 'iz' atomic number:
  //  - nothing happens if it has already been initialised for that Z.
  void InitialiseForZ(std::size_t iz);

  // Computes the elastic, first and second cross sections for the given kinetic
  // energy and target atom.
  // Cross sections are zero ff ekin is below/above the kinetic energy grid
  void  ComputeCSPerAtom(G4int iz, G4double ekin, G4double& elcs, G4double& tr1cs,
                        G4double& tr2cs, G4double mumin=0.0, G4double mumax=1.0);


  // samples cos(theta) i.e. cosine of the polar angle of scattering in elastic
  // interaction (Coulomb scattering) of the projectile (e- or e+ depending on
  // fIsElectron) with kinetic energy of exp('lekin'), target atom with atomic
  // muber of 'iz'. See the 'SampleCosineThetaRestricted' for obtain samples on
  // a restricted inteval.
  G4double SampleCosineTheta(std::size_t iz, G4double lekin, G4double r1,
                             G4double r2, G4double r3);

  // samples cos(theta) i.e. cosine of the polar angle of scattering in elastic
  // interaction (Coulomb scattering) of the projectile (e- or e+ depending on
  // fIsElectron) with kinetic energy of exp('lekin'), target atom with atomic
  // muber of 'iz'.
  // The cosine theta will be in the [costMin, costMax] interval where costMin
  // corresponds to a maximum allowed polar scattering angle thetaMax while
  // costMin corresponds to minimum allowed polar scatterin angle thetaMin.
  // See the 'SampleCosineTheta' for obtain samples on the entire [-1,1] range.
  G4double SampleCosineThetaRestricted(std::size_t iz, G4double lekin,
                                       G4double r1, G4double r2,
                                       G4double costMax, G4double costMin);

  // interpolate scattering power correction form table buit at init.
  G4double ComputeScatteringPowerCorrection(const G4MaterialCutsCouple *matcut,
                                            G4double ekin);

  // build scattering power correction table at init.
  void     InitSCPCorrection(G4double lowEnergyLimit, G4double highEnergyLimit);


private:

  // data structure to store one sampling table: combined Alias + RatIn
  // NOTE: when Alias is used, sampling on a resctricted interval is not possible
  //       However, Alias makes possible faster sampling. Alias is used in case
  //       of single scattering model while it's not used in case of mixed-model
  //       when restricted interval sampling is needed. This is controlled by
  //       the fIsRestrictedSamplingRequired flag (false by default).
  struct OneSamplingTable {
    OneSamplingTable () = default;
    void SetSize(std::size_t nx, G4bool useAlias)  {
      fN = nx;
      // Alias
      if (useAlias) {
        fW.resize(nx);
        fI.resize(nx);
      }
      // Ratin
      fCum.resize(nx);
      fA.resize(nx);
      fB.resize(nx);
    }

    // members
    std::size_t           fN;            // # data points
    G4double              fScreenParA;   // the screening parameter
    std::vector<G4double> fW;
    std::vector<G4double> fCum;
    std::vector<G4double> fA;
    std::vector<G4double> fB;
    std::vector<G4int>    fI;
  };


  // loads the kinetic energy and theta grids for the DCS data (first init step)
  // should be called only by the master
  void LoadGrid();

  // load DCS data for a given Z
  void LoadDCSForZ(G4int iz);

  // loads sampling table for the given Z over the enrgy grid
  void BuildSmplingTableForZ(G4int iz);

  G4double SampleMu(std::size_t izet, std::size_t ie, G4double r1, G4double r2);

  G4double FindCumValue(G4double u, const OneSamplingTable& stable,
                        const std::vector<G4double>& uvect);

  // muMin and muMax : no checks on these
  G4double SampleMu(std::size_t izet, std::size_t ie, G4double r1, G4double muMin,
                    G4double muMax);

  // set the DCS data directory path
  const G4String& FindDirectoryPath();

  // uncompress one data file into the input string stream
  void ReadCompressedFile(G4String fname, std::istringstream &iss);

  // compute Molier material dependent parameters
  void ComputeMParams(const G4Material* mat, G4double& theBc, G4double& theXc2);


// members
private:

  // indicates if the object is for mixed-simulation (single scatterin otherwise)
  G4bool                           fIsRestrictedSamplingRequired;
  // indicates if the object is for e- (for e+ otherwise)
  G4bool                           fIsElectron;
  // indicates if the ekin, mu grids has already been loaded (only once)
  static G4bool                    gIsGridLoaded;
  // data directory
  static G4String                  gDataDirectory;
  // max atomic number (Z) for which DCS has been computed (103)
  static constexpr std::size_t     gMaxZ = 103;
  // energy and theta grid(s) relaed variables: loaded from gridinfo by LoadGrid
  static std::size_t               gNumEnergies;
  static std::size_t               gIndxEnergyLim;// the energy index just above 2 [keV]
  static std::size_t               gNumThetas1;   // used for e- below 2 [keV]
  static std::size_t               gNumThetas2;   // used for e+ and for e- bove 2 [keV]
  static std::vector<G4double>     gTheEnergies;  // log-kinetic energy grid
  static std::vector<G4double>     gTheMus1;      // mu(theta) = 0.5[1-cos(theta)]
  static std::vector<G4double>     gTheMus2;
  static std::vector<G4double>     gTheU1;        // u(mu; A'=0.01) = (A'+1)mu/(mu+A')
  static std::vector<G4double>     gTheU2;
  static G4double                  gLogMinEkin;   // log(gTheEnergies[0])
  static G4double                  gInvDelLogEkin;// 1./log(gTheEnergies[i+1]/gTheEnergies[i])
  // abscissas and weights of an 8 point Gauss-Legendre quadrature
  // for numerical integration on [0,1]
  static const G4double            gXGL[8];
  static const G4double            gWGL[8];
  //
  std::vector<G4Physics2DVector*>  fDCS;    // log(DCS) data per Z
  std::vector<G4Physics2DVector*>  fDCSLow; // only for e- E < 2keV
  // sampling tables: only one of the followings will be utilized
  std::vector< std::vector<OneSamplingTable>* > fSamplingTables;
  //
  // scattering power correction: to account sub-threshold inelastic deflections
  const G4int                      fNumSPCEbinPerDec = 3;
  struct SCPCorrection {
    G4bool   fIsUse;               //
    G4double fPrCut;               // sec. e- production cut energy
    G4double fLEmin;               // log min energy
    G4double fILDel;               // inverse log delta kinetic energy
    std::vector<G4double> fVSCPC;  // scattering power correction vector
  };
  std::vector<SCPCorrection*>      fSCPCPerMatCuts;


};

#endif
