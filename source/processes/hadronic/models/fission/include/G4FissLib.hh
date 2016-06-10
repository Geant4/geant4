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
// $Id: G4FissLib.hh 67966 2013-03-13 09:38:38Z gcosmo $
//
// ********************************************************************
// !                A neutron-induced fission package                 !
// !                 J.M. Verbeke, Dec-2006 / LLNL                    !
// !                                                                  !
// ! G4NeutronFissionModule.cc                                        !
// !                                                                  !
// ! Classes to simulate neutron-induced fissions, emitting neutrons  !
// ! and gamma-rays. Algorithm uses data whenever available, and      !
// ! models such as the Terrell approximation, the Watt spectrum      !
// ! otherwise.                                                       !
// !                                                                  !
// ! The complete list of references used is shown below:             !
// !                                                                  !
// ! J. Terrell, "Distributions of Fission Neutron Numbers", Phys.    !
// !   Rev. 108, 783 (1957).                                          !
// ! M.S. Zucker, N.E. Holden, "Energy Dependence of Neutron          !
// !   Multiplicity P_nu in Fast-Neutron-Induced Fission for U-235,   !
// !   U-238 and Pu-239," BNL-38491 (1986).                           !
// ! T.E. Valentine, "MCNP-DSP Users Manual," ORNL/TM-13334, R2, Oak  !
// !   Ridge National Laboratory (2000).                              !
// ! T.E. Valentine, J.T. Mihalczo, "MCNP-DSP: A Neutron and Gamma    !
// !   Ray Monte Carlo Calculation of Source-Driven Noise-Measured    !
// !   Parameters ," Ann. of Nucl. Eng., 23, 16, p. 1271 (1996).      !
// ! R. Gwin, R.R. Spencer, R.W. Ingle, "Measurements of the Energy   !
// !   Dependence of Prompt Neutron Emission from U-233, U-235,       !
// !   Pu-239, and Pu-241 for E_n=0.005 to 10 eV Relative to Emission !
// !   from Spontaneous Fission of Cf-252," Nucl. Sci. Eng., 87, 381  !
// !   (1984).                                                        !
// ! J. Frehaut, "Neutron Multiplicity Distribution in Fast           !
// !   Neutron-Induced Fission," Proc. of IAEA Consultant's Meeting   !
// !   on Physics of Neutron Emission in Fission, Mito, Japan (1988). !
// ! R.R. Spencer, R. Gwin, R.W. Ingle, "A measurement of the Average !
// !   Number of Prompt Neutrons from Spontaneous Fission of          !
// !   Californium-252," Nucl. Sci. Eng. 80, 603 (1982).              !
// ! J.W. Boldeman, M.G. Hines, "Prompt Neutron Emission              !
// !   Probabilities Following Spontaneous and Thermal Neutron        !
// !   Fission," Nucl. Sci. Eng., 91, 114 (1985).                     !
// ! N.E. Holden, M.S. Zucker, "A Reevaluation of the Average Prompt  !
// !   Neutron Emission Multiplicity (nubar) Values from Fission of   !
// !   Uranium and Transuranium Nuclides," BNL-NCS-35513, Brookhaven  !
// !   National Laboratory).                                          !
// ! R.J. Howerton, et al, "The LLL Evaluated Nuclear Data Library    !
// !   (ENDL): Evaluation Techniques, Reaction Index, and Description !
// !   of Individual Evaluations," UCRL-50400, V. 15, Part A,         !
// !   Lawrence Livermore National Laboratory (1975).                 !
// ! D.E. Cullen, "Sampling ENDL Watt Fission Spectra,"               !
// !   UCRL-TR-203251, Lawrence Livermore National Laboratory (2004). !
// ! C.J. Everett, E.D. Cashwell, "A Third Monte Carlo Sampler,"      !
// !   LA-9721-MS, Los Alamos National Laboratory (1983).             !
// ! D.E. Cullen, "TART 2002: A Couple Neutron-Photon 3-D,            !
// !   Combinatorial Geometry, Time Dependent Monte-Carlo Transport   !
// !   Code," UCRL-ID-126455, Rev. 4, Lawrence Livermore National     !
// !   Laboratory (2003).                                             !
// ! W. Mannhart, "Evaluation of the Cf-252 Fission Neutron Spectrum  !
// !   Between 0 MeV and 20 MeV," Proc. Advisory Group Mtg. Neutron   !
// !   Sources, Leningrad, USSR, 1986 (IAEA-TECDOC-410), Vienna       !
// !   (1987).                                                        !
// ! D.G. Madland, J.R. Nix, "Prompt Fission Neutron Spectra and      !
// !   Average Prompt Neutron Multiplicities,"NEANDC Specialist's     !
// !   Meeting on Yields and Decay Data of Fission Products,          !
// !   Brookhaven National Laboratory, BNL 51778 (1984).              !
// ! F.H. Froehner, "Evaluation of Cf-252 Prompt Fission Neutron Data !
// !   from 0 to 20 MeV by Watt Spectrum Fit," Nucl. Sci. Eng. 106,   !
// !   345 (1990).                                                    !
// ! G.S. Brunson, Jr., "Multiplicity and Correlated Energy of Gamma  !
// !   Rays Emitted in the Spontaneous Fission of Californium-252,"   !
// !   Ph.D. Thesis, University of Utah (1982).                       !
// ! T.E. Valentine, "Evaluation of Prompt Fission Gamma Rays for Use !
// !   in Simulating Nuclear Safeguard Measurements," Ann. Nucl.      !
// !   Eng., 28, 191 (2001).                                          !
// ! C. Wagemans, "The Nuclear Fission Process," CRC Press, Inc., Boca!
// !   Raton, Florida (1991).                                         !
// ! F.C. Maienschein, R.W. Peelle, T.A. Love, Neutron Phys. Ann.     !
// !   Prog. Rep. for Sept. 1, 1958, ORNL-2609, Oak Ridge National    !
// !   Laboratory (1958).                                             !
// ! "Fundamental Aspects of Reactor Shielding," Addison-Wesley       !
// !   Publishing Company, Inc. Reading, Massachussetts (1959).       !
// !                                                                  !
// ********************************************************************
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by J.M. Verbeke, LLNL, 5-Jan-07
 // Builds and has the Cross-section data for one material.
  
#ifndef G4FissLib_h
#define G4FissLib_h 1

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron induced fission below 10 MeV.
// Note that this model (by intent of avoiding the possibility of heating studies) does
// not provide the nuclear fragments.
//
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "globals.hh"
#include "G4ParticleHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4FissionLibrary.hh"
// #include "G4FissLib.hh"

class G4FissLib : public G4HadronicInteraction
{
  public: 
    G4FissLib();
    ~G4FissLib();
  
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& aTargetNucleus);
    const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

  private:
    G4FissionLibrary theLibrary;
  
  private:
    G4double* xSec;
    G4ParticleHPChannel* theFission;
    G4String dirName;
    G4int numEle;
};

#endif
