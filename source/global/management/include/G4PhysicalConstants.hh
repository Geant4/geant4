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
// G4PhysicalConstants
//
// Import CLHEP constants on global namespace.
// Restricted to internal use -only- in source code

// Author: G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4PhysicalConstants_hh
#define G4PhysicalConstants_hh 1

#include <CLHEP/Units/PhysicalConstants.h>

using CLHEP::alpha_rcl2;
using CLHEP::amu;
using CLHEP::amu_c2;
using CLHEP::Avogadro;
using CLHEP::Bohr_radius;
using CLHEP::c_light;
using CLHEP::c_squared;
using CLHEP::classic_electr_radius;
using CLHEP::e_squared;
using CLHEP::electron_charge;
using CLHEP::electron_Compton_length;
using CLHEP::electron_mass_c2;
using CLHEP::elm_coupling;
using CLHEP::epsilon0;
using CLHEP::fine_structure_const;
using CLHEP::h_Planck;
using CLHEP::halfpi;
using CLHEP::hbar_Planck;
using CLHEP::hbarc;
using CLHEP::hbarc_squared;
using CLHEP::k_Boltzmann;
using CLHEP::kGasThreshold;
using CLHEP::mu0;
using CLHEP::neutron_mass_c2;
using CLHEP::pi;
using CLHEP::pi2;
using CLHEP::proton_mass_c2;
using CLHEP::STP_Pressure;
using CLHEP::STP_Temperature;
using CLHEP::twopi;
using CLHEP::twopi_mc2_rcl2;
using CLHEP::Bohr_magneton;
using CLHEP::nuclear_magneton;
using CLHEP::universe_mean_density;

#endif
