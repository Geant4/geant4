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
/*
 * File:   G4WattFissionSpectrumValues.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 11, 2011, 11:32 AM
 */

/* * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * *
 *                                                                            *
 *  1.  MCNP - A General Monte carlo N-Particle Transport Code, Version 5,    *
 *      X-5 Monte Carlo Team, Volume I: Overview and Theory, April, 2005      *
 *                                                                            *
 * * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * */

#ifndef G4WATTFISSIONSPECTRUMVALUES_HH
#define	G4WATTFISSIONSPECTRUMVALUES_HH

#include "globals.hh"

#include "G4FFGDefaultValues.hh"
#include "G4FFGEnumerations.hh"

// TODO Migrate to existing neutron_hp watt constants in G4NeutronHPWattSpectrum.hh
//      and then remove this file from the repo and sources.cmake

/** WattSpectrumConstants contains constants and other variables for use in
 *  sampling the Watt fission spectrum.
 */
struct WattSpectrumConstants
{
    /** Isotope code in ZZZAAA format for which the Watt fission
     *  spectrum is being sampled
     */
    G4int Product;
    /** Fission cause for which the Watt fission spectrum is being
     * sampled
     */
    G4FFGEnumerations::FissionCause Cause;
    /** Energy, if any, of the incident particle that cause the fission */
    G4double Energy;

    /** Sampling constant. Calculated as:
     * \f[
     *  L = \frac{[K + (K^2 - 1)^\frac{1}{2}]}{a}
     * \f]
     * \f[
     * K = 1 + \frac{b}{8a}
     * \f]
     */
    G4double L;
    /** Sampling constant. Calculated as:
     * \f[
     * M = a*L-1
     * \f]
     */
    G4double M;
    /** Sampling constant taken from the data tables. */
    G4double B;
};

/** These are the energy values in MeV for the neutron induced Watt fission
 *  spectrum constants.
 */
static const G4double IncidentEnergyBins[] = 
{
    G4FFGDefaultValues::ThermalNeutronEnergy,
    1.0 * CLHEP::MeV,
    14.0 * CLHEP::MeV,
    -1 // End of array
};

/** Watt fission spectrum constants for neutron induced fission.
 *  \n <b> Constants </b>
 *  \n Column 1: 'a' value
 *  \n Column 2: 'b' value
 *
 *  \n <b> Incident Neutron Energies </b>
 *  \n Row 1: Thermal (~0.025 eV)
 *  \n Row 2: 1 MeV
 *  \n Row 3: 14 MeV
 */
static const G4double NeutronInducedWattConstants[][3][2] =
{
// Default
    { {0.95,    2.7},
      {1.0,     2.5},
      {1.05,    2.4}, },
// Thorium
    // 90232
    { {1.0888,  1.6871},
      {1.1096,  1.6316},
      {1.1700,  1.4610}, },
// Uranium
    // 92233
    { {0.977,   2.546},
      {0.977,   2.249},
      {1.0036,  2.6377}, },
    // 92235
    { {0.988,   2.249},
      {0.988,   2.249},
      {1.028,   2.084}, },
    // 92238
    { {0.88111, 3.4005},
      {0.89506, 3.2953},
      {0.96534, 2.8330}, },
// Plutonium
    // 94239
    { {0.966,   2.842},
      {0.966,   2.842},
      {1.055,   2.383}, }
};

/** This table provides the indexing for NeutronInducedWattConstants_. The
 *  index of an isotope in this table is the index for the Watt fission spectrum
 *  constants in NeutronInducedWattConstants_. The isotopes are listed in ZZZAAA
 *  format.
 */
static const G4int NeutronInducedWattIsotopesIndex[] =
{
// Default
    0,
// Thorium
    90232,
// Uranium
    92233,
    92235,
    92238,
// Plutonium
    94239,
// End of array
    -1
};

/** Watt fission spectrum constants for spontaneous fission.
 *  \n Column 1: 'a' value
 *  \n Column 2: 'b' value
 */
static const G4double SpontaneousWattConstants[][2] =
{
// Default
    {0.8,      4.0},
// Plutonium
    // 94240
    {0.799,    4.903},
    // 94242
    {0.833668, 4.431658},
// Curium
    // 96242
    {0.891,    4.046},
    // 96244
    {0.906,    3.848},
// Californium
    // 98252
    {1.025,    2.926}
};

/** This table provides the indexing for SpontaneousWattConstants_. The index of
 *  an isotope in this table is the index for the Watt fission spectrum constants
 *  in SpontaneousWattConstants_. The isotopes are listed in ZZZAAA format.
 */
static const G4int SpontaneousWattIsotopesIndex[] =
{
// Default
    0,
// Plutonium
    94240,
    94242,
// Curium
    96242,
    96244,
// Californium
    98252,
// End of array
    -1
};

#endif	/* G4WATTFISSIONSPECTRUMVALUES_HH */

