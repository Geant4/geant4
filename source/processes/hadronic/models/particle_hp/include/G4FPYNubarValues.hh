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
 * File:   G4FPYNubarValues.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 11, 2011, 11:23 AM
 */

/* * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * *
 *                                                                            *
 *  1.  "Neutron Multiplicity Measurements of Cf and Fm Isotopes",            *
 *      D. C. Hoffman, G. P. Ford, J. P. Balagna, and L. R. Veeser, Physical  *
 *      Review C, 21.2, February 1980                                         *
 *  2.  "Spontaneous Fission Properties and Lifetime Systematics",            *
 *      D. C. Hoffman, International Conference on Fifty Years Research in    *
 *      Nuclear Fission, Berlin, Germany, 1989                                *
 *  3.  MCNP - A General Monte carlo N-Particle Transport Code, Version 5,    *
 *      X-5 Monte Carlo Team, Volume I: Overview and Theory, April, 2005      *
 *  4.  ENDF database, Search parameters: ENDF/B-VII.0; MF=1; MAT=452         *
 *  5.  Prompt Neutron Multiplicities for the Transplutonium Nuclides,        *
 *      N. E. Holden and M. S. Zucker, Brookhaven National Lab, January, 1985 *
 *                                                                            *
 * * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * */

#ifndef G4FPYNUBARVALUES_HH
#define	G4FPYNUBARVALUES_HH

#include "G4Types.hh"

/** Evaluated nubar values for neutron induced fission. The data are represented
 *  as linear regression fits from experimental data in the form: v = mx + b,
 *  where x is the incident neutron energy.
 *  \n \n Value 1: The isotope in ZZZAAA format
 *  \n Value 2: 'x' of the linear regression fit multiplied by 10^13
 *  \n Value 3: 'b' of the linear regression fit multiplied by 10^4
 *  \n \n See Reference 4 for further information
 */
static const G4int NeutronInducedNubar_[][3] = {
    // Default
        {0, 1476773, 29881},
    // Thorium
        {90227, 818327, 20647}, {90229, 936480, 20872}, {90232, 1252285, 21061},
    // Proactinum
        {91231, 1030251, 22480},
    // Uranium
        {92232, 1070000, 31300}, {92233, 1218868, 25043},
        {92234, 1230006, 23673}, {92235, 1392331, 24336},
        {92236, 1353849, 23717}, {92237, 1502160, 24911},
        {92238, 1343405, 24911},
    // Neptunium
        {93237, 1485790, 26362}, {93238, 1476892, 25150},
    // Plutonium
        {94238, 1479368, 28950}, {94239, 1477500, 28505},
        {94240, 1518313, 28030}, {94241, 1496809, 29443},
        {94242, 1575406, 28100},
    // Americium
        {95241, 1425195, 30827}, {95242, 1337871, 32717},
        {95243, 1305739, 32737},
    // Curium
        {96242, 1720000, 34400}, {96243, 1220064, 34329},
        {96244, 1839123, 32444}, {96245, 1166218, 35968},
        {96246, 1223446, 36153}, {96248, 2080000, 34900},
    // Californium
        {98249, 1796326, 38876}, {98251, 2416144, 41400},
    // Einsteinium
        {99254, 2438676, 40832},
    // Fermium
        {100255, 2499442, 43924},
    // End of array
        {-1, -1, -1} };

/** Recommended Gaussian widths for neutron induced fission.
 *  \n Column 1: The isotope in ZZZAAA format
 *  \n Column 2: The width multiplied by 10^6
 *  \n \n See Reference 3 for further information
 */
// Still need: a lot
static const G4int NeutronInducedNubarWidth_[][2] = {
    // Default
        {0, 1210000},
    // Uranium
        {92233, 1144900}, {92235, 1183744}, {92238, 1245456},
    // Plutonium
        {94239, 1299600}, {94241, 1322500},
    // End of array
        {-1, -1} };

/** Evaluated nubar values for neutron induced fission. The data are represented
 *  as linear regression fits from experimental data in the form: v = mx + b,
 *  where x is the incident neutron energy.
 *  \n \n Value 1: The isotope in ZZZAAA format
 *  \n Value 2: 'x' of the linear regression fit multiplied by 10^13
 *  \n Value 3: 'b' of the linear regression fit multiplied by 10^4
 *  \n \n See References 1, 2, 4, & 5 for further information
 */
static const G4int SpontaneousNubar_[][3] = {
    // Default
        {0, 25000},
    // Uranium
        {92238, 0, 20000},
    // Curium
        {96244, 0, 26875}, {96246, 0, 29480}, {96248, 0, 31500},
    // Californium
        {98250, 0, 25200}, {98252, 0, 37676},
    // Einsteinium
        {99253, 0, 47000},
    // Fermium
        {100254, 0, 39800}, {100256, 0, 37300},
    // End of array
        {-1, -1, -1} };

/** Recommended Gaussian widths for neutron induced fission.
 *  See Reference 3 for further information.
 *  \n Column 1: The isotope in ZZZAAA format
 *  \n Column 2: The width multiplied by 10^6
 */
// Still need: 92238, 99253
static const G4int SpontaneousNubarWidth_[][2] = {
    // Default
        {0, 1210000},
    // Plutonium
        {94238, 1288225}, {94240, 1324801}, {94242, 1347921},
    // Curium
        {96242, 1190281}, {96244, 1216609}, {96246, 1205604}, {96248, 1227664},
    // Californium
        {98250, 1488400}, {98252, 1550025}, {98254, 1476225},
    // Fermium
        {100254, 1552516}, {100256, 1349007}, {100257, 1178983},
    // End of array
        {-1, -1} };

#endif	/* G4FPYNUBARVALUES_HH */

