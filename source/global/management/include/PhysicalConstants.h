// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PhysicalConstants.h,v 1.4 1999-11-19 09:19:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// ----------------------------------------------------------------------
//
// Class description:
//
// HEP coherent Physical Constants
//
// This file has been provided to CLHEP by Geant4 (simulation toolkit for HEP).
//
// The basic units are :
//  		millimeter  
// 		nanosecond  
// 		Mega electron Volt  
// 		positon charge 
// 		degree Kelvin
//              amount of substance (mole)
//              luminous intensity (candela)
// 		radian  
//              steradian 
//
// Below is a non exhaustive List of Physical CONSTANTS,
// computed in the Internal HEP System Of Units.
//
// Most of them are extracted from the Particle Data Book :
//        Phys. Rev. D  volume 50 3-1 (1994) page 1233
// 
//        ...with a meaningful (?) name ...
//
// You can add your own constants.
//

// Author: M.Maire
//
// History:
//
// 23.02.96 Created
// 26.03.96 Added constants for standard conditions of temperature
//          and pressure; also added Gas threshold.

#ifndef HEP_PHYSICAL_CONSTANTS_H
#define HEP_PHYSICAL_CONSTANTS_H

#include "SystemOfUnits.h"
///#include "G4UnitsTest.hh"

//
//
//
static const HepDouble     pi  = 3.14159265358979323846;
static const HepDouble  twopi  = 2*pi;
static const HepDouble halfpi  = pi/2;
static const HepDouble     pi2 = pi*pi;

//
// 
//
static const HepDouble Avogadro = 6.0221367e+23/mole;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2 
//
static const HepDouble c_light   = 2.99792458e+8 * m/s;
static const HepDouble c_squared = c_light * c_light;

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
static const HepDouble h_Planck      = 6.6260755e-34 * joule*s;
static const HepDouble hbar_Planck   = h_Planck/twopi;
static const HepDouble hbarc         = hbar_Planck * c_light;
static const HepDouble hbarc_squared = hbarc * hbarc;

//
//
//
static const HepDouble electron_charge = - eplus; // see SystemOfUnits.h
static const HepDouble e_squared = eplus * eplus;

//
// amu_c2 - atomic equivalent mass unit
// amu    - atomic mass unit
//
static const HepDouble electron_mass_c2 = 0.51099906 * MeV;
static const HepDouble   proton_mass_c2 = 938.27231 * MeV;
static const HepDouble  neutron_mass_c2 = 939.56563 * MeV;
static const HepDouble           amu_c2 = 931.49432 * MeV;
static const HepDouble              amu = amu_c2/c_squared;

//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
static const HepDouble mu0      = 4*pi*1.e-7 * henry/m;
static const HepDouble epsilon0 = 1./(c_squared*mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
static const HepDouble elm_coupling           = e_squared/(4*pi*epsilon0);
static const HepDouble fine_structure_const   = elm_coupling/hbarc;
static const HepDouble classic_electr_radius  = elm_coupling/electron_mass_c2;
static const HepDouble electron_Compton_length = hbarc/electron_mass_c2;
static const HepDouble Bohr_radius = electron_Compton_length/fine_structure_const;

static const HepDouble alpha_rcl2 = fine_structure_const
                                   *classic_electr_radius
                                   *classic_electr_radius;

static const HepDouble twopi_mc2_rcl2 = twopi*electron_mass_c2
                                             *classic_electr_radius
                                             *classic_electr_radius;
//
//
//
static const HepDouble k_Boltzmann = 8.617385e-11 * MeV/kelvin;

//
//
//
static const HepDouble STP_Temperature = 273.15*kelvin;
static const HepDouble STP_Pressure    = 1.*atmosphere;
static const HepDouble kGasThreshold   = 10.*mg/cm3;

//
//
//
static const HepDouble universe_mean_density = 1.e-25*g/cm3;

#endif /* HEP_PHYSICAL_CONSTANTS_H */
