// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: SystemOfUnits.h,v 1.1 1999-01-07 16:09:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -*- C++ -*-
// 
// ----------------------------------------------------------------------
// HEP coherent system of Units
//
// This file has been provided by Geant4 (simulation toolkit for HEP).
//
// The basic units are :
//  		millimeter              (millimeter)
// 		nanosecond              (nanosecond)
// 		Mega electron Volt      (MeV)
// 		positron charge         (eplus)
// 		degree Kelvin           (kelvin)
//              the amount of substance (mole)
//              luminous intensity      (candela)
// 		radian                  (radian)
//              steradian               (steradian)
//
// Below is a non exhaustive List of derived and pratical units
// (i.e. mostly the SI units).
// You can add your own units.
//
// The SI numerical value of the positron charge is defined here,
// as it is needed for conversion factor : positron charge = e_SI (coulomb)
//
// The others physical constants are defined in the header file :
//			PhysicalConstants.h
//
// Authors: M.Maire, S.Giani
//
// History:
//
// 06.02.96   Created.
// 28.03.96   Added miscellaneous constants.
// 05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
// 20.05.98   names: meter, second, gram, radian, degree. (from Blasiuk (STAR))
//            Added luminous units.
// 05.08.98   angstrom,picobarn,microsecond,picosecond,petaelectronvolt    

#ifndef HEP_SYSTEM_OF_UNITS_H
#define HEP_SYSTEM_OF_UNITS_H

#include "CLHEP/config/CLHEP.h"

// 
// Length [L]
//
static const HepDouble millimeter  = 1.;                        
static const HepDouble millimeter2 = millimeter*millimeter;
static const HepDouble millimeter3 = millimeter*millimeter*millimeter;

static const HepDouble centimeter  = 10.*millimeter;   
static const HepDouble centimeter2 = centimeter*centimeter;
static const HepDouble centimeter3 = centimeter*centimeter*centimeter;

static const HepDouble meter  = 1000.*millimeter;                  
static const HepDouble meter2 = meter*meter;
static const HepDouble meter3 = meter*meter*meter;

static const HepDouble kilometer = 1000.*meter;                   
static const HepDouble kilometer2 = kilometer*kilometer;
static const HepDouble kilometer3 = kilometer*kilometer*kilometer;

static const HepDouble micrometer = 1.e-6 *meter;             
static const HepDouble  nanometer = 1.e-9 *meter;
static const HepDouble  angstrom  = 1.e-10*meter;
static const HepDouble  fermi     = 1.e-15*meter;

static const HepDouble      barn = 1.e-28*meter2;
static const HepDouble millibarn = 1.e-3 *barn;
static const HepDouble microbarn = 1.e-6 *barn;
static const HepDouble  nanobarn = 1.e-9 *barn;
static const HepDouble  picobarn = 1.e-12*barn;

// symbols
static const HepDouble mm  = millimeter;                        
static const HepDouble mm2 = millimeter2;
static const HepDouble mm3 = millimeter3;

static const HepDouble cm  = centimeter;   
static const HepDouble cm2 = centimeter2;
static const HepDouble cm3 = centimeter3;

static const HepDouble m  = meter;                  
static const HepDouble m2 = meter2;
static const HepDouble m3 = meter3;

static const HepDouble km  = kilometer;                   
static const HepDouble km2 = kilometer2;
static const HepDouble km3 = kilometer3;

//
// Angle
//
static const HepDouble radian      = 1.;                  
static const HepDouble milliradian = 1.e-3*radian;
static const HepDouble degree = (3.14159265358979323846/180.0)*radian;

static const HepDouble   steradian = 1.;
	
// symbols
static const HepDouble rad  = radian;	
static const HepDouble mrad = milliradian;
static const HepDouble sr   = steradian;
static const HepDouble deg  = degree;

//
// Time [T]
//
static const HepDouble nanosecond  = 1.;
static const HepDouble second      = 1.e+9 *nanosecond;
static const HepDouble millisecond = 1.e-3 *second;
static const HepDouble microsecond = 1.e-6 *second;
static const HepDouble  picosecond = 1.e-12*second;

static const HepDouble hertz = 1./second;
static const HepDouble kilohertz = 1.e+3*hertz;
static const HepDouble megahertz = 1.e+6*hertz;

// symbols
static const HepDouble ns = nanosecond;			
static const HepDouble  s = second;
static const HepDouble ms = millisecond;

//
// Electric charge [Q]
//
static const HepDouble eplus = 1. ;		// positron charge
static const HepDouble e_SI  = 1.60217733e-19;	// positron charge in coulomb
static const HepDouble coulomb = eplus/e_SI;	// coulomb = 6.24150 e+18 * eplus

//
// Energy [E]
//
static const HepDouble megaelectronvolt = 1. ;
static const HepDouble     electronvolt = 1.e-6*megaelectronvolt;
static const HepDouble kiloelectronvolt = 1.e-3*megaelectronvolt;
static const HepDouble gigaelectronvolt = 1.e+3*megaelectronvolt;
static const HepDouble teraelectronvolt = 1.e+6*megaelectronvolt;
static const HepDouble petaelectronvolt = 1.e+9*megaelectronvolt;

static const HepDouble joule = electronvolt/e_SI;	// joule = 6.24150 e+12 * MeV

// symbols
static const HepDouble MeV = megaelectronvolt;
static const HepDouble  eV = electronvolt;
static const HepDouble keV = kiloelectronvolt;
static const HepDouble GeV = gigaelectronvolt;
static const HepDouble TeV = teraelectronvolt;
static const HepDouble PeV = petaelectronvolt;

//
// Mass [E][T^2][L^-2]
//
static const HepDouble  kilogram = joule*second*second/(meter*meter);   
static const HepDouble      gram = 1.e-3*kilogram;
static const HepDouble milligram = 1.e-3*gram;

// symbols
static const HepDouble  kg = kilogram;
static const HepDouble   g = gram;
static const HepDouble  mg = milligram;

//
// Power [E][T^-1]
//
static const HepDouble watt = joule/second;	// watt = 6.24150 e+3 * MeV/ns

//
// Force [E][L^-1]
//
static const HepDouble newton = joule/meter;	// newton = 6.24150 e+9 * MeV/mm

//
// Pressure [E][L^-3]
//
#define pascal hep_pascal                          // a trick to avoid warnings 
static const HepDouble hep_pascal = newton/m2;	   // pascal = 6.24150 e+3 * MeV/mm3
static const HepDouble bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
static const HepDouble atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

//
// Electric current [Q][T^-1]
//
static const HepDouble      ampere = coulomb/second; // ampere = 6.24150 e+9 * eplus/ns
static const HepDouble milliampere = 1.e-3*ampere;
static const HepDouble microampere = 1.e-6*ampere;
static const HepDouble  nanoampere = 1.e-9*ampere;

//
// Electric potential [E][Q^-1]
//
static const HepDouble megavolt = megaelectronvolt/eplus;
static const HepDouble kilovolt = 1.e-3*megavolt;
static const HepDouble     volt = 1.e-6*megavolt;

//
// Electric resistance [E][T][Q^-2]
//
static const HepDouble ohm = volt/ampere;	// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

//
// Electric capacitance [Q^2][E^-1]
//
static const HepDouble farad = coulomb/volt;	// farad = 6.24150e+24 * eplus/Megavolt
static const HepDouble millifarad = 1.e-3*farad;
static const HepDouble microfarad = 1.e-6*farad;
static const HepDouble  nanofarad = 1.e-9*farad;
static const HepDouble  picofarad = 1.e-12*farad;

//
// Magnetic Flux [T][E][Q^-1]
//
static const HepDouble weber = volt*second;	// weber = 1000*megavolt*ns

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static const HepDouble tesla     = volt*second/meter2;	// tesla =0.001*megavolt*ns/mm2

static const HepDouble gauss     = 1.e-4*tesla;
static const HepDouble kilogauss = 1.e-1*tesla;

//
// Inductance [T^2][E][Q^-2]
//
static const HepDouble henry = weber/ampere;	// henry = 1.60217e-7*MeV*(ns/eplus)**2

//
// Temperature
//
static const HepDouble kelvin = 1.;

//
// Amount of substance
//
static const HepDouble mole = 1.;

//
// Activity [T^-1]
//
static const HepDouble becquerel = 1./second ;
static const HepDouble curie = 3.7e+10 * becquerel;

//
// Absorbed dose [L^2][T^-2]
//
static const HepDouble gray = joule/kilogram ;

//
// Luminous intensity [I]
//
static const HepDouble candela = 1.;

//
// Luminous flux [I]
//
static const HepDouble lumen = candela*steradian;

//
// Illuminance [I][L^-2]
//
static const HepDouble lux = lumen/meter2;

//
// Miscellaneous
//
static const HepDouble perCent     = 0.01 ;
static const HepDouble perThousand = 0.001;
static const HepDouble perMillion  = 0.000001;

#endif /* HEP_SYSTEM_OF_UNITS_H */
