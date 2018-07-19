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
// $Id: G4SIunits.hh 96706 2016-05-02 09:31:38Z gcosmo $
// 
// ----------------------------------------------------------------------
//
// Class description:
//
// This file is a modified version of SystemOfUnits.h
// It is provided for checking the overall 'units coherence' of the
// Geant4 kernel.
// -------
// Warning: if you use it, do not forget to recompile the whole Geant4 kernel
// ------- 
// The basic units are those of the International System:
//
//  	meter             
// 		second             
// 		kilogram      
// 		ampere         
// 		degree kelvin          
//      the amount of substance (mole)
//      luminous intensity      (candela)
// 		radian                  
//      steradian              
//
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
// 10.03.99   created
// 01.03.01   parsec
// 11.06.15   upgrate. Equivalent to SystemOfUnits.h
// 08.08.15   add decimeter, liter  (mma)      
// 12.01.16   added symbols for microsecond (us) and picosecond (ps) (mma)

#ifndef SI_SYSTEM_OF_UNITS_HH
#define SI_SYSTEM_OF_UNITS_HH

//
//
//
static constexpr double     pi  = 3.14159265358979323846;
static constexpr double  twopi  = 2*pi;
static constexpr double halfpi  = pi/2;
static constexpr double     pi2 = pi*pi;
// 
// Length [L]
//
static constexpr double meter  = 1.;                  
static constexpr double meter2 = meter*meter;
static constexpr double meter3 = meter*meter*meter;

static constexpr double millimeter  = 0.001*meter;                        
static constexpr double millimeter2 = millimeter*millimeter;
static constexpr double millimeter3 = millimeter*millimeter*millimeter;

static constexpr double centimeter  = 10.*millimeter;   
static constexpr double centimeter2 = centimeter*centimeter;
static constexpr double centimeter3 = centimeter*centimeter*centimeter;
  
static constexpr double kilometer = 1000.*meter;                   
static constexpr double kilometer2 = kilometer*kilometer;
static constexpr double kilometer3 = kilometer*kilometer*kilometer;

static constexpr double parsec = 3.0856775807e+16*meter;

static constexpr double micrometer = 1.e-6 *meter;             
static constexpr double  nanometer = 1.e-9 *meter;
static constexpr double  angstrom  = 1.e-10*meter;
static constexpr double  fermi     = 1.e-15*meter;

static constexpr double      barn = 1.e-28*meter2;
static constexpr double millibarn = 1.e-3 *barn;
static constexpr double microbarn = 1.e-6 *barn;
static constexpr double  nanobarn = 1.e-9 *barn;
static constexpr double  picobarn = 1.e-12*barn;

// symbols
static constexpr double nm  = nanometer;                        
static constexpr double um  = micrometer;                        

static constexpr double mm  = millimeter;                        
static constexpr double mm2 = millimeter2;
static constexpr double mm3 = millimeter3;

static constexpr double cm  = centimeter;   
static constexpr double cm2 = centimeter2;
static constexpr double cm3 = centimeter3;

static constexpr double liter = 1.e+3*cm3;
static constexpr double  L = liter;
static constexpr double dL = 1.e-1*liter;
static constexpr double cL = 1.e-2*liter;
static constexpr double mL = 1.e-3*liter;
           
static constexpr double m  = meter;                  
static constexpr double m2 = meter2;
static constexpr double m3 = meter3;

static constexpr double km  = kilometer;                   
static constexpr double km2 = kilometer2;
static constexpr double km3 = kilometer3;

static constexpr double pc = parsec;

//
// Angle
//
static constexpr double radian      = 1.;                  
static constexpr double milliradian = 1.e-3*radian;
static constexpr double degree = (pi/180.0)*radian;

static constexpr double   steradian = 1.;
	
// symbols
static constexpr double rad  = radian;	
static constexpr double mrad = milliradian;
static constexpr double sr   = steradian;
static constexpr double deg  = degree;

//
// Time [T]
//
static constexpr double second      = 1.;
static constexpr double nanosecond  = 1.e-9 *second;
static constexpr double millisecond = 1.e-3 *second;
static constexpr double microsecond = 1.e-6 *second;
static constexpr double  picosecond = 1.e-12*second;

static constexpr double hertz = 1./second;
static constexpr double kilohertz = 1.e+3*hertz;
static constexpr double megahertz = 1.e+6*hertz;

// symbols
static constexpr double ns = nanosecond;			
static constexpr double  s = second;
static constexpr double ms = millisecond;
static constexpr double us = microsecond;
static constexpr double ps = picosecond;

//
// Mass [E][T^2][L^-2]
//
static constexpr double  kilogram = 1.;   
static constexpr double      gram = 1.e-3*kilogram;
static constexpr double milligram = 1.e-3*gram;

// symbols
static constexpr double  kg = kilogram;
static constexpr double   g = gram;
static constexpr double  mg = milligram;

//
// Electric current [Q][T^-1]
//
static constexpr double      ampere = 1.;
static constexpr double milliampere = 1.e-3*ampere;
static constexpr double microampere = 1.e-6*ampere;
static constexpr double  nanoampere = 1.e-9*ampere;

//
// Electric charge [Q]
//
static constexpr double coulomb = ampere*second;
static constexpr double e_SI  = 1.602176487e-19;	// positron charge in coulomb
static constexpr double eplus = e_SI*coulomb ;		// positron charge

//
// Energy [E]
//
static constexpr double joule = kg*m*m/(s*s);

static constexpr double     electronvolt = e_SI*joule;
static constexpr double kiloelectronvolt = 1.e+3*electronvolt;
static constexpr double megaelectronvolt = 1.e+6*electronvolt; 
static constexpr double gigaelectronvolt = 1.e+9*electronvolt;
static constexpr double teraelectronvolt = 1.e+12*electronvolt;
static constexpr double petaelectronvolt = 1.e+15*electronvolt;

// symbols
static constexpr double MeV = megaelectronvolt;
static constexpr double  eV = electronvolt;
static constexpr double keV = kiloelectronvolt;
static constexpr double GeV = gigaelectronvolt;
static constexpr double TeV = teraelectronvolt;
static constexpr double PeV = petaelectronvolt;

//
// Power [E][T^-1]
//
static constexpr double watt = joule/second;	// watt = 6.24150 e+3 * MeV/ns

//
// Force [E][L^-1]
//
static constexpr double newton = joule/meter;	// newton = 6.24150 e+9 * MeV/mm

//
// Pressure [E][L^-3]
//
#define pascal hep_pascal                          // a trick to avoid warnings 
static constexpr double hep_pascal = newton/m2;	   // pascal = 6.24150 e+3 * MeV/mm3
static constexpr double bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
static constexpr double atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

//
// Electric potential [E][Q^-1]
//
static constexpr double megavolt = megaelectronvolt/eplus;
static constexpr double kilovolt = 1.e-3*megavolt;
static constexpr double     volt = 1.e-6*megavolt;

//
// Electric resistance [E][T][Q^-2]
//
static constexpr double ohm = volt/ampere;	// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

//
// Electric capacitance [Q^2][E^-1]
//
static constexpr double farad = coulomb/volt;	// farad = 6.24150e+24 * eplus/Megavolt
static constexpr double millifarad = 1.e-3*farad;
static constexpr double microfarad = 1.e-6*farad;
static constexpr double  nanofarad = 1.e-9*farad;
static constexpr double  picofarad = 1.e-12*farad;

//
// Magnetic Flux [T][E][Q^-1]
//
static constexpr double weber = volt*second;	// weber = 1000*megavolt*ns

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static constexpr double tesla     = volt*second/meter2;	// tesla =0.001*megavolt*ns/mm2

static constexpr double gauss     = 1.e-4*tesla;
static constexpr double kilogauss = 1.e-1*tesla;

//
// Inductance [T^2][E][Q^-2]
//
static constexpr double henry = weber/ampere;	// henry = 1.60217e-7*MeV*(ns/eplus)**2

//
// Temperature
//
static constexpr double kelvin = 1.;

//
// Amount of substance
//
static constexpr double mole = 1.;

//
// Activity [T^-1]
//
static constexpr double becquerel = 1./second ;
static constexpr double curie = 3.7e+10 * becquerel;
static constexpr double kilobecquerel = 1.e+3*becquerel;
static constexpr double megabecquerel = 1.e+6*becquerel;
static constexpr double gigabecquerel = 1.e+9*becquerel;
static constexpr double millicurie = 1.e-3*curie;
static constexpr double microcurie = 1.e-6*curie;
static constexpr double Bq = becquerel;
static constexpr double kBq = kilobecquerel;
static constexpr double MBq = megabecquerel;
static constexpr double GBq = gigabecquerel;
static constexpr double Ci = curie;
static constexpr double mCi = millicurie;
static constexpr double uCi = microcurie;

//
// Absorbed dose [L^2][T^-2]
//
static constexpr double      gray = joule/kilogram;
static constexpr double  kilogray = 1.e+3*gray;
static constexpr double milligray = 1.e-3*gray;
static constexpr double microgray = 1.e-6*gray;

//
// Luminous intensity [I]
//
static constexpr double candela = 1.;

//
// Luminous flux [I]
//
static constexpr double lumen = candela*steradian;

//
// Illuminance [I][L^-2]
//
static constexpr double lux = lumen/meter2;

//
// Miscellaneous
//
static constexpr double perCent     = 0.01 ;
static constexpr double perThousand = 0.001;
static constexpr double perMillion  = 0.000001;


#endif /* SI_SYSTEM_OF_UNITS_HH */
