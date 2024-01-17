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

// Parameters comes from JQMD 
// Niita et al., JAERI-Data/Code 99-042
//
// 230307 Skyrme-QMD parameters added by Y-H. Sato and A. Haga

#include "G4LightIonQMDParameters.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"

G4ThreadLocal G4LightIonQMDParameters* G4LightIonQMDParameters::parameters = NULL;

G4LightIonQMDParameters::G4LightIonQMDParameters()
{
   G4Pow* pow=G4Pow::GetInstance();

   wl = 1.26; // width of wave packet [fm]
   hbc = 0.19732857;      //   h-bar c in GeVfm

   //Pauli
   cpw = 1.0 / 2.0 / wl;

   cph = 2.0 * wl / (hbc*hbc);

   cpc = 4.0;

   epsx = -20.0 ;

    
// JQMD
/*
   rho0 = 0.168;     // satulation density
   G4double rpot = 1.0/3.0;

   G4double ebinm = -16.0; // bounding energy [MeV]
   G4double ebin = ebinm * 0.001;

   G4double pfer  = hbc * pow->A13 ( 3./2. *pi*pi * rho0 );

   G4double rmass = 0.938;

   G4double efer  = pfer*pfer / 2. / rmass;

   G4double t3 = 8. / 3. / rpot / pow->powA( rho0 , ( 1.+rpot ) ) * ( efer / 5. - ebin );

   G4double t0 = -16./15. * efer / rho0 - ( 1.+rpot ) * t3 * pow->powA( rho0 , rpot );
      

   G4double aaa = 3./4. * t0 * rho0;
   G4double bbb = 3./8. * t3 * ( 2.+rpot ) * pow->powA( rho0 , ( 1.+rpot ) );
   G4double esymm = 25 * 0.001; // symetric potential 25 [MeV] -> GeV

   gamm = rpot + 1.0;
*/
    
    
// Skyrme-QMD
// Ref. Y. Zhang and Z. Li, Elliptic flow and system size dependence of transition energies at intermediate energies, Phys.Rev. C74 (2006) 014602.
        
// ImQMD-SLy4
/*
   rho0 = 0.159546;
   G4double aaa = -297.82 * 0.001;
   G4double bbb = 219.21 * 0.001;
   gamm = 7.0/6;
   eta = 5.0/3;
   kappas = 0.08;
   g0 = 24.569/(2 * rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001;
   g0iso = 4.557/(rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001;
   gtau0 = 9.70/(pow->powA( rho0 , eta ) * pow->powA ( (4.0*pi*wl) , (1.5*eta) ))  * 0.001;
   G4double esymm = 32 * 0.001;
*/
        
// ImQMD-SkMstar
    
   //rho0 = 0.165;
   rho0 = 0.1603;     // satulation density
   G4double aaa = -318.0 * 0.001;
   G4double bbb = 249.5 * 0.001;
   gamm = 7.0/6;
   eta = 5.0/3;
   kappas = 0.08;
   g0 = 21.86/(2 * rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001;
   //g0iso = -5.485/(rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001; -> kappas
   gtau0 = 5.9357/(pow->powA( rho0 , eta ) * pow->powA ( (4.0*pi*wl) , (1.5*eta) ))  * 0.001;
   G4double esymm = 32 * 0.001;
    
// ImQMD-SIII
/*
   rho0 = 0.1452;     // satulation density
   G4double aaa = -122.921 * 0.001;
   G4double bbb = 55.343 * 0.001;
   gamm = 2;
   eta = 5.0/3;
   kappas = 0.08;
   g0 = 18.286/(2 * rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001;
   //g0iso = -5.485/(rho0 * pow->powA( 4 * pi * wl , 1.5 )) * 0.001; -> kappas
   gtau0 = 6.439/(pow->powA( rho0 , eta ) * pow->powA ( (4.0*pi*wl) , (1.5*eta) ))  * 0.001;
   G4double esymm = 28.17 * 0.001;
*/

// Local Potenials
   c0 = aaa / ( rho0 * pow->powA( 4 * pi * wl , 1.5 ) * 2.0 );

   c3 = bbb / ( pow->powA( rho0 , gamm ) * pow->powA ( (4.0*pi*wl) , (1.5*gamm) ) * ( gamm+1.0) );

   cs = esymm / ( rho0 * pow->powA( (4.0*pi*wl) , 1.5 ) * 2.0 );

   G4double ccoul = 0.001439767;
   cl = ccoul/2.0 * 1;  // Include Coulomb interaction
   //cl = ccoul/2.0 * 0;  // Not Include Coulomb interaction
 


// GroundStateNucleus
   cdp = 1.0 / pow->powA ( ( 4.0 * pi * wl ) , 1.5 );
   c0p = c0 * 2.0;
   c3p = c3 * ( gamm + 1.0 );
   csp = cs * 2.0;
   clp = cl * 2.0;
    
   g0p = g0 * 2.0; // Skyrme-QMD
   g0isop = g0iso * 2.0; // Skyrme-QMD
   gtau0p = gtau0 * ( eta + 1.0 ); // Skyrme-QMD

}



G4LightIonQMDParameters::~G4LightIonQMDParameters()
{ 
   ;
}

