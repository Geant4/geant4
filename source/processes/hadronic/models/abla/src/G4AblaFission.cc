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
// ABLAXX statistical de-excitation model
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4AblaFission.hh"
#include <time.h>
#include <cmath>

G4AblaFission::G4AblaFission()
{
}

G4AblaFission::~G4AblaFission()
{
}

void G4AblaFission::doFission(G4double &A, G4double &Z, G4double &E,
			      G4double &A1, G4double &Z1, G4double &E1, G4double &K1,
			      G4double &A2, G4double &Z2, G4double &E2, G4double &K2)
{
  fissionDistri(A,Z,E,A1,Z1,E1,K1,A2,Z2,E2,K2);
}

void G4AblaFission::even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out)     
{
  // Procedure to calculate I_OUT from R_IN in a way that
  // on the average a flat distribution in R_IN results in a
  // fluctuating distribution in I_OUT with an even-odd effect as
  // given by R_EVEN_ODD

  //     /* ------------------------------------------------------------ */
  //     /* EXAMPLES :                                                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD = 0 :                                       */
  //     /*           CEIL(R_IN)  ----                                   */
  //     /*                                                              */
  //     /*              R_IN ->                                         */
  //     /*            (somewhere in between CEIL(R_IN) and FLOOR(R_IN)) */                                            */
  //     /*                                                              */
  //     /*           FLOOR(R_IN) ----       --> I_OUT                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD > 0 :                                       */
  //     /*      The interval for the above treatment is                 */
  //     /*         larger for FLOOR(R_IN) = even and                    */
  //     /*         smaller for FLOOR(R_IN) = odd                        */
  //     /*    For R_EVEN_ODD < 0 : just opposite treatment              */
  //     /* ------------------------------------------------------------ */

  //     /* ------------------------------------------------------------ */
  //     /* On input:   R_ORIGIN    nuclear charge (real number)         */
  //     /*             R_EVEN_ODD  requested even-odd effect            */
  //     /* Intermediate quantity: R_IN = R_ORIGIN + 0.5                 */
  //     /* On output:  I_OUT       nuclear charge (integer)             */
  //     /* ------------------------------------------------------------ */

  //      G4double R_ORIGIN,R_IN,R_EVEN_ODD,R_REST,R_HELP;
  G4double r_in = 0.0, r_rest = 0.0, r_help = 0.0;
  G4double r_floor = 0.0;
  G4double r_middle = 0.0;
  //      G4int I_OUT,N_FLOOR;
  G4int n_floor = 0;

  r_in = r_origin + 0.5;
  r_floor = (float)((int)(r_in));
  if (r_even_odd < 0.001) {
    i_out = (int)(r_floor);
  } 
  else {
    r_rest = r_in - r_floor;
    r_middle = r_floor + 0.5;
    n_floor = (int)(r_floor);
    if (n_floor%2 == 0) {
      // even before modif.
      r_help = r_middle + (r_rest - 0.5) * (1.0 - r_even_odd);
    } 
    else {
      // odd before modification
      r_help = r_middle + (r_rest - 0.5) * (1.0 + r_even_odd);
    }
    i_out = (int)(r_help);
  }
}

G4double G4AblaFission::umass(G4double z,G4double n,G4double beta)
{
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
  // pure liquid drop, without pairing and shell effects

  // On input:    Z     nuclear charge of nucleus
  //              N     number of neutrons in nucleus
  //              beta  deformation of nucleus
  // On output:   binding energy of nucleus

  G4double a = 0.0, fumass = 0.0;
  G4double alpha = 0.0;
  G4double xcom = 0.0, xvs = 0.0, xe = 0.0;
  const G4double pi = 3.1416;

  a = n + z;
  alpha = ( std::sqrt(5.0/(4.0*pi)) ) * beta;
  
  xcom = 1.0 - 1.7826 * ((a - 2.0*z)/a)*((a - 2.0*z)/a);
  // factor for asymmetry dependence of surface and volume term
  xvs = - xcom * ( 15.4941 * a - 
		   17.9439 * std::pow(a,2.0/3.0) * (1.0+0.4*alpha*alpha) );
  // sum of volume and surface energy
  xe = z*z * (0.7053/(std::pow(a,1.0/3.0)) * (1.0-0.2*alpha*alpha) - 1.1529/a);
  fumass = xvs + xe;
  
  return fumass;
}

G4double G4AblaFission::ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d)
{
  // Coulomb potential between two nuclei
  // surfaces are in a distance of d
  // in a tip to tip configuration

  // approximate formulation
  // On input: Z1      nuclear charge of first nucleus
  //           N1      number of neutrons in first nucleus
  //           beta1   deformation of first nucleus
  //           Z2      nuclear charge of second nucleus
  //           N2      number of neutrons in second nucleus
  //           beta2   deformation of second nucleus
  //           d       distance of surfaces of the nuclei

  //      G4double Z1,N1,beta1,Z2,N2,beta2,d,ecoul;
  G4double fecoul = 0;
  G4double dtot = 0;
  const G4double r0 = 1.16;

  dtot = r0 * ( std::pow((z1+n1),1.0/3.0) * (1.0+0.6666*beta1)
		+ std::pow((z2+n2),1.0/3.0) * (1.0+0.6666*beta2) ) + d;
  fecoul = z1 * z2 * 1.44 / dtot;

  return fecoul;
}

void G4AblaFission::fissionDistri(G4double &a,G4double &z,G4double &e,
				  G4double &a1,G4double &z1,G4double &e1,G4double &v1,
				  G4double &a2,G4double &z2,G4double &e2,G4double &v2)
{
  //  // G4cout <<"Fission: a = " << a << " z = " << z << " e = " << e << G4endl;
  //  On input: A, Z, E (mass, atomic number and exc. energy of compound nucleus
  //                     before fission)
  //  On output: Ai, Zi, Ei (mass, atomic number and exc. energy of fragment 1 and 2
  //                     after fission)

  //  Additionally calculated but not put in the parameter list:
  //  Kinetic energy of prefragments EkinR1, EkinR2

  //  Translation of SIMFIS18.PLI (KHS, 2.1.2001)

  // This program calculates isotopic distributions of fission fragments
  // with a semiempirical model                      
  // Copy from SIMFIS3, KHS, 8. February 1995        
  // Modifications made by Jose Benlliure and KHS in August 1996
  // Energy counted from lowest barrier (J. Benlliure, KHS 1997)
  // Some bugs corrected (J. Benlliure, KHS 1997)         
  // Version used for thesis S. Steinhaueser (August 1997)
  // (Curvature of LD potential increased by factor of 2!)

  // Weiter veraendert mit der Absicht, eine Version zu erhalten, die
  // derjenigen entspricht, die von J. Benlliure et al.
  // in Nucl. Phys. A 628 (1998) 458 verwendet wurde,
  // allerdings ohne volle Neutronenabdampfung.

  // The excitation energy was calculate now for each fission channel
  // separately. The dissipation from saddle to scission was taken from
  // systematics, the deformation energy at scission considers the shell
  // effects in a simplified way, and the fluctuation is included. 
  // KHS, April 1999 

  // The width in N/Z was carefully adapted to values given by Lang et al.

  // The width and eventually a shift in N/Z (polarization) follows the
  // following rules:                                        

  // The line N/Z following UCD has an angle of std::atan(Zcn/Ncn)
  // to the horizontal axis on a chart of nuclides.
  // (For 238U the angle is 32.2 deg.)

  // The following relations hold: (from Armbruster)
  // 
  // sigma(N) (A=const) = sigma(Z) (A=const)
  // sigma(A) (N=const) = sigma(Z) (N=const)
  // sigma(A) (Z=const) = sigma(N) (Z=const)
  // 
  // From this we get:
  // sigma(Z) (N=const) * N = sigma(N) (Z=const) * Z
  // sigma(A) (Z=const) = sigma(Z) (A=const) * A/Z
  // sigma(N) (Z=const) = sigma(Z) (A=const) * A/Z
  // Z*sigma(N) (Z=const) = N*sigma(Z) (N=const) = A*sigma(Z) (A=const)

  // Excitation energy now calculated above the lowest potential point
  // Inclusion of a distribution of excitation energies             

  // Several modifications, starting from SIMFIS12: KHS November 2000 
  // This version seems to work quite well for 238U.                  
  // The transition from symmetric to asymmetric fission around 226Th 
  // is reasonably well reproduced, although St. I is too strong and St. II
  // is too weak. St. I and St. II are also weakly seen for 208Pb.

  // Extensions for an event generator of fission events (21.11.2000,KHS)

  // Defalt parameters (IPARS) rather carefully adjusted to
  // pre-neutron mass distributions of Vives et al. (238U + n)
  // Die Parameter Fgamma1 und Fgamma2 sind kleiner als die resultierenden
  // Breiten der Massenverteilungen!!!  
  // Fgamma1 und Fgamma2 wurden angepa�, so da�
  // Sigma-A(ST-I) = 3.3, Sigma-A(St-II) = 5.8 (nach Vives) 

  // Parameters of the model carefully adjusted by KHS (2.2.2001) to
  // 238U + 208Pb, 1000 A MeV, Timo Enqvist et al.         

  G4double     n = 0.0;
  G4double     aheavy1 = 0.0, aheavy2 = 0.0;
  //  G4double     eheavy1 = 0.0, elight1 = 0.0, eheavy2 = 0.0, elight2 = 0.0;
  G4double     zheavy1_shell = 0.0, zheavy2_shell = 0.0;
  G4double     masscurv = 0.0;
  G4double     sasymm1 = 0.0, sasymm2 = 0.0, ssymm = 0.0, ysum = 0.0;
  G4double     ssymm_mode1 = 0.0, ssymm_mode2 = 0.0;
  // Curvature at saddle, modified by ld-potential
  G4double     wzasymm1_saddle, wzasymm2_saddle, wzsymm_saddle  = 0.0;
  G4double     wzasymm1_scission = 0.0, wzasymm2_scission = 0.0, wzsymm_scission = 0.0;
  G4double     wzasymm1 = 0.0, wzasymm2 = 0.0, wzsymm = 0.0;
  G4int  imode = 0;
  G4double     rmode = 0.0;
  G4double     z1mean = 0.0, z1width = 0.0;
  //      G4double     Z1,Z2,N1R,N2R,A1R,A2R,N1,N2,A1,A2;
  G4double     n1r = 0.0, n2r = 0.0;

  G4double     zsymm = 0.0, nsymm = 0.0, asymm = 0.0;
  G4double     n1mean = 0.0, n1width = 0.0;
  // effective shell effect at lowest barrier
  // Excitation energy with respect to ld barrier
  G4double     re1 = 0.0, re2 = 0.0, re3 = 0.0;
  G4double     n1ucd = 0.0, n2ucd = 0.0;

  // shift of most probable neutron number for given Z,
  // according to polarization
  G4int  i_help = 0;

  //   /* Parameters of the semiempirical fission model */
  G4double a_levdens = 0.0;
  //           /* level-density parameter */
  G4double a_levdens_heavy1 = 0.0, a_levdens_heavy2 = 0.0;
  const G4double r_null = 1.16;
  //          /* radius parameter */
  G4double epsilon_1_saddle = 0.0, epsilon0_1_saddle = 0.0;
  G4double epsilon_2_saddle = 0.0, epsilon0_2_saddle = 0.0, epsilon_symm_saddle = 0.0;
  G4double epsilon_1_scission = 0.0, epsilon0_1_scission = 0.0;
  G4double epsilon_2_scission = 0.0, epsilon0_2_scission = 0.0;
  G4double epsilon_symm_scission = 0.0;
  //                                   /* modified energy */
  G4double e_eff1_saddle = 0.0, e_eff2_saddle = 0.0;
  G4int icz = 0, k = 0;

  G4int i_inter = 0;
  G4double ne_min = 0;
  G4double ed1_low = 0.0, ed2_low = 0.0, ed1_high = 0.0, ed2_high = 0.0, ed1 = 0.0, ed2 = 0.0;
  G4double atot = 0.0;

  //   Input parameters:
  //OMMENT(Nuclear charge number);
  //      G4double Z;
  //OMMENT(Nuclear mass number);
  //      G4double A;
  //OMMENT(Excitation energy above fission barrier);
  //      G4double E;

  //   Model parameters:
  //OMMENT(position of heavy peak valley 1);
  const G4double nheavy1 = 82.0;
  const G4double nheavy2 = 89.0;
  //OMMENT(position of heavy peak valley 2);
  const G4double e_crit = 5; // Critical pairing energy :::PK
  //OMMENT(Shell effect for valley 1);
  //        Parameter (Delta_U2_shell = -3.2)
  //OMMENT(I: used shell effect);
  G4double delta_u1 = 0.0;
  //omment(I: used shell effect);
  G4double delta_u2 = 0.0;
  const G4double delta_u1_shell = -2.5;
  //        Parameter (Delta_U1_shell = -2)
  //OMMENT(Shell effect for valley 2);
  const G4double delta_u2_shell = -5.5;
  const G4double el = 30.0;
  //OMMENT(Curvature of asymmetric valley 1);
  const G4double cz_asymm1_shell = 0.7;
  //OMMENT(Curvature of asymmetric valley 2);
  const G4double cz_asymm2_shell = 0.08;
  //OMMENT(Factor for width of distr. valley 1);
  const G4double fwidth_asymm1 = 1.2;
  //OMMENT(Factor for width of distr. valley 2);
  const G4double fwidth_asymm2 = 1.0;
  //       Parameter (CZ_asymm2_scission = 0.12)
  //OMMENT(Factor to gamma_heavy1);
  const G4double fgamma1 = 1.0;
  //OMMENT(I: fading of shells (general));
  G4double gamma = 0.0;
  //OMMENT(I: fading of shell 1);
  G4double gamma_heavy1 = 0.0;
  //OMMENT(I: fading of shell 2);
  G4double gamma_heavy2 = 0.0;
  //OMMENT(Zero-point energy at saddle);
  const G4double e_zero_point = 0.5;
  G4int i_eva = 0; // Calculate  A = 1 or Aprime = 0 :::PK
  //OMMENT(I: friction from saddle to scission);
  G4double e_saddle_scission = 10.0;
  //OMMENT(Friction factor);
  const G4double friction_factor = 1.0;
  //OMMENT(I: Internal counter for different modes); INIT(0,0,0)
  //      Integer*4 I_MODE(3)
  //OMMENT(I: Yield of symmetric mode);
  G4double ysymm = 0.0;
  //OMMENT(I: Yield of asymmetric mode 1);
  G4double yasymm1 = 0.0;
  //OMMENT(I: Yield of asymmetric mode 2);
  G4double yasymm2 = 0.0;
  //OMMENT(I: Effective position of valley 1);
  G4double nheavy1_eff = 0.0;
  //OMMENT(I: position of heavy peak valley 1);
  G4double zheavy1 = 0.0;
  //omment(I: Effective position of valley 2);
  G4double nheavy2_eff = 0.0;
  //OMMENT(I: position of heavy peak valley 2);
  G4double zheavy2 = 0.0;
  //omment(I: Excitation energy above saddle 1);
  G4double eexc1_saddle = 0.0;
  //omment(I: Excitation energy above saddle 2);
  G4double eexc2_saddle = 0.0;
  //omment(I: Excitation energy above lowest saddle);
  G4double eexc_max = 0.0;
  //OMMENT(I: Even-odd effect in Z);
  G4double r_e_o = 0.0;
  G4double r_e_o_max = 0.0;
  G4double e_pair = 0.0;
  //OMMENT(I: Curveture of symmetric valley);
  G4double cz_symm = 0.0;
  //OMMENT(I: Curvature of mass distribution for fixed Z);
  G4double cn = 0.0;
  //OMMENT(=1: test output, =0: no test output);
  const G4int itest = 0;
  //      G4double UMASS, ECOUL, reps1, reps2, rn1_pol;
  G4double reps1 = 0.0, reps2 = 0.0, rn1_pol = 0.0;
  //      Float_t HAZ,GAUSSHAZ;
  //G4int kkk = 0;
  G4int kkk = 10;
  G4double Bsym = 0.0;
  G4double Basym_1 = 0.0;
  G4double Basym_2 = 0.0;
  //     I_MODE = 0;

  //     /* average Z of asymmetric and symmetric components: */
  n = a - z;  /* neutron number of the fissioning nucleus */

  k = 0;
  icz = 0;
  if ( (std::pow(z,2)/a < 25.0) || (n < nheavy2) || (e > 500.0) ) {
    icz = -1;
    //          GOTO 1002;
    goto milledeux;
  }

  //    /* Polarisation assumed for standard I and standard II:
  //      Z - Zucd = cpol (for A = const);
  //      from this we get (see Armbruster)
  //      Z - Zucd =  Acn/Ncn * cpol (for N = const)                        */

  //  zheavy1_shell = ((nheavy1/n) * z) - ((a/n) * cpol1); // Simfis18  PK:::
  zheavy1_shell = ((nheavy1/n) * z) - 0.8;                 // Simfis3  PK:::
  //zheavy2_shell = ((nheavy2/n) * z) - ((a/n) * cpol2);   // Simfis18  PK:::
  zheavy2_shell = ((nheavy2/n) * z) - 0.8;                 // Simfis3  PK:::

  //  p(zheavy1_shell, zheavy2_shell);

  //  e_saddle_scission = 
  //    (-24.0 + 0.02227 * (std::pow(z,2))/(std::pow(a,0.33333)) ) * friction_factor;
  e_saddle_scission = (3.535 * std::pow(z,2)/a - 121.1) * friction_factor;  // PK:::

  //      /* Energy dissipated from saddle to scission                        */
  //      /* F. Rejmund et al., Nucl. Phys. A 678 (2000) 215, fig. 4 b        */
  //      E_saddle_scission = DMAX1(0.,E_saddle_scission);
  // Heavy Ion Induced Reactions, Schroeder W. ed., Harwood, 1986, 101
  if (e_saddle_scission < 0.0) {
    e_saddle_scission = 0.0;
  }

  //     /* Semiempirical fission model: */

  //    /* Fit to experimental result on curvature of potential at saddle */
  //           /* reference:                                              */
  //    /* IF Z**2/A < 33.15E0 THEN
  //       MassCurv = 30.5438538E0 - 4.00212049E0 * Z**2/A
  //                               + 0.11983384E0 * Z**4 / (A**2) ;
  //     ELSE
  //       MassCurv = 10.E0 ** (7.16993332E0 - 0.26602401E0 * Z**2/A
  //                               + 0.00283802E0 * Z**4 / (A**2)) ;  */
  //  /* New parametrization of T. Enqvist according to Mulgin et al. 1998 NPA 640(1998) 375 */
  if ( (std::pow(z,2))/a < 33.9186) {
    masscurv =  std::pow( 10.0,(-1.093364 + 0.082933 * (std::pow(z,2)/a)
				- 0.0002602 * (std::pow(z,4)/std::pow(a,2))) );
  } else {
    masscurv = std::pow( 10.0,(3.053536 - 0.056477 * (std::pow(z,2)/a)
			       + 0.0002454 * (std::pow(z,4)/std::pow(a,2))) );
  }

  cz_symm = (8.0/std::pow(z,2)) * masscurv;

  if(itest == 1) {
    //    // G4cout << "cz_symmetry= " << cz_symm << G4endl;
  }

  icz = 0;
  if (cz_symm < 0) {
    icz = -1;
    //          GOTO 1002;
    goto milledeux;
  }

  //  /* proton number in symmetric fission (centre) */
  zsymm  = z/2.0;
  nsymm  = n/2.0;
  asymm = nsymm + zsymm;

  zheavy1 = (cz_symm*zsymm + cz_asymm1_shell*zheavy1_shell)/(cz_symm + cz_asymm1_shell);
  zheavy2 = (cz_symm*zsymm + cz_asymm2_shell*zheavy2_shell)/(cz_symm + cz_asymm2_shell);

  //            /* position of valley due to influence of liquid-drop potential */
  nheavy1_eff = (zheavy1 + 0.8)*(n/z);
  nheavy2_eff = (zheavy2 + 0.8)*(n/z);
  aheavy1 = nheavy1_eff + zheavy1;
  aheavy2 = nheavy2_eff + zheavy2;
  // Eheavy1 = E * Aheavy1 / A
  // Eheavy2 = E * Aheavy2 / A
  // Elight1 = E * Alight1 / A
  // Elight2 = E * Alight2 / A
  a_levdens = a / 8.0;
  a_levdens_heavy1 = aheavy1 / 8.0;
  a_levdens_heavy2 = aheavy2 / 8.0;
  gamma = a_levdens / (0.4 * (std::pow(a,1.3333)) );
  gamma_heavy1 = ( a_levdens_heavy1 / (0.4 * (std::pow(aheavy1,1.3333)) ) ) * fgamma1;
  gamma_heavy2 = a_levdens_heavy2 / (0.4 * (std::pow(aheavy2,1.3333)) );

  // Up to here: Ok! Checked CS 10/10/05      	   

  cn = umass(zsymm,(nsymm+1.),0.0) + umass(zsymm,(nsymm-1.),0.0)
    + 1.44 * (std::pow(zsymm,2))/
    ( (std::pow(r_null,2)) * 
      ( std::pow((asymm+1.0),1.0/3.0) + std::pow((asymm-1.0),1.0/3.0) ) *
      ( std::pow((asymm+1.0),1.0/3.0) + std::pow((asymm-1.0),1.0/3.0) ) )
    - 2.0 * umass(zsymm,nsymm,0.0)
    - 1.44 * (std::pow(zsymm,2))/
    ( ( 2.0 * r_null * (std::pow(asymm,1.0/3.0)) ) * 
      ( 2.0 * r_null * (std::pow(asymm,1.0/3.0)) ) );

  // /* shell effect in valley of mode 1 */
  delta_u1 = delta_u1_shell + (std::pow((zheavy1_shell-zheavy1),2))*cz_asymm1_shell;
  // /* shell effect in valley of mode 2 */
  delta_u2 = delta_u2_shell + (std::pow((zheavy2_shell-zheavy2),2))*cz_asymm2_shell;

  Bsym = 0.0;
  Basym_1 = Bsym + std::pow((zheavy1-zsymm), 2) * cz_symm + delta_u1;
  Basym_2 = Bsym + std::pow((zheavy2-zsymm), 2) * cz_symm + delta_u2;
  if(Bsym < Basym_1 && Bsym < Basym_2) {
    // Excitation energies at the saddle point
    // without and with shell effect
    epsilon0_1_saddle = (e + e_zero_point - std::pow((zheavy1 - zsymm), 2) * cz_symm);
    epsilon0_2_saddle = (e + e_zero_point - std::pow((zheavy2 - zsymm), 2) * cz_symm);

    epsilon_1_saddle = epsilon0_1_saddle - delta_u1;
    epsilon_2_saddle = epsilon0_2_saddle - delta_u2;

    epsilon_symm_saddle = e + e_zero_point;
    eexc1_saddle = epsilon_1_saddle;
    eexc2_saddle = epsilon_2_saddle;

    // Excitation energies at the scission point
    // heavy fragment without and with shell effect
    epsilon0_1_scission = (e + e_saddle_scission - std::pow((zheavy1 - zsymm), 2) * cz_symm) * aheavy1/a;
    epsilon_1_scission = epsilon0_1_scission - delta_u1*(aheavy1/a);

    epsilon0_2_scission = (e + e_saddle_scission - std::pow((zheavy2 - zsymm), 2) * cz_symm) * aheavy2/a;
    epsilon_2_scission = epsilon0_2_scission - delta_u2*(aheavy2/a);

    epsilon_symm_scission = e + e_saddle_scission;
  } else if (Basym_1 < Bsym && Basym_1 < Basym_2) {
    // Excitation energies at the saddle point
    // without and with shell effect
    epsilon_symm_saddle = (e + e_zero_point + delta_u1 + std::pow((zheavy1-zsymm), 2) * cz_symm);
    epsilon0_2_saddle = (epsilon_symm_saddle - std::pow((zheavy2-zsymm), 2) * cz_symm);
    epsilon_2_saddle = epsilon0_2_saddle - delta_u2;
    epsilon0_1_saddle = e + e_zero_point + delta_u1;
    epsilon_1_saddle = e + e_zero_point;
    eexc1_saddle = epsilon_1_saddle;
    eexc2_saddle = epsilon_2_saddle;

    // Excitation energies at the scission point
    // heavy fragment without and with shell effect
    epsilon_symm_scission = (e + e_saddle_scission + std::pow((zheavy1-zsymm), 2) * cz_symm + delta_u1);
    epsilon0_2_scission = (epsilon_symm_scission - std::pow((zheavy2-zsymm), 2) * cz_symm) * aheavy2/a;
    epsilon_2_scission = epsilon0_2_scission - delta_u2*aheavy2/a;
    epsilon0_1_scission = (e + e_saddle_scission + delta_u1) * aheavy1/a;
    epsilon_1_scission = (e + e_saddle_scission) * aheavy1/a;
  } else if (Basym_2 < Bsym && Basym_2 < Basym_1) {
    // Excitation energies at the saddle point
    // without and with shell effect
    epsilon_symm_saddle = (e + e_zero_point + delta_u2 + std::pow((zheavy2-zsymm), 2) * cz_symm);
    epsilon0_1_saddle = (epsilon_symm_saddle - std::pow((zheavy1-zsymm), 2) * cz_symm);
    epsilon_1_saddle = epsilon0_1_saddle - delta_u1;
    epsilon0_2_saddle = e + e_zero_point + delta_u2;
    epsilon_2_saddle = e + e_zero_point;
    eexc1_saddle = epsilon_1_saddle;
    eexc2_saddle = epsilon_2_saddle;

    // Excitation energies at the scission point
    // heavy fragment without and with shell effect
    epsilon_symm_scission = (e + e_saddle_scission + std::pow((zheavy2-zsymm), 2) * cz_symm + delta_u2);
    epsilon0_1_scission = (epsilon_symm_scission - std::pow((zheavy1-zsymm), 2) * cz_symm) * aheavy1/a;
    epsilon_1_scission = epsilon0_1_scission - delta_u1*aheavy1/a;
    epsilon0_2_scission = (e + e_saddle_scission + delta_u2) * aheavy2/a;
    epsilon_2_scission = (e + e_saddle_scission) * aheavy2/a;

  } else {
    // G4cout <<"G4AblaFission: " << G4endl;
  }
  if(epsilon_1_saddle < 0.0) epsilon_1_saddle = 0.0;
  if(epsilon_2_saddle < 0.0) epsilon_2_saddle = 0.0;
  if(epsilon0_1_saddle < 0.0) epsilon0_1_saddle = 0.0;
  if(epsilon0_2_saddle < 0.0) epsilon0_2_saddle = 0.0;
  if(epsilon_symm_saddle < 0.0) epsilon_symm_saddle = 0.0;

  if(epsilon_1_scission < 0.0) epsilon_1_scission = 0.0;
  if(epsilon_2_scission < 0.0) epsilon_2_scission = 0.0;
  if(epsilon0_1_scission < 0.0) epsilon0_1_scission = 0.0;
  if(epsilon0_2_scission < 0.0) epsilon0_2_scission = 0.0;
  if(epsilon_symm_scission < 0.0) epsilon_symm_scission = 0.0;

  if(itest == 1) {
    // G4cout <<"E, E1, E2, Es" << e << epsilon_1_saddle << epsilon_2_saddle << epsilon_symm_saddle << G4endl;
  }

  e_eff1_saddle = epsilon0_1_saddle - delta_u1 * (std::exp((-epsilon_1_saddle*gamma)));
   
  if (e_eff1_saddle > 0.0) {
    wzasymm1_saddle = std::sqrt( (0.5) * 
				 (std::sqrt(1.0/a_levdens*e_eff1_saddle)) /
				 (cz_asymm1_shell * std::exp((-epsilon_1_saddle*gamma)) + cz_symm) );
  } else {
    wzasymm1_saddle = 1.0;
  }

  e_eff2_saddle = epsilon0_2_saddle - delta_u2 * std::exp((-epsilon_2_saddle*gamma));
  if (e_eff2_saddle > 0.0) {
    wzasymm2_saddle = std::sqrt( (0.5 * 
 				  (std::sqrt(1.0/a_levdens*e_eff2_saddle)) /
 				  (cz_asymm2_shell * std::exp((-epsilon_2_saddle*gamma)) + cz_symm) ) );
  } else {
    wzasymm2_saddle = 1.0;
  }

  if(e - e_zero_point > 0.0) {
    wzsymm_saddle = std::sqrt( (0.5 * 
 				(std::sqrt(1.0/a_levdens*(e+epsilon_symm_saddle))) / cz_symm ) );
  } else {
    wzsymm_saddle = 1.0;
  }

//   if (itest == 1) {
//     // G4cout << "wz1(saddle) = " << wzasymm1_saddle << G4endl;
//     // G4cout << "wz2(saddle) = " << wzasymm2_saddle << G4endl;
//     // G4cout << "wzsymm(saddle) = " << wzsymm_saddle << G4endl;
//   }
     
  //     /* Calculate widhts at the scission point: */
  //     /* fits of ref. Beizin 1991 (Plots brought to GSI by Sergei Zhdanov) */

  wzsymm_scission = wzsymm_saddle;

  if (e_saddle_scission == 0.0) {
    wzasymm1_scission = wzasymm1_saddle;
    wzasymm2_scission = wzasymm2_saddle;
  } else {
    if (nheavy1_eff > 75.0) {
      wzasymm1_scission = (std::sqrt(21.0)) * z/a;
      double RR = (70.0-28.0)/3.0*(z*z/a-35.0)+28.0;
      if(RR > 0.0) {
	wzasymm2_scission = std::sqrt(RR)*(z/a);
      } else {
	wzasymm2_scission = 0.0;
      }
      wzasymm2_scission = (std::sqrt (max( (70.0-28.0)/3.0*(z*z/a-35.0)+28.,0.0 )) ) * z/a;
    } else {
      wzasymm1_scission = wzasymm1_saddle;
      wzasymm2_scission = wzasymm2_saddle;
    }
  }

  wzasymm1_scission = max(wzasymm1_scission,wzasymm1_saddle);
  wzasymm2_scission = max(wzasymm2_scission,wzasymm2_saddle);

  wzasymm1 = wzasymm1_scission * fwidth_asymm1;
  wzasymm2 = wzasymm2_scission * fwidth_asymm2;
  wzsymm = wzsymm_scission;

//   /*      if (ITEST == 1) {
// 	  // G4cout << "WZ1(scission) = " << WZasymm1_scission << G4endl;
// 	  // G4cout << "WZ2(scission) = " << WZasymm2_scission << G4endl;
// 	  // G4cout << "WZsymm(scission) = " << WZsymm_scission << G4endl;
// 	  }
// 	  if (ITEST == 1) {
// 	  // G4cout << "WZ1(scission) final= " << WZasymm1 << G4endl;
// 	  // G4cout << "WZ2(scission) final= " << WZasymm2 << G4endl;
// 	  // G4cout << "WZsymm(scission) final= " << WZsymm << G4endl;
// 	  } */
      

  //  // G4cout <<"al, e, es, cn " << a_levdens << e << e_saddle_scission << cn << G4endl;

  //   if (itest == 1) {
  //     // G4cout << "wasymm = " << wzsymm << G4endl;
  //     // G4cout << "waheavy1 = " << waheavy1 << G4endl;
  //     // G4cout << "waheavy2 = " << waheavy2 << G4endl;
  //   }
            
  // sig_0 = quantum fluctuation = 0.45 z units for A=cte
  //                               0.45*2.58 = 1.16 n units for Z=cte
  //                     sig_0^2 = 1.16*2 = 1.35    n units for Z=cte
  n1width = std::sqrt(0.5 * std::sqrt(1.0/a_levdens*(e + e_saddle_scission)) / cn + 1.35);
  if ( (epsilon0_1_saddle - delta_u1*std::exp((-epsilon_1_saddle*gamma_heavy1))) < 0.0) {
    sasymm1 = -10.0;
  } else {
    sasymm1 = 2.0 * std::sqrt( a_levdens * (epsilon0_1_saddle - 
					    delta_u1*(std::exp((-epsilon_1_saddle*gamma_heavy1))) ) );
  }

  if ( (epsilon0_2_saddle - delta_u2*std::exp((-epsilon_2_saddle*gamma_heavy2))) < 0.0) {
    sasymm2 = -10.0;
  } else {
    sasymm2 = 2.0 * std::sqrt( a_levdens * (epsilon0_2_saddle - 
					    delta_u2*(std::exp((-epsilon_2_saddle*gamma_heavy2))) ) );
  }
              
  if (epsilon_symm_saddle > 0.0) {
    ssymm = 2.0 * std::sqrt( a_levdens*(epsilon_symm_saddle) );
  } else {
    ssymm = -10.0;
  }
      
  if (ssymm > -10.0) {
    ysymm = 1.0;
    if (epsilon0_1_saddle < 0.0) { //  /* low energy */
      yasymm1 = std::exp((sasymm1-ssymm)) * wzasymm1_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } else { //        /* high energy */
      ssymm_mode1 = 2.0 * std::sqrt( a_levdens*(epsilon0_1_saddle) );
      yasymm1 = ( std::exp((sasymm1-ssymm)) - std::exp((ssymm_mode1 - ssymm)) )  
 	* wzasymm1_saddle / wzsymm_saddle * 2.0;
    }

    if (epsilon0_2_saddle < 0.0) { //  /* low energy */
      yasymm2 = std::exp((sasymm2-ssymm)) * wzasymm2_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } else { //        /* high energy */
      ssymm_mode2 = 2.0 * std::sqrt( a_levdens*(epsilon0_2_saddle) );
      yasymm2 = ( std::exp((sasymm2-ssymm)) - std::exp((ssymm_mode2 - ssymm)) )  
 	* wzasymm2_saddle / wzsymm_saddle * 2.0;
    }       
    //                            /* difference in the exponent in order */
    //                            /* to avoid numerical overflow         */
   } 
  else {
    if ( (sasymm1 > -10.0) && (sasymm2 > -10.0) ) {
      ysymm = 0.0;
      yasymm1 = std::exp(sasymm1) * wzasymm1_saddle * 2.0;
      yasymm2 = std::exp(sasymm2) * wzasymm2_saddle * 2.0;
    }
  }

  //  /* normalize */
  ysum = ysymm + yasymm1 + yasymm2;
  if (ysum > 0.0) {
    ysymm = ysymm / ysum;
    yasymm1 = yasymm1 / ysum;
    yasymm2 = yasymm2 / ysum;
  } else {
    ysymm = 0.0;
    yasymm1 = 0.0;
    yasymm2 = 0.0;
    //        /* search minimum threshold and attribute all events to this mode */
    if ( (epsilon_symm_saddle < epsilon_1_saddle) && (epsilon_symm_saddle < epsilon_2_saddle) ) {
      ysymm = 1.0;
    } else {
      if (epsilon_1_saddle < epsilon_2_saddle) {
	yasymm1 = 1.0;
      } else {
	yasymm2 = 1.0;
      }
    }
  }

//   if (itest == 1) {
//     // G4cout << "ysymm normalized= " << ysymm  << G4endl;
//     // G4cout << "yasymm1 normalized= " << yasymm1  << G4endl;
//     // G4cout << "yasymm2 normalized= " << yasymm2  << G4endl;
//   }
      
  //      /* even-odd effect */
  //      /* simple parametrization KHS, Nov. 2000. From Rejmund et al. */
  eexc_max = max(eexc1_saddle, eexc2_saddle);
  eexc_max = max(eexc_max, e);
  //  // G4cout << "mod(z, 2)" << iz%2 << G4endl;
  if ((G4int)(z) % 2 == 0) {
    r_e_o_max = 0.3 * (1.0 - 0.2 * (std::pow(z, 2)/a - std::pow(92.0, 2)/238.0));
    e_pair = 2.0 * 12.0 / std::sqrt(a);
    if(eexc_max > (e_crit + e_pair)) {
      r_e_o = 0.0;
    } else {
      if(eexc_max < e_pair) {
	r_e_o = r_e_o_max;
      } else {
	r_e_o = std::pow((eexc_max - e_crit - e_pair)/e_crit, 2) * r_e_o_max;
      }
    }
  } else {
    r_e_o = 0.0;
  }

  // // G4cout <<"rmax " << r_e_o_max << G4endl;
  // if(r_e_o > 0.0) // G4cout <<"e_crit, r_e_o" << e_crit << r_e_o << G4endl;
  //      $LOOP;    /* event loop */
  //     I_COUNT = I_COUNT + 1;

  /* random decision: symmetric or asymmetric */
  /* IMODE = 3 means asymmetric fission, mode 1,
     IMODE = 2 means asymmetric fission, mode 2,
     IMODE = 1 means symmetric  */
  //  RMODE = dble(HAZ(k));
  //      rmode = rnd.rndm();  
//   // Safety check added to make sure we always select well defined
//   // fission mode.
    rmode = haz(k);
    // Cast for test CS 11/10/05
    //      RMODE = 0.54;    
    //  rmode = 0.54;
    if (rmode < ysymm) {
      imode = 1;
    } else if (rmode < (ysymm + yasymm1)) {
      imode = 2;
    } else {
      imode = 3;
    }
    //     /* determine parameters of the Z distribution */
    // force imode (for testing, PK)
    // imode = 3;
    
    if (imode == 1) {
      z1mean = zsymm;
      z1width = wzsymm;
    } else if (imode == 2) {
      z1mean = zheavy1;
      z1width = wzasymm1;
    } else if (imode == 3) {
      z1mean = zheavy2;
      z1width = wzasymm2;
    }

    if (itest == 1) {
      // G4cout << "nbre aleatoire tire " << rmode << G4endl;
      // G4cout << "fission mode " << imode << G4endl;
      // G4cout << "z1mean= " << z1mean << G4endl;
      // G4cout << "z1width= " << z1width << G4endl;
    }
		      
    //     /* random decision: Z1 and Z2 at scission: */
    z1 = 1.0;
    z2 = 1.0;

    while  ( (z1<5.0) || (z2<5.0) ) {
      //         Z1 = dble(GAUSSHAZ(K,sngl(Z1mean),sngl(Z1width)));
      //	 z1 = rnd.gaus(z1mean,z1width);
      //      z1 = 48.26; // gausshaz(k, z1mean, z1width);
      z1 = gausshaz(k, z1mean, z1width);
      even_odd(z1, r_e_o, i_help);
      z1 = double(i_help);
      z2 = z - z1;
    }

    if (itest == 1) {
      // G4cout << "ff charge sample " << G4endl;
      // G4cout << "z1 =  " << z1 << G4endl;
      // G4cout << "z2 = " << z2 << G4endl;
    }

//   //     CALL EVEN_ODD(Z1,R_E_O,I_HELP);
//   //         /* Integer proton number with even-odd effect */
//   //     Z1 = REAL(I_HELP)
//   //      /* Z1 = INT(Z1+0.5E0); */
//   z2 = z - z1;

    //     /* average N of both fragments: */
    if (imode == 1) {
      n1ucd = z1 * n/z;
      n2ucd = z2 * n/z;
      re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) + ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,2.0); // umass == massdef
      re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) + ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,2.0);
      re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) + ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,2.0);
      reps2 = (re1-2.0*re2+re3) / 2.0;
      reps1 = re2 - re1 - reps2;
      rn1_pol = - reps1 / (2.0 * reps2);
      n1mean = n1ucd + rn1_pol;
    } else {
      n1mean = (z1 + 0.5) * n/z;
    }
      
// n1mean nsymm + (z1 - zsymm) * 1.6 from 238 U(nth, f)
// n1width = 0.9 + E * 0.002 KHS

// random decision: N1R and N2R at scission, before evaporation
    n1r = 1.0;
    n2r = 1.0;
    while (n1r < 5 || n2r < 5) {
      //      n1r = 76.93; gausshaz(kkk,n1mean,n1width);
      n1r = gausshaz(kkk,n1mean,n1width);
      // modification to have n1r as integer, and n=n1r+n2r rigorously a.b. 19/4/2001
      i_inter = int(n1r + 0.5);
      n1r = double(i_inter);
      n2r = n - n1r;
    }

    // neutron evaporation from fragments
    if (i_eva > 0) {
      // treatment sz
      ne_min = 0.095e0 * a - 20.4e0;                                  
      if (ne_min < 0) ne_min = 0.0;                                 
      ne_min = ne_min + e / 8.e0; // 1 neutron per 8 mev */           
    }

    // excitation energy due to deformation                                 

    a1 = z1 + n1r;         // mass of first fragment */        
    a2 = z2 + n2r;      // mass of second fragment */        
    if (a1 < 80) {
      ed1_low = 0.0;
    } else if (a1 >= 80 && a1 < 110) {
      ed1_low = (a1-80.)*20./30.;
    } else if (a1 >= 110 && a1 < 130) {
      ed1_low = -(a1-110.)*20./20. + 20.;
    } else if (a1 >= 130) {
      ed1_low = (a1-130.)*20./30.;
    }
    
    if (a2 < 80) {
      ed2_low = 0.0;
    } else if (a2 >= 80 && a2 < 110) {
      ed2_low = (a2-80.)*20./30.;
    } else if (a2 >= 110 && a2 < 130) {
      ed2_low = -(a2-110.)*20./20. + 20.;
    } else if (a2 >= 130) {
      ed2_low = (a2-130.)*20./30.;
    }

    ed1_high = 20.0*a1/(a1+a2);
    ed2_high = 20.0 - ed1_high;
    ed1 = ed1_low*std::exp(-e/el) + ed1_high*(1-std::exp(-e/el));
    ed2 = ed2_low*std::exp(-e/el) + ed2_high*(1-std::exp(-e/el));

    //  write(6,101)e,a1,a2,ed1,ed2,ed1+ed2
    //  write(6,102)ed1_low,ed1_high,ed2_low,ed2_high
    e1 = e*a1/(a1+a2) + ed1;
    e2 = e - e*a1/(a1+a2) + ed2;
    atot = a1+a2;
    if (atot > a+1) {
      // write(6,*)'a,,a1,a2,atot',a,a1,a2,atot
      // write(6,*)'n,n1r,n2r',n,n1r,n2r
      // write(6,*)'z,z1,z2',z,z1,z2
    }

 milledeux:       
    // only symmetric fission
    // Symmetric fission: Ok! Checked CS 10/10/05
    if ( (icz == -1) || (a1 < 0.0) || (a2 < 0.0) ) {
      //           IF (z.eq.92) THEN
      //              write(6,*)'symmetric fission'
      //              write(6,*)'Z,A,E,A1,A2,icz,Atot',Z,A,E,A1,A2,icz,Atot
      //           END IF

      if (itest == 1) {
	// G4cout << "milledeux: liquid-drop option "  << G4endl;
      }

      n = a-z;
      //  proton number in symmetric fission (centre) *
      zsymm  = z / 2.0;
      nsymm  = n / 2.0;
      asymm = nsymm + zsymm;

      a_levdens = a / 8.0;

      masscurv = 2.0;
      cz_symm = 8.0 / std::pow(z,2) * masscurv;

      wzsymm = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cz_symm) ) ;

      if (itest == 1) {
	// G4cout << " symmetric high energy fission " << G4endl;
	// G4cout << "wzsymm " << wzsymm << G4endl;
      }

      z1mean = zsymm;
      z1width = wzsymm;

      // random decision: Z1 and Z2 at scission: */
      z1 = 1.0;
      z2 = 1.0;
      while  ( (z1 < 5.0) || (z2 < 5.0) ) {
	//           z1 = dble(gausshaz(kkk,sngl(z1mean),sngl(z1width)));
	//	   z1 = rnd.gaus(z1mean,z1width);
	//	z1 = 24.8205585; //gausshaz(kkk, z1mean, z1width);
	z1 = gausshaz(kkk, z1mean, z1width);
	z2 = z - z1;
      }

      if (itest == 1) {
	// G4cout << " z1 " << z1 << G4endl;
	// G4cout << " z2 " << z2 << G4endl;
      }
      if (itest == 1) {
	// G4cout << " zsymm " << zsymm << G4endl;
	// G4cout << " nsymm " << nsymm << G4endl;
	// G4cout << " asymm " << asymm << G4endl;
      }

      cn = umass(zsymm, nsymm+1.0, 0.0) + umass(zsymm, nsymm-1.0, 0.0)
	+ 1.44 * std::pow(zsymm, 2)/
	(std::pow(r_null, 2) * std::pow(std::pow(asymm+1.0, 1.0/3.0) + std::pow(asymm-1.0, 1.0/3.0), 2))
	- 2.0 * umass(zsymm, nsymm, 0.0) - 1.44 * std::pow(zsymm, 2) /
	std::pow(r_null * 2.0 *std::pow(asymm, 1.0/3.0), 2);
      //      This is an approximation! Coulomb energy is neglected.

      n1width = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cn) + 1.35);
      if (itest == 1) {
	// G4cout << " cn " << cn << G4endl;
	// G4cout << " n1width " << n1width << G4endl;
      }
	
    //     /* average N of both fragments: */
      n1ucd = z1 * n/z;
      n2ucd = z2 * n/z;
      re1 = umass(z1,n1ucd,   0.6) + umass(z2,n2ucd,   0.6) + ecoul(z1,n1ucd,   0.6,z2,n2ucd,   0.6,2.0);
      re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) + ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,2.0);
      re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) + ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,2.0);
      reps2 = (re1-2.0*re2+re3) / 2.0;
      reps1 = re2 - re1 - reps2;
      rn1_pol = - reps1 / (2.0 * reps2);
      n1mean = n1ucd + rn1_pol;

      // random decision: N1R and N2R at scission, before evaporation: */
      //       N1R = dfloat(NINT(GAUSSHAZ(KKK,sngl(N1mean),sngl(N1width))));
      //      n1r = (float)( (int)(rnd.gaus(n1mean,n1width)) );
      //      n1r = 34.0; //(float)( (int)(gausshaz(k, n1mean,n1width)) );
      n1r = (float)( (int)(gausshaz(k, n1mean,n1width)) );
      n2r = n - n1r;
      // Mass of first and second fragment */
      a1 = z1 + n1r;
      a2 = z2 + n2r;

      e1 = e*a1/(a1+a2);
      e2 = e - e*a1/(a1+a2);
    }
    v1 = 0.0; // These are not calculated in SimFis3.
    v2 = 0.0; 
    if (itest == 1) {
      // G4cout << " n1r " << n1r << G4endl;
      // G4cout << " n2r " << n2r << G4endl;
    }

    if (itest == 1) {
      // G4cout << " a1 " << a1 << G4endl;
      // G4cout << " z1 " << z1 << G4endl;
      // G4cout << " a2 " << a2 << G4endl;
      // G4cout << " z2 " << z2 << G4endl;
      // G4cout << " e1 " << e1 << G4endl;
      // G4cout << " e2 " << e << G4endl;
    }
}

G4double G4AblaFission::haz(G4int k)
{
  const G4int pSize = 110;
  static G4ThreadLocal G4double p[pSize];
  static G4ThreadLocal G4long ix = 0, i = 0;
  static G4ThreadLocal G4double x = 0.0, y = 0.0, a = 0.0, fhaz = 0.0;
  //  k =< -1 on initialise                                        
  //  k = -1 c'est reproductible                                   
  //  k < -1 || k > -1 ce n'est pas reproductible

  // Zero is invalid random seed. Set proper value from our random seed collection:
  if(ix == 0) {
    //    ix = hazard->ial;
  }

  if (k <= -1) { //then                                             
    if(k == -1) { //then                                            
      ix = 0;
    }
    else {
      x = 0.0;
      y = secnds(int(x));
      ix = int(y * 100 + 43543000);
      if(mod(ix,2) == 0) {
	ix = ix + 1;
      }
    }

    // Here we are using random number generator copied from INCL code
    // instead of the CERNLIB one! This causes difficulties for
    // automatic testing since the random number generators, and thus
    // the behavior of the routines in C++ and FORTRAN versions is no
    // longer exactly the same!
    x = G4AblaRandom::flat();
    //    standardRandom(&x, &ix);
    for(G4int iRandom = 0; iRandom < pSize; iRandom++) { //do i=1,110                                                 
      p[iRandom] = G4AblaRandom::flat();
      //      standardRandom(&(p[i]), &ix);
    }
    a = G4AblaRandom::flat();
    //standardRandom(&a, &ix);
    k = 0;
  }

  i = nint(100*a)+1;
  fhaz = p[i];
  a = G4AblaRandom::flat();
  //  standardRandom(&a, &ix);
  p[i] = a;

  //  hazard->ial = ix;
  //  haz=0.4;
  return fhaz;
}

G4double G4AblaFission::gausshaz(int k, double xmoy, double sig)
{
  // Gaussian random numbers:

  //   1005       C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
  static G4ThreadLocal G4int  iset = 0;
  static G4ThreadLocal G4double v1,v2,r,fac,gset,fgausshaz;

  if(iset == 0) { //then                                              
    do {
      v1 = 2.0*haz(k) - 1.0;
      v2 = 2.0*haz(k) - 1.0;
      r = std::pow(v1,2) + std::pow(v2,2);
    } while(r >= 1);

    fac = std::sqrt(-2.*std::log(r)/r);
    gset = v1*fac;
    fgausshaz = v2*fac*sig+xmoy;
    iset = 1;
  }
  else {
    fgausshaz=gset*sig+xmoy;
    iset=0;
  }
  return fgausshaz;                                                         
}

// Utilities

G4double G4AblaFission::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFission::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4AblaFission::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFission::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFission::nint(G4double number)
{
  G4double intpart = 0.0;
  G4double fractpart = 0.0;
  fractpart = std::modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return int(std::floor(number));
    }
    else {
      return int(std::ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return int(std::floor(number));
    }
    else {
      return int(std::ceil(number));
    }
  }

  return int(std::floor(number));
}

G4int G4AblaFission::secnds(G4int x)
{
  time_t mytime;
  tm *mylocaltime;

  time(&mytime);
  mylocaltime = localtime(&mytime);

  if(x == 0) {
    return(mylocaltime->tm_hour*60*60 + mylocaltime->tm_min*60 + mylocaltime->tm_sec);
  }
  else {
    return(mytime - x);
  }
}

G4int G4AblaFission::mod(G4int a, G4int b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0;
  } 
}

G4double G4AblaFission::dmod(G4double a, G4double b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0.0;
  } 
}

G4double G4AblaFission::dint(G4double a)
{
  G4double value = 0.0;

  if(a < 0.0) {
    value = double(std::ceil(a));
  }
  else {
    value = double(std::floor(a));
  }

  return value;
}

G4int G4AblaFission::idint(G4double a)
{
  G4int value = 0;

  if(a < 0) {
    value = int(std::ceil(a));
  }
  else {
    value = int(std::floor(a));
  }

  return value;
}

G4double G4AblaFission::dmin1(G4double a, G4double b, G4double c)
{
  if(a < b && a < c) {
    return a;
  }
  if(b < a && b < c) {
    return b;
  }
  if(c < a && c < b) {
    return c;
  }
  return a;
}

G4double G4AblaFission::utilabs(G4double a)
{
  if(a > 0) {
    return a;
  }
  if(a < 0) {
    return (-1*a);
  }
  if(a == 0) {
    return a;
  }

  return a;
}

