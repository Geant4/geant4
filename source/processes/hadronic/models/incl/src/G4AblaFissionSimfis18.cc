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
// $Id: G4AblaFissionSimfis18.cc,v 1.5 2010-11-17 20:19:09 kaitanie Exp $
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4AblaFissionSimfis18.hh"
#include <time.h>

G4AblaFissionSimfis18::G4AblaFissionSimfis18()
{
  hazard = 0;
  randomGenerator = 0;
}

G4AblaFissionSimfis18::G4AblaFissionSimfis18(G4Hazard *hzr, G4InclRandomInterface *rndm)
{
  hazard = hzr;
  randomGenerator = rndm;
  setAboutString("Fission model: Based on ABLA with SimFis18");
}

G4AblaFissionSimfis18::~G4AblaFissionSimfis18()
{

}

void G4AblaFissionSimfis18::doFission(G4double &A, G4double &Z, G4double &E,
		 G4double &A1, G4double &Z1, G4double &E1, G4double &K1,
		 G4double &A2, G4double &Z2, G4double &E2, G4double &K2)
{
  fissionDistri(A,Z,E,A1,Z1,E1,K1,A2,Z2,E2,K2);
}

void G4AblaFissionSimfis18::even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out)     
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

G4double G4AblaFissionSimfis18::umass(G4double z,G4double n,G4double beta)
{
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
  // pure liquid drop, without pairing and shell effects

  // On input:    Z     nuclear charge of nucleus
  //              N     number of neutrons in nucleus
  //              beta  deformation of nucleus
  // On output:   binding energy of nucleus

  G4double a = 0.0, umass = 0.0;
  G4double alpha = 0.0;
  G4double xcom = 0.0, xvs = 0.0, xe = 0.0;
  const G4double pi = 3.1416;
     
  a = n + z;
  alpha = ( std::sqrt(5.0/(4.0*pi)) ) * beta;
  
  xcom = 1.0 - 1.7826 * ((a - 2.0*z)/a)*((a - 2.0*z)/a);
  // factor for asymmetry dependence of surface and volume term
  xvs = - xcom * ( 15.4941 * a - 
		   17.9439 * std::pow(a,0.66667) * (1.0+0.4*alpha*alpha) );
  // sum of volume and surface energy
  xe = z*z * (0.7053/(std::pow(a,0.33333)) * (1.0-0.2*alpha*alpha) - 1.1529/a);
  umass = xvs + xe;
  
  return umass;
}

G4double G4AblaFissionSimfis18::ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d)
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
  G4double ecoul = 0;
  G4double dtot = 0;
  const G4double r0 = 1.16;

  dtot = r0 * ( std::pow((z1+n1),0.33333) * (1.0+(2.0/3.0)*beta1)
		+ std::pow((z2+n2),0.33333) * (1.0+(2.0/3.0)*beta2) ) + d;
  ecoul = z1 * z2 * 1.44 / dtot;

  return ecoul;
}

void G4AblaFissionSimfis18::fissionDistri(G4double &a,G4double &z,G4double &e,
			   G4double &a1,G4double &z1,G4double &e1,G4double &v1,
			   G4double &a2,G4double &z2,G4double &e2,G4double &v2)
{
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
  G4double     nlight1 = 0.0, nlight2 = 0.0;
  G4double     aheavy1 = 0.0,alight1 = 0.0, aheavy2 = 0.0, alight2 = 0.0;
  G4double     eheavy1 = 0.0, elight1 = 0.0, eheavy2 = 0.0, elight2 = 0.0;
  G4double     zheavy1_shell = 0.0, zheavy2_shell = 0.0;
  G4double     zlight1 = 0.0, zlight2 = 0.0;
  G4double     masscurv = 0.0;
  G4double     sasymm1 = 0.0, sasymm2 = 0.0, ssymm = 0.0, ysum = 0.0, yasymm = 0.0;
  G4double     ssymm_mode1 = 0.0, ssymm_mode2 = 0.0;
  G4double     cz_asymm1_saddle = 0.0, cz_asymm2_saddle = 0.0;
  // Curvature at saddle, modified by ld-potential
  G4double     wzasymm1_saddle, wzasymm2_saddle, wzsymm_saddle  = 0.0;
  G4double     wzasymm1_scission = 0.0, wzasymm2_scission = 0.0, wzsymm_scission = 0.0;
  G4double     wzasymm1 = 0.0, wzasymm2 = 0.0, wzsymm = 0.0;
  G4double     nlight1_eff = 0.0, nlight2_eff = 0.0;
  G4int  imode = 0;
  G4double     rmode = 0.0;
  G4double     z1mean = 0.0, z2mean = 0.0, z1width = 0.0, za1width = 0.0;
  //      G4double     Z1,Z2,N1R,N2R,A1R,A2R,N1,N2,A1,A2;
  G4double     n1r = 0.0, n2r = 0.0, a1r = 0.0, a2r = 0.0, n1 = 0.0, n2 = 0.0;

  G4double     zsymm = 0.0, nsymm = 0.0, asymm = 0.0;
  G4double     n1mean = 0.0, n2mean, n1width;
  G4double     dueff = 0.0;
  // effective shell effect at lowest barrier
  G4double     eld = 0.0;
  // Excitation energy with respect to ld barrier
  G4double     re1 = 0.0, re2 = 0.0, re3 = 0.0;
  G4double     eps1 = 0.0, eps2 = 0.0;
  G4double     n1ucd = 0.0, n2ucd = 0.0, z1ucd = 0.0, z2ucd = 0.0;
  G4double     beta = 0.0, beta1 = 0.0, beta2 = 0.0;

  G4double     dn1_pol = 0.0;
  // shift of most probable neutron number for given Z,
  // according to polarization
  G4int  i_help = 0;

  //   /* Parameters of the semiempirical fission model */
  G4double a_levdens = 0.0;
  //           /* level-density parameter */
  G4double a_levdens_light1 = 0.0, a_levdens_light2 = 0.0;
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
  G4double epot0_mode1_saddle = 0.0, epot0_mode2_saddle = 0.0, epot0_symm_saddle = 0.0;
  G4double epot_mode1_saddle = 0.0, epot_mode2_saddle = 0.0, epot_symm_saddle = 0.0;
  G4double e_defo = 0.0, e_defo1 = 0.0, e_defo2 = 0.0, e_scission = 0.0, e_asym = 0.0;
  G4double e1exc = 0.0, e2exc = 0.0;
  G4double e1exc_sigma = 0.0, e2exc_sigma = 0.0;
  G4double e1final = 0.0, e2final = 0.0;

  const G4double r0 = 1.16;
  G4double tker = 0.0;
  G4double ekin1 = 0.0, ekin2 = 0.0;
  //      G4double EkinR1,EkinR2,E1,E2,V1,V2;
  G4double ekinr1 = 0.0, ekinr2 = 0.0;
  G4int icz = 0, k = 0;

  //   Input parameters:
  //OMMENT(Nuclear charge number);
  //      G4double Z;
  //OMMENT(Nuclear mass number);
  //      G4double A;
  //OMMENT(Excitation energy above fission barrier);
  //      G4double E;

  //   Model parameters:
  //OMMENT(position of heavy peak valley 1);
  const G4double nheavy1 = 83.0;
  //OMMENT(position of heavy peak valley 2);
  const G4double nheavy2 = 90.0;
  //OMMENT(Shell effect for valley 1);
  const G4double delta_u1_shell = -2.65;
  //        Parameter (Delta_U1_shell = -2)
  //OMMENT(Shell effect for valley 2);
  const G4double delta_u2_shell = -3.8;
  //        Parameter (Delta_U2_shell = -3.2)
  //OMMENT(I: used shell effect);
  G4double delta_u1 = 0.0;
  //omment(I: used shell effect);
  G4double delta_u2 = 0.0;
  //OMMENT(Curvature of asymmetric valley 1);
  const G4double cz_asymm1_shell = 0.7;
  //OMMENT(Curvature of asymmetric valley 2);
  const G4double cz_asymm2_shell = 0.15;
  //OMMENT(Factor for width of distr. valley 1);
  const G4double fwidth_asymm1 = 0.63;
  //OMMENT(Factor for width of distr. valley 2);
  const G4double fwidth_asymm2 = 0.97;
  //       Parameter (CZ_asymm2_scission = 0.12)
  //OMMENT(Parameter x: a = A/x);
  const G4double xlevdens = 12.0;
  //OMMENT(Factor to gamma_heavy1);
  const G4double fgamma1 = 2.0;
  //OMMENT(I: fading of shells (general));
  G4double gamma = 0.0;
  //OMMENT(I: fading of shell 1);
  G4double gamma_heavy1 = 0.0;
  //OMMENT(I: fading of shell 2);
  G4double gamma_heavy2 = 0.0;
  //OMMENT(Zero-point energy at saddle);
  const G4double e_zero_point = 0.5;
  //OMMENT(I: friction from saddle to scission);
  G4double e_saddle_scission = 0.0;
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
  //omment(I: Effective mass mode 1);
  G4double aheavy1_mean = 0.0;
  //omment(I: Effective mass mode 2);
  G4double aheavy2_mean = 0.0;
  //omment(I: Width of symmetric mode);
  G4double wasymm_saddle = 0.0;
  //OMMENT(I: Width of asymmetric mode 1);
  G4double waheavy1_saddle = 0.0;
  //OMMENT(I: Width of asymmetric mode 2);
  G4double waheavy2_saddle = 0.0;
  //omment(I: Width of symmetric mode);
  G4double wasymm = 0.0;
  //OMMENT(I: Width of asymmetric mode 1);
  G4double waheavy1 = 0.0;
  //OMMENT(I: Width of asymmetric mode 2);
  G4double waheavy2 = 0.0;
  //OMMENT(I: Even-odd effect in Z);
  G4double r_e_o = 0.0, r_e_o_exp = 0.0;
  //OMMENT(I: Curveture of symmetric valley);
  G4double cz_symm = 0.0;
  //OMMENT(I: Curvature of mass distribution for fixed Z);
  G4double cn = 0.0;
  //OMMENT(I: Curvature of Z distribution for fixed A);
  G4double cz = 0.0;
  //OMMENT(Minimum neutron width for constant Z);
  const G4double sigzmin = 1.16;
  //OMMENT(Surface distance of scission configuration);
  const G4double d = 2.0;

  //   /* Charge polarisation from Wagemanns p. 397: */
  //OMMENT(Charge polarisation standard I);
  const G4double cpol1 = 0.65;
  //OMMENT(Charge polarisation standard II);
  const G4double cpol2 = 0.55;
  //OMMENT(=1: Polarisation simult. in N and Z);
  const G4int nzpol = 1;
  //OMMENT(=1: test output, =0: no test output);
  const G4int itest = 0;
      
  //      G4double UMASS, ECOUL, reps1, reps2, rn1_pol;
  G4double reps1 = 0.0, reps2 = 0.0, rn1_pol = 0.0;
  //      Float_t HAZ,GAUSSHAZ;
  G4int kkk = 0;
  //  G4int kkk = 10; // PK
  
  //     I_MODE = 0;

  if(itest == 1) {
    G4cout << " cn mass " << a << G4endl;
    G4cout << " cn charge " << z << G4endl;
    G4cout << " cn energy " << e << G4endl;
  }

  //     /* average Z of asymmetric and symmetric components: */
  n = a - z;  /* neutron number of the fissioning nucleus */

  k = 0;
  icz = 0;
  if ( (std::pow(z,2)/a < 25.0) || (n < nheavy2) || (e > 500.0) ) {
    icz = -1;
    //          GOTO 1002;
    goto milledeux;
  }

  nlight1 = n - nheavy1;
  nlight2 = n - nheavy2;
	
  //    /* Polarisation assumed for standard I and standard II:
  //      Z - Zucd = cpol (for A = const);
  //      from this we get (see Armbruster)
  //      Z - Zucd =  Acn/Ncn * cpol (for N = const)                        */

  zheavy1_shell = ((nheavy1/n) * z) - ((a/n) * cpol1);
  zheavy2_shell = ((nheavy2/n) * z) - ((a/n) * cpol2);

  e_saddle_scission = 
    (-24.0 + 0.02227 * (std::pow(z,2))/(std::pow(a,0.33333)) ) * friction_factor;
    
  //      /* Energy dissipated from saddle to scission                        */
  //      /* F. Rejmund et al., Nucl. Phys. A 678 (2000) 215, fig. 4 b        */
  //      E_saddle_scission = DMAX1(0.,E_saddle_scission);
  if (e_saddle_scission > 0.) {
    e_saddle_scission = e_saddle_scission;
  }
  else {
    e_saddle_scission = 0.;
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
  //  /* New parametrization of T. Enqvist according to Mulgin et al. 1998 */
  if ( (std::pow(z,2))/a < 34.0) {
    masscurv =  std::pow( 10.0,(-1.093364 + 0.082933 * (std::pow(z,2)/a)
			   - 0.0002602 * (std::pow(z,4)/std::pow(a,2))) );
  } else {
    masscurv = std::pow( 10.0,(3.053536 - 0.056477 * (std::pow(z,2)/a)
			  + 0.0002454 * (std::pow(z,4)/std::pow(a,2))) );
  }

  cz_symm = (8.0/std::pow(z,2)) * masscurv;

  if(itest == 1) {
    G4cout << "cz_symmetry= " << cz_symm << G4endl;
  }

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
  nheavy1_eff = (zheavy1 + (a/n * cpol1))*(n/z);
  nheavy2_eff = (zheavy2 + (a/n * cpol2))*(n/z);
  nlight1_eff = n - nheavy1_eff;
  nlight2_eff = n - nheavy2_eff;
  //  /* proton number of light fragments (centre) */
  zlight1 = z - zheavy1;
  //  /* proton number of light fragments (centre) */
  zlight2 = z - zheavy2;
  aheavy1 = nheavy1_eff + zheavy1;
  aheavy2 = nheavy2_eff + zheavy2;
  aheavy1_mean = aheavy1;
  aheavy2_mean = aheavy2;
  alight1 = nlight1_eff + zlight1;
  alight2 = nlight2_eff + zlight2;

  a_levdens = a / xlevdens;
  a_levdens_heavy1 = aheavy1 / xlevdens;
  a_levdens_heavy2 = aheavy2 / xlevdens;
  a_levdens_light1 = alight1 / xlevdens;
  a_levdens_light2 = alight2 / xlevdens;
  gamma = a_levdens / (0.4 * (std::pow(a,1.3333)) );
  gamma_heavy1 = ( a_levdens_heavy1 / (0.4 * (std::pow(aheavy1,1.3333)) ) ) * fgamma1;
  gamma_heavy2 = a_levdens_heavy2 / (0.4 * (std::pow(aheavy2,1.3333)) );

  cz_asymm1_saddle = cz_asymm1_shell + cz_symm;
  cz_asymm2_saddle = cz_asymm2_shell + cz_symm;
	
  // Up to here: Ok! Checked CS 10/10/05      	   

  cn = umass(zsymm,(nsymm+1.),0.0) + umass(zsymm,(nsymm-1.),0.0)
    + 1.44 * (std::pow(zsymm,2))/
    ( (std::pow(r_null,2)) * 
      ( std::pow((asymm+1.0),0.33333) + std::pow((asymm-1.0),0.33333) ) *
      ( std::pow((asymm+1.0),0.33333) + std::pow((asymm-1.0),0.33333) ) )
    - 2.0 * umass(zsymm,nsymm,0.0)
    - 1.44 * (std::pow(zsymm,2))/
    ( ( 2.0 * r_null * (std::pow(asymm,0.33333)) ) * 
      ( 2.0 * r_null * (std::pow(asymm,0.33333)) ) );
	
  // /* shell effect in valley of mode 1 */
  delta_u1 = delta_u1_shell + (std::pow((zheavy1_shell-zheavy1),2))*cz_asymm1_shell;
  // /* shell effect in valley of mode 2 */
  delta_u2 = delta_u2_shell + (std::pow((zheavy2_shell-zheavy2),2))*cz_asymm2_shell;

  //     /* liquid drop energies
  //        at the centres of the different shell effects
  //        with respect to liquid drop at symmetry: */
  epot0_mode1_saddle = (std::pow((zheavy1-zsymm),2)) * cz_symm;
  epot0_mode2_saddle = (std::pow((zheavy2-zsymm),2)) * cz_symm;
  epot0_symm_saddle = 0.0;
      
  if (itest == 1) {
    G4cout << "check zheavy1 = " << zheavy1  << G4endl;
    G4cout << "check zheavy2 = " << zheavy2  << G4endl;
    G4cout << "check zsymm = " << zsymm  << G4endl;
    G4cout << "check czsymm = " << cz_symm  << G4endl;
    G4cout << "check epot0_mode1_saddle = " << epot0_mode1_saddle  << G4endl;
    G4cout << "check epot0_mode2_saddle = " << epot0_mode2_saddle  << G4endl;
    G4cout << "check epot0_symm_saddle = " << epot0_symm_saddle  << G4endl;
    G4cout << "delta_u1 = " << delta_u1 << G4endl;
    G4cout << "delta_u2 = " << delta_u2 << G4endl;
  }
      
  //     /* energies including shell effects
  //        at the centres of the different shell effects
  //        with respect to liquid drop at symmetry: */
  epot_mode1_saddle = epot0_mode1_saddle + delta_u1;
  epot_mode2_saddle = epot0_mode2_saddle + delta_u2;
  epot_symm_saddle = epot0_symm_saddle;
  if (itest == 1) {
    G4cout << "check epot_mode1_saddle = " << epot_mode1_saddle  << G4endl;
    G4cout << "check epot_mode2_saddle = " << epot_mode2_saddle  << G4endl;
    G4cout << "check epot_symm_saddle = " << epot_symm_saddle  << G4endl;
  }

  //     /* Minimum of potential with respect to ld potential at symmetry */
  dueff = min(epot_mode1_saddle,epot_mode2_saddle);
  dueff = min(dueff,epot_symm_saddle);
  dueff = dueff - epot_symm_saddle;

  eld = e + dueff + e_zero_point;
      
  if (itest == 1) {
    G4cout << "check dueff = " << dueff  << G4endl;
    G4cout << "check e = " << e  << G4endl;
    G4cout << "check e_zero_point = " << e_zero_point  << G4endl;
    G4cout << "check eld = " << eld  << G4endl;
  }
  // Up to here: Ok! Checked CS 10/10/05 
     
  //          /* E = energy above lowest effective barrier */
  //          /* Eld = energy above liquid-drop barrier */

  //     /* Due to this treatment the energy E on input means the excitation  */
  //     /* energy above the lowest saddle.                                   */

  //  /* These energies are not used */
  eheavy1 = e * aheavy1 / a;
  eheavy2 = e * aheavy2 / a;
  elight1 = e * alight1 / a;
  elight2 = e * alight2 / a;

  epsilon0_1_saddle = eld - e_zero_point - epot0_mode1_saddle;
  //            /* excitation energy at saddle mode 1 without shell effect */
  epsilon0_2_saddle = eld - e_zero_point - epot0_mode2_saddle;
  //            /* excitation energy at saddle mode 2 without shell effect */

  epsilon_1_saddle = eld - e_zero_point - epot_mode1_saddle;
  //            /* excitation energy at saddle mode 1 with shell effect */
  epsilon_2_saddle = eld - e_zero_point - epot_mode2_saddle;
  //            /* excitation energy at saddle mode 2 with shell effect */
  epsilon_symm_saddle = eld - e_zero_point - epot_symm_saddle;

  //  /* global parameters */
  eexc1_saddle = epsilon_1_saddle;
  eexc2_saddle = epsilon_2_saddle;
  eexc_max = max(eexc1_saddle,eexc2_saddle);
  eexc_max = max(eexc_max,eld);
       
  //         /* EEXC_MAX is energy above the lowest saddle */


  epsilon0_1_scission = eld + e_saddle_scission - epot0_mode1_saddle;
  //                    /* excitation energy without shell effect */
  epsilon0_2_scission = eld + e_saddle_scission - epot0_mode2_saddle;
  //                    /* excitation energy without shell effect */

  epsilon_1_scission = eld + e_saddle_scission - epot_mode1_saddle;
  //                    /* excitation energy at scission */
  epsilon_2_scission = eld+ e_saddle_scission - epot_mode2_saddle;
  //                    /* excitation energy at scission */
  epsilon_symm_scission = eld + e_saddle_scission - epot_symm_saddle;
  //           /* excitation energy of symmetric fragment at scission */

  //     /* Calculate widhts at the saddle:                */

  e_eff1_saddle = epsilon0_1_saddle - delta_u1 * (std::exp((-epsilon_1_saddle*gamma)));
   
  if (e_eff1_saddle > 0.0) {
    wzasymm1_saddle = std::sqrt( (0.5 * 
			     (std::sqrt(1.0/a_levdens*e_eff1_saddle)) /
			     (cz_asymm1_shell * std::exp((-epsilon_1_saddle*gamma)) + cz_symm) ) );
  } 
  else {
    wzasymm1_saddle = 1.0;
  }

  e_eff2_saddle = epsilon0_2_saddle - delta_u2 * (std::exp((-epsilon_2_saddle*gamma)));
  if (e_eff2_saddle > 0.0) {
    wzasymm2_saddle = std::sqrt( (0.5 * 
			     (std::sqrt(1.0/a_levdens*e_eff2_saddle)) /
			     (cz_asymm2_shell * std::exp((-epsilon_2_saddle*gamma)) + cz_symm) ) );
  } 
  else {
    wzasymm2_saddle = 1.0;
  }

  if (eld > e_zero_point) {
    if ( (eld + epsilon_symm_saddle) < 0.0)  {
      G4cout << "<e> eld + epsilon_symm_saddle < 0" << G4endl;
    }
    wzsymm_saddle = std::sqrt( (0.5 * 
			   (std::sqrt(1.0/a_levdens*(eld+epsilon_symm_saddle))) / cz_symm ) );
  } else {
    wzsymm_saddle = 1.0;
  }

  if (itest == 1) {
    G4cout << "wz1(saddle) = " << wzasymm1_saddle << G4endl;
    G4cout << "wz2(saddle) = " << wzasymm2_saddle << G4endl;
    G4cout << "wzsymm(saddle) = " << wzsymm_saddle << G4endl;
  }
     
  //     /* Calculate widhts at the scission point: */
  //     /* fits of ref. Beizin 1991 (Plots brought to GSI by Sergei Zhdanov) */

  wzsymm_scission = wzsymm_saddle;

  if (e_saddle_scission == 0.0) {

    wzasymm1_scission = wzasymm1_saddle;
    wzasymm2_scission = wzasymm2_saddle;

  } 
  else {

    if (nheavy1_eff > 75.0) {
      wzasymm1_scission = (std::sqrt(21.0)) * z/a;
      wzasymm2_scission = (std::sqrt (max( (70.0-28.0)/3.0*(z*z/a-35.0)+28.,0.0 )) ) * z/a;
    } 
    else {
      wzasymm1_scission = wzasymm1_saddle;
      wzasymm2_scission = wzasymm2_saddle;
    }

  }

  wzasymm1_scission = max(wzasymm1_scission,wzasymm1_saddle);
  wzasymm2_scission = max(wzasymm2_scission,wzasymm2_saddle);

  wzasymm1 = wzasymm1_scission * fwidth_asymm1;
  wzasymm2 = wzasymm2_scission * fwidth_asymm2;
  wzsymm = wzsymm_scission;

  /*      if (ITEST == 1) {
	  G4cout << "WZ1(scission) = " << WZasymm1_scission << G4endl;
	  G4cout << "WZ2(scission) = " << WZasymm2_scission << G4endl;
	  G4cout << "WZsymm(scission) = " << WZsymm_scission << G4endl;
	  }
	  if (ITEST == 1) {
	  G4cout << "WZ1(scission) final= " << WZasymm1 << G4endl;
	  G4cout << "WZ2(scission) final= " << WZasymm2 << G4endl;
	  G4cout << "WZsymm(scission) final= " << WZsymm << G4endl;
	  } */
      
  wasymm = wzsymm * a/z;
  waheavy1 = wzasymm1 * a/z;
  waheavy2 = wzasymm2 * a/z;

  wasymm_saddle = wzsymm_saddle * a/z;
  waheavy1_saddle = wzasymm1_saddle * a/z;
  waheavy2_saddle = wzasymm2_saddle * a/z;

  if (itest == 1) {
    G4cout << "wasymm = " << wzsymm << G4endl;
    G4cout << "waheavy1 = " << waheavy1 << G4endl;
    G4cout << "waheavy2 = " << waheavy2 << G4endl;
  }
  // Up to here: Ok! Checked CS 11/10/05
            
  if ( (epsilon0_1_saddle - delta_u1*std::exp((-epsilon_1_saddle*gamma_heavy1))) < 0.0) {
    sasymm1 = -10.0;
  } 
  else {
    sasymm1 = 2.0 * std::sqrt( a_levdens * (epsilon0_1_saddle - 
				       delta_u1*(std::exp((-epsilon_1_saddle*gamma_heavy1))) ) );
  }

  if ( (epsilon0_2_saddle - delta_u2*std::exp((-epsilon_2_saddle*gamma_heavy2))) < 0.0) {
    sasymm2 = -10.0;
  } 
  else {
    sasymm2 = 2.0 * std::sqrt( a_levdens * (epsilon0_2_saddle - 
				       delta_u2*(std::exp((-epsilon_2_saddle*gamma_heavy2))) ) );
  }
              
  if (epsilon_symm_saddle > 0.0) {
    ssymm = 2.0 * std::sqrt( a_levdens*(epsilon_symm_saddle) );
  } 
  else {
    ssymm = -10.0;
  }
      
  if (ssymm > -10.0) {
    ysymm = 1.0;

    if (epsilon0_1_saddle < 0.0) {
      //  /* low energy */
      yasymm1 = std::exp((sasymm1-ssymm)) * wzasymm1_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } 
    else {
      //        /* high energy */
      ssymm_mode1 = 2.0 * std::sqrt( a_levdens*(epsilon0_1_saddle) );
      yasymm1 = ( std::exp((sasymm1-ssymm)) - std::exp((ssymm_mode1 - ssymm)) )  
	* wzasymm1_saddle / wzsymm_saddle * 2.0;
    }

    if (epsilon0_2_saddle < 0.0) {
      //  /* low energy */
      yasymm2 = std::exp((sasymm2-ssymm)) * wzasymm2_saddle / wzsymm_saddle * 2.0;
      //           /* factor of 2 for symmetry classes */
    } 
    else {
      //        /* high energy */
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
    yasymm = yasymm1 + yasymm2;
  } 
  else {
    ysymm = 0.0;
    yasymm1 = 0.0;
    yasymm2 = 0.0;
    //        /* search minimum threshold and attribute all events to this mode */
    if ( (epsilon_symm_saddle < epsilon_1_saddle) && (epsilon_symm_saddle < epsilon_2_saddle) ) {
      ysymm = 1.0;
    } 
    else {
      if (epsilon_1_saddle < epsilon_2_saddle) {
	yasymm1 = 1.0;
      } 
      else {
	yasymm2 = 1.0;
      }
    }
  }

  if (itest == 1) {
    G4cout << "ysymm normalized= " << ysymm  << G4endl;
    G4cout << "yasymm1 normalized= " << yasymm1  << G4endl;
    G4cout << "yasymm2 normalized= " << yasymm2  << G4endl;
  }
  // Up to here: Ok! Ckecked CS 11/10/05      
      
  //      /* even-odd effect */
  //      /* simple parametrization KHS, Nov. 2000. From Rejmund et al. */
  if ((int)(z) % 2 == 0) {
    r_e_o_exp = -0.017 * (e_saddle_scission + eld) * (e_saddle_scission + eld);
    if ( r_e_o_exp < -307.0) {
      r_e_o_exp = -307.0;
      r_e_o = std::pow(10.0,r_e_o_exp);
    }
    else {
      r_e_o = std::pow(10.0,r_e_o_exp);
    }
  } 
  else {
    r_e_o = 0.0;
  }
    
  //      $LOOP;    /* event loop */
  //     I_COUNT = I_COUNT + 1;

  //     /* random decision: symmetric or asymmetric */
  //     /* IMODE = 1 means asymmetric fission, mode 1,
  //        IMODE = 2 means asymmetric fission, mode 2,
  //        IMODE = 3 means symmetric  */
  //      RMODE = dble(HAZ(k));
  //      rmode = rnd.rndm();  

  // Safety check added to make sure we always select well defined
  // fission mode.
  do {
    rmode = haz(k);
    // Cast for test CS 11/10/05
    //      RMODE = 0.54;    
    //  rmode = 0.54;
    if (rmode < yasymm1) {
      imode = 1;
    } else if ( (rmode > yasymm1) && (rmode < (yasymm1+yasymm2)) ) {
      imode = 2;
    } else if ( (rmode > yasymm1) && (rmode > (yasymm1+yasymm2)) ) {
      imode = 3;
    }
  } while(imode == 0);

  //     /* determine parameters of the Z distribution */
  // force imode (for testing, PK)
  // imode = 3;
  if (imode == 1) {
    z1mean = zheavy1;
    z1width = wzasymm1;
  }
  if (imode == 2) {
    z1mean = zheavy2;
    z1width = wzasymm2;
  }
  if (imode == 3) {
    z1mean = zsymm;
    z1width = wzsymm;
  }

  if (itest == 1) {
    G4cout << "nbre aleatoire tire " << rmode << G4endl;
    G4cout << "fission mode " << imode << G4endl;
    G4cout << "z1mean= " << z1mean << G4endl;
    G4cout << "z1width= " << z1width << G4endl;
  }
		      
  //     /* random decision: Z1 and Z2 at scission: */
  z1 = 1.0;
  z2 = 1.0;
  while  ( (z1<5.0) || (z2<5.0) ) {
    //         Z1 = dble(GAUSSHAZ(K,sngl(Z1mean),sngl(Z1width)));
    //	 z1 = rnd.gaus(z1mean,z1width);
    z1 = gausshaz(k, z1mean, z1width);
    z2 = z - z1;
  }
  if (itest == 1) {
    G4cout << "ff charge sample " << G4endl;
    G4cout << "z1 =  " << z1 << G4endl;
    G4cout << "z2 = " << z2 << G4endl;
  }

  //     CALL EVEN_ODD(Z1,R_E_O,I_HELP);
  //         /* Integer proton number with even-odd effect */
  //     Z1 = REAL(I_HELP)
  //      /* Z1 = INT(Z1+0.5E0); */
  z2 = z - z1;

  //     /* average N of both fragments: */
  if (imode == 1) {
    n1mean = (z1 + cpol1 * a/n) * n/z;
  }
  if (imode == 2) { 
    n1mean = (z1 + cpol2 * a/n) * n/z;
  }
  /*       CASE(99)   ! only for testing;
	   N1UCD = Z1 * N/Z;
	   N2UCD = Z2 * N/Z;
	   re1 = UMASS(Z1,N1UCD,0.6) +;
	   &         UMASS(Z2,N2UCD,0.6) +;
	   &         ECOUL(Z1,N1UCD,0.6,Z2,N2UCD,0.6,d);
	   re2 = UMASS(Z1,N1UCD+1.,0.6) +;
	   &         UMASS(Z2,N2UCD-1.,0.6) +;
	   &         ECOUL(Z1,N1UCD+1.,0.6,Z2,N2UCD-1.,0.6,d);
	   re3 = UMASS(Z1,N1UCD+2.,0.6) +;
	   &         UMASS(Z2,N2UCD-2.,0.6) +;
	   &         ECOUL(Z1,N1UCD+2.,0.6,Z2,N2UCD-2.,0.6,d);
	   eps2 = (re1-2.0*re2+re3) / 2.0;
	   eps1 = re2 - re1 - eps2;
	   DN1_POL = - eps1 / (2.0 * eps2);
	   N1mean = N1UCD + DN1_POL; */
  if (imode == 3) {
    n1ucd = z1 * n/z;
    n2ucd = z2 * n/z;
    re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) + ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,d);
    re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) + ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,d);
    re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) + ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,d);
    eps2 = (re1-2.0*re2+re3) / 2.0;
    eps1 = re2 - re1 - eps2;
    dn1_pol = - eps1 / (2.0 * eps2);
    n1mean = n1ucd + dn1_pol;
  }
  // all fission modes features have been checked CS 11/10/05
  n2mean = n - n1mean;
  z2mean = z - z1mean;
      
  //   /* Excitation energies */
  //     /* formulated in energies in close consistency with the fission model */

  //     /* E_defo = UMASS(Z*0.5E0,N*0.5E0,0.6E0) -
  //                 UMASS(Z*0.5E0,N*0.5E0,0);   */
  //       /* calculates the deformation energy of the liquid drop for
  //          deformation beta = 0.6 which is most probable at scission */

  //    /* N1R and N2R provisionaly taken without fluctuations in
  //       polarisation:                                              */
  n1r = n1mean;
  n2r = n2mean;
  a1r = n1r + z1;
  a2r = n2r + z2;

  if (imode == 1) { /* N = 82 */;
    //! /* Eexc at scission */
    e_scission = max(epsilon_1_scission,1.0);
    if (n1mean > (n * 0.5) ) {
      //! /* 1. fragment is spherical */
      beta1 = 0.0;
      beta2 = 0.6;
      e1exc = epsilon_1_scission * a1r / a;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_1_scission * a2r / a + e_defo;
    }
    else {
      //!                       /* 2. fragment is spherical */
      beta1 = 0.6;
      beta2 = 0.0;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_1_scission * a1r / a + e_defo;
      e2exc = epsilon_1_scission * a2r / a;
    }
  }
	   
  if (imode == 2) { 
    //! /* N appr. 86 */
    e_scission = max(epsilon_2_scission,1.0);
    if (n1mean >  (n * 0.5) ) { 
      //! /* 2. fragment is spherical */
      beta1 = (n1r - nheavy2) * 0.034 + 0.3;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_2_scission * a1r / a + e_defo;
      beta2 = 0.6 - beta1;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_2_scission * a2r / a + e_defo;
    }
    else {
      //!                      /* 1. fragment is spherical */
      beta2 = (n2r - nheavy2) * 0.034 + 0.3;
      e_defo = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
      e2exc = epsilon_2_scission * a2r / a + e_defo;
      beta1 = 0.6 - beta2;
      e_defo = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
      e1exc = epsilon_2_scission * a1r / a + e_defo;
    }
  }
         
  if (imode == 3) { 
    // ! /* Symmetric fission channel */

    //             /* the fit function for beta is the deformation for
    //                optimum energy at the scission point, d = 2 */
    //             /* beta  : deformation of symmetric fragments */
    //             /* beta1 : deformation of first fragment */
    //             /* beta2 : deformation of second fragment */
    beta =  0.177963 + 0.0153241 * zsymm - 0.000162037 * zsymm*zsymm;
    beta1 = 0.177963 + 0.0153241 * z1 - 0.000162037 * z1*z1;
    //            beta1 = 0.6
    e_defo1 = umass(z1,n1r,beta1) - umass(z1,n1r,0.0);
    beta2 = 0.177963 + 0.0153241 * z2 - 0.000162037 * z2*z2;
    //            beta2 = 0.6
    e_defo2 = umass(z2,n2r,beta2) - umass(z2,n2r,0.0);
    e_asym = umass(z1 , n1r, beta1) + umass(z2, n2r ,beta2)
      + ecoul(z1,n1r,beta1,z2,n2r,beta2,2.0)
      - 2.0 * umass(zsymm,nsymm,beta)
      - ecoul(zsymm,nsymm,beta,zsymm,nsymm,beta,2.0);
    //            E_asym = CZ_symm * (Z1 - Zsymm)**2
    e_scission = max((epsilon_symm_scission - e_asym),1.0);
    //         /*  $LIST(Z1,N1R,Z2,N2R,E_asym,E_scission); */
    e1exc = e_scission * a1r / a + e_defo1;
    e2exc = e_scission * a2r / a + e_defo2;
  }
  // Energies checked for all the modes CS 11/10/05
	  
  //    /* random decision: N1R and N2R at scission, before evaporation: */
  //    /* CN =       UMASS(Zsymm , Nsymm + 1.E0,0) +
  //                UMASS(Zsymm, Nsymm - 1.E0,0)
  //                 + 1.44E0 * (Zsymm)**2 /
  //                 (r_null**2 * ((Asymm+1)**1/3 + (Asymm-1)**1/3)**2 )
  //                 - 2.E0 * UMASS(Zsymm,Nsymm,0)
  //                 - 1.44E0 * (Zsymm)**2 / (r_null * 2.E0 * (Asymm)**1/3)**2; */


  //    /* N1width = std::sqrt(0.5E0 * std::sqrt(1.E0/A_levdens*(Eld+E_saddle_scission)) / CN); */
  //    /* 8. 9. 1998: KHS (see also consideration in the first comment block)
  //       sigma_N(Z=const) = A/Z  * sigma_Z(A=const)
  //       sigma_Z(A=const) = 0.4 to 0.5  (from Lang paper Nucl Phys. A345 (1980) 34)
  //       sigma_N(Z=const) = 0.45 * A/Z  (= 1.16 for 238U)
  //       therefore: SIGZMIN = 1.16                                              */

  if ( (imode == 1) || (imode == 2) ) {
    cn=(umass(z1,n1mean+1.,beta1) + umass(z1,n1mean-1.,beta1)
	+ umass(z2,n2mean+1.,beta2) + umass(z2,n2mean-1.,beta2)
	+ ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
	+ ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
	- 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	- 2.0 * umass(z1, n1mean, beta1)
	- 2.0 * umass(z2, n2mean, beta2) ) * 0.5;
    //          /* Coulomb energy neglected for the moment! */
    //           IF (E_scission.lt.0.) Then
    //             write(6,*)'<E> E_scission < 0, MODE 1,2'
    //           ENDIF
    //           IF (CN.lt.0.) Then
    //             write(6,*)'CN < 0, MODE 1,2'
    //           ENDIF
    n1width=std::sqrt( (0.5 * (std::sqrt(1.0/a_levdens*(e_scission)))/cn) );
    n1width=max(n1width, sigzmin);

    //          /* random decision: N1R and N2R at scission, before evaporation: */
    n1r = 1.0;
    n2r = 1.0;
    while  ( (n1r<5.0) || (n2r<5.0) ) {
      //             n1r = dble(gausshaz(k,sngl(n1mean),sngl(n1width)));
      //	   	n1r = rnd.gaus(n1mean,n1width);
      n1r = gausshaz(k, n1mean, n1width);
      n2r = n - n1r;
    }	   
    //           N1R = GAUSSHAZ(K,N1mean,N1width)
    if (itest == 1) {
      G4cout << "after neutron sample " << n1r << G4endl;
    }	   
    n1r = (float)( (int)((n1r+0.5)) );
    n2r = n - n1r;
 
    even_odd(z1,r_e_o,i_help);
    //         /*  proton number with even-odd effect */
    z1 = (float)(i_help);
    z2 = z - z1;

    a1r = z1 + n1r;
    a2r = z2 + n2r;
  }
	
  if (imode == 3) {
    //!  /* When(3) */
    if (nzpol > 0.0) {
      //           /* We treat a simultaneous split in Z and N to determine polarisation */
      cz = ( umass(z1-1., n1mean+1.,beta1)
	     +  umass(z2+1., n2mean-1.,beta1) 
	     +  umass(z1+1., n1mean-1.,beta2)
	     +  umass(z2 - 1., n2mean + 1.,beta2)
	     +  ecoul(z1-1.,n1mean+1.,beta1,z2+1.,n2mean-1.,beta2,2.0)
	     +  ecoul(z1+1.,n1mean-1.,beta1,z2-1.,n2mean+1.,beta2,2.0)
	     -  2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	     -  2.0 * umass(z1, n1mean,beta1)
	     -  2.0 * umass(z2, n2mean,beta2) ) * 0.5;
      //           IF (E_scission.lt.0.) Then
      //             write(6,*) '<E> E_scission < 0, MODE 1,2'
      //           ENDIF
      //           IF (CZ.lt.0.) Then
      //             write(6,*) 'CZ < 0, MODE 1,2'
      //           ENDIF
      za1width=std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*(e_scission)) / cz) );
      za1width=std::sqrt( (max((za1width*za1width-(1.0/12.0)),0.1)) );
      //                        /* Check the value of 0.1 ! */
      //                        /* Shephard correction */
      a1r = z1 + n1mean;
      a1r = (float)((int)((a1r+0.5)));
      a2r = a - a1r;
      //           /* A1R and A2R are integer numbers now */
      //        /* $LIST(A1R,A2R,ZA1WIDTH); */

      n1ucd = n/a * a1r;
      n2ucd = n/a * a2r;
      z1ucd = z/a * a1r;
      z2ucd = z/a * a2r;

      re1 = umass(z1ucd-1.,n1ucd+1.,beta1) + umass(z2ucd+1.,n2ucd-1.,beta2)
	+ ecoul(z1ucd-1.,n1ucd+1.,beta1,z2ucd+1.,n2ucd-1.,beta2,d);
      re2 = umass(z1ucd,n1ucd,beta1) + umass(z2ucd,n2ucd,beta2)
	+ ecoul(z1ucd,n1ucd,beta1,z2ucd,n2ucd,beta2,d);
      re3 = umass(z1ucd+1.,n1ucd-1.,beta1) + umass(z2ucd-1.,n2ucd+1.,beta2) +
	+ ecoul(z1ucd+1.,n1ucd-1.,beta1,z2ucd-1.,n2ucd+1.,beta2,d);
		 
      eps2 = (re1-2.0*re2+re3) / 2.0;
      eps1 = (re3 - re1)/2.0;
      dn1_pol = - eps1 / (2.0 * eps2);
      z1 = z1ucd + dn1_pol;
      if (itest == 1) {
	G4cout << "before proton sample " << z1 << G4endl;
      }  
      //           Z1 = dble(GAUSSHAZ(k,sngl(Z1),sngl(ZA1width)));
      //	   z1 = rnd.gaus(z1,za1width);
      z1 = gausshaz(k, z1, za1width);
      if (itest == 1) {
	G4cout << "after proton sample " << z1 << G4endl;
      }	   
      even_odd(z1,r_e_o,i_help);
      //         /* proton number with even-odd effect */
      z1 = (float)(i_help);
      z2 = (float)((int)( (z - z1 + 0.5)) );

      n1r = a1r - z1;
      n2r = n - n1r;
    } 
    else {
      //           /* First division of protons, then adjustment of neutrons */
      cn = ( umass(z1, n1mean+1.,beta1) + umass(z1, n1mean-1., beta1)
	     + umass(z2, n2mean+1.,beta2) + umass(z2, n2mean-1., beta2)
	     + ecoul(z1,n1mean+1.,beta1,z2,n2mean-1.,beta2,2.0)
	     + ecoul(z1,n1mean-1.,beta1,z2,n2mean+1.,beta2,2.0)
	     - 2.0 * ecoul(z1,n1mean,beta1,z2,n2mean,beta2,2.0)
	     - 2.0 * umass(z1, n1mean, 0.6)
	     - 2.0 * umass(z2, n2mean, 0.6) ) * 0.5;
      //          /* Coulomb energy neglected for the moment! */
      //           IF (E_scission.lt.0.) Then
      //             write(6,*) '<E> E_scission < 0, MODE 1,2'
      //           Endif
      //           IF (CN.lt.0.) Then
      //             write(6,*) 'CN < 0, MODE 1,2'
      //           Endif
      n1width=std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*(e_scission)) / cn) );
      n1width=max(n1width, sigzmin);

      //          /* random decision: N1R and N2R at scission, before evaporation: */
      //           N1R = dble(GAUSSHAZ(k,sngl(N1mean),sngl(N1width)));
      //	   n1r = rnd.gaus(n1mean,n1width);
      n1r = gausshaz(k, n1mean, n1width);
      n1r = (float)( (int)((n1r+0.5)) );
      n2r = n - n1r;

      even_odd(z1,r_e_o,i_help);
      //         /* Integer proton number with even-odd effect */
      z1 = (float)(i_help);
      z2 = z - z1;

      a1r = z1 + n1r;
      a2r = z2 + n2r;
          
    }
  }

  if (itest == 1) {
    G4cout << "remid imode = " << imode << G4endl;
    G4cout << "n1width =  " << n1width << G4endl;
    G4cout << "n1r = " << n1r << G4endl;
    G4cout << "a1r = " << a1r << G4endl;
    G4cout << "n2r = " << n2r << G4endl;
    G4cout << "a2r = " << a2r << G4endl;
  }
  // Up to here: checked CS 11/10/05	
	
  //      /* Extracted from Lang et al. Nucl. Phys. A 345 (1980) 34 */
  e1exc_sigma = 5.5;
  e2exc_sigma = 5.5;

 neufcentquatrevingtsept:
  //       E1final = dble(Gausshaz(k,sngl(E1exc),sngl(E1exc_sigma)));
  //       E2final = dble(Gausshaz(k,sngl(E2exc),sngl(E2exc_sigma)));
  //        e1final = rnd.gaus(e1exc,e1exc_sigma);
  //        e2final = rnd.gaus(e2exc,e2exc_sigma);
  e1final = gausshaz(k, e1exc, e1exc_sigma);
  e2final = gausshaz(k, e2exc, e2exc_sigma);
  if ( (e1final < 0.0) || (e2final < 0.0) ) goto neufcentquatrevingtsept;
  if (itest == 1) {
    G4cout << "sampled exc 1 " << e1final << G4endl;
    G4cout << "sampled exc 2 " << e2final << G4endl;
  }

  //      /* OUTPUT QUANTITIES OF THE EVENT GENERATOR:        */

  //      /* Quantities before neutron evaporation            */

  //      /* Neutron number of prefragments: N1R and N2R      */
  //      /* Atomic number of fragments: Z1 and Z2            */
  //      /* Kinetic energy of fragments: EkinR1, EkinR2      *7

  //      /* Quantities after neutron evaporation:            */

  //      /* Neutron number of fragments: N1 and N2           */
  //      /* Mass number of fragments: A1 and A2              */
  //      /* Atomic number of fragments: Z1 and Z2            */
  //      /* Number of evaporated neutrons: N1R-N1 and N2R-N2 */
  //      /* Kinetic energy of fragments: EkinR1*A1/A1R and
  //                                      EkinR2*A2/A2R       */

  n1 = n1r;
  n2 = n2r;
  a1 = n1 + z1;
  a2 = n2 + z2;
  e1 = e1final;
  e2 = e2final;

  //       /* Pre-neutron-emission total kinetic energy: */
  tker = (z1 * z2 * 1.44) /
    ( r0 * std::pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
      r0 * std::pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
  //       /* Pre-neutron-emission kinetic energy of 1. fragment: */
  ekinr1 = tker * a2 / a;
  //       /* Pre-neutron-emission kinetic energy of 2. fragment: */
  ekinr2 = tker * a1 / a;

  v1 = std::sqrt( (ekinr1/a1) ) * 1.3887;
  v2 = std::sqrt( (ekinr2/a2) ) * 1.3887;

  if (itest == 1) {
    G4cout << "ekinr1 " << ekinr1 << G4endl;
    G4cout << "ekinr2 " << ekinr2 << G4endl;
  }

 milledeux:       
  //**************************
  //*** only symmetric fission
  //**************************
  // Symmetric fission: Ok! Checked CS 10/10/05
  if ( (icz == -1) || (a1 < 0.0) || (a2 < 0.0) ) {
    //           IF (z.eq.92) THEN
    //              write(6,*)'symmetric fission'
    //              write(6,*)'Z,A,E,A1,A2,icz,Atot',Z,A,E,A1,A2,icz,Atot
    //           END IF

    if (itest == 1) {
      G4cout << "milledeux: liquid-drop option "  << G4endl;
    }

    n = a-z;
    //  proton number in symmetric fission (centre) *
    zsymm  = z / 2.0;
    nsymm  = n / 2.0;
    asymm = nsymm + zsymm;

    a_levdens = a / xlevdens;

    masscurv = 2.0;
    cz_symm = 8.0 / std::pow(z,2) * masscurv;

    wzsymm = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cz_symm) ) ;

    if (itest == 1) {
      G4cout << " symmetric high energy fission " << G4endl;
      G4cout << "wzsymm " << wzsymm << G4endl;
    }

    z1mean = zsymm;
    z1width = wzsymm;

    // random decision: Z1 and Z2 at scission: */
    z1 = 1.0;
    z2 = 1.0;
    while  ( (z1 < 5.0) || (z2 < 5.0) ) {
      //           z1 = dble(gausshaz(kkk,sngl(z1mean),sngl(z1width)));
      //	   z1 = rnd.gaus(z1mean,z1width);
      z1 = gausshaz(kkk, z1mean, z1width);
      z2 = z - z1;
    }

    if (itest == 1) {
      G4cout << " z1 " << z1 << G4endl;
      G4cout << " z2 " << z2 << G4endl;
    }
    if (itest == 1) {
      G4cout << " zsymm " << zsymm << G4endl;
      G4cout << " nsymm " << nsymm << G4endl;
      G4cout << " asymm " << asymm << G4endl;
    }
    //    	CN =  UMASS(Zsymm , Nsymm + 1.E0) + UMASS(Zsymm, Nsymm - 1.E0)
    //    #            + 1.44E0 * (Zsymm)**2 /
    //    #            (r_null**2 * ((Asymm+1)**(1./3.) +
    //    #            (Asymm-1)**(1./3.))**2 )
    //    #            - 2.E0 * UMASS(Zsymm,Nsymm)
    //    #            - 1.44E0 * (Zsymm)**2 /
    //    #            (r_null * 2.E0 * (Asymm)**(1./3.))**2

    n1ucd = z1 * n/z;
    n2ucd = z2 * n/z;
    re1 = umass(z1,n1ucd,0.6) + umass(z2,n2ucd,0.6) +
      ecoul(z1,n1ucd,0.6,z2,n2ucd,0.6,2.0);
    re2 = umass(z1,n1ucd+1.,0.6) + umass(z2,n2ucd-1.,0.6) +
      ecoul(z1,n1ucd+1.,0.6,z2,n2ucd-1.,0.6,2.0);
    re3 = umass(z1,n1ucd+2.,0.6) + umass(z2,n2ucd-2.,0.6) +
      ecoul(z1,n1ucd+2.,0.6,z2,n2ucd-2.,0.6,2.0);
    reps2 = (re1-2.0*re2+re3)/2.0;
    reps1 = re2 - re1 -reps2;
    rn1_pol = -reps1/(2.0*reps2);
    n1mean = n1ucd + rn1_pol;
    n2mean = n - n1mean;

    if (itest == 1) {
      G4cout << " n1mean " << n1mean << G4endl;
      G4cout << " n2mean " << n2mean << G4endl;
    }

    cn = (umass(z1,n1mean+1.,0.0) + umass(z1,n1mean-1.,0.0) +
	  + umass(z2,n2mean+1.,0.0) + umass(z2,n2mean-1.,0.0)
	  - 2.0 * umass(z1,n1mean,0.0) +
	  - 2.0 * umass(z2,n2mean,0.0) ) * 0.5;
    //      This is an approximation! Coulomb energy is neglected.

    n1width = std::sqrt( (0.5 * std::sqrt(1.0/a_levdens*e) / cn) );

    if (itest == 1) {
      G4cout << " cn " << cn << G4endl;
      G4cout << " n1width " << n1width << G4endl;
    }
	
    // random decision: N1R and N2R at scission, before evaporation: */
    //       N1R = dfloat(NINT(GAUSSHAZ(KKK,sngl(N1mean),sngl(N1width))));
    //      n1r = (float)( (int)(rnd.gaus(n1mean,n1width)) );
    n1r = (float)( (int)(gausshaz(k, n1mean,n1width)) );
    n2r = n - n1r;
    // Mass of first and second fragment */
    a1 = z1 + n1r;
    a2 = z2 + n2r;

    e1 = e*a1/(a1+a2);
    e2 = e - e*a1/(a1+a2);
    if (itest == 1) {
      G4cout << " n1r " << n1r << G4endl;
      G4cout << " n2r " << n2r << G4endl;
    }

  }

  if (itest == 1) {
    G4cout << " a1 " << a1 << G4endl;
    G4cout << " z1 " << z1 << G4endl;
    G4cout << " a2 " << a2 << G4endl;
    G4cout << " z2 " << z2 << G4endl;
    G4cout << " e1 " << e1 << G4endl;
    G4cout << " e2 " << e << G4endl;
  }

  //       /* Pre-neutron-emission total kinetic energy: */
  tker = (z1 * z2 * 1.44) /
    ( r0 * std::pow(a1,0.33333) * (1.0 + 2.0/3.0 * beta1) +
      r0 * std::pow(a2,0.33333) * (1.0 + 2.0/3.0 * beta2) + 2.0 );
  //       /* Pre-neutron-emission kinetic energy of 1. fragment: */
  ekin1 = tker * a2 / a;
  //       /* Pre-neutron-emission kinetic energy of 2. fragment: */
  ekin2 = tker * a1 / a;

  v1 = std::sqrt( (ekin1/a1) ) * 1.3887;
  v2 = std::sqrt( (ekin2/a2) ) * 1.3887;

  if (itest == 1) {
    G4cout << " kinetic energies " << G4endl;
    G4cout << " ekin1 " << ekin1 << G4endl;
    G4cout << " ekin2 " << ekin2 << G4endl;
  }
}

void G4AblaFissionSimfis18::standardRandom(G4double *rndm, G4long*)
{
  // Use Geant4 G4UniformRand
  (*rndm) = randomGenerator->getRandom();
}

G4double G4AblaFissionSimfis18::haz(G4int k)
{
  const G4int pSize = 110;
  static G4double p[pSize];
  static G4long ix = 0, i = 0;
  static G4double x = 0.0, y = 0.0, a = 0.0, haz = 0.0;
  //  k =< -1 on initialise                                        
  //  k = -1 c'est reproductible                                   
  //  k < -1 || k > -1 ce n'est pas reproductible

  // Zero is invalid random seed. Set proper value from our random seed collection:
  if(ix == 0) {
    ix = hazard->ial;
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
    x = randomGenerator->getRandom();
    //    standardRandom(&x, &ix);
    for(G4int i = 0; i < pSize; i++) { //do i=1,110                                                 
      p[i] = randomGenerator->getRandom();
      //      standardRandom(&(p[i]), &ix);
    }
    a = randomGenerator->getRandom();
    standardRandom(&a, &ix);
    k = 0;
  }

  i = nint(100*a)+1;
  haz = p[i];
  a = randomGenerator->getRandom();
  //  standardRandom(&a, &ix);
  p[i] = a;

  hazard->ial = ix;
  return haz;
}


G4double G4AblaFissionSimfis18::gausshaz(int k, double xmoy, double sig)
{
  // Gaussian random numbers:

  //   1005       C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
  static G4int  iset = 0;
  static G4double v1,v2,r,fac,gset,gausshaz;

  if(iset == 0) { //then                                              
    do {
      v1 = 2.0*haz(k) - 1.0;
      v2 = 2.0*haz(k) - 1.0;
      r = std::pow(v1,2) + std::pow(v2,2);
    } while(r >= 1);

    fac = std::sqrt(-2.*std::log(r)/r);
    gset = v1*fac;
    gausshaz = v2*fac*sig+xmoy;
    iset = 1;
  }
  else {
    gausshaz=gset*sig+xmoy;
    iset=0;
  }
  return gausshaz;                                                         
}


// Utilities

G4double G4AblaFissionSimfis18::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFissionSimfis18::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4AblaFissionSimfis18::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFissionSimfis18::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4AblaFissionSimfis18::nint(G4double number)
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

G4int G4AblaFissionSimfis18::secnds(G4int x)
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

G4int G4AblaFissionSimfis18::mod(G4int a, G4int b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0;
  } 
}

G4double G4AblaFissionSimfis18::dmod(G4double a, G4double b)
{
  if(b != 0) {
    return (a - (a/b)*b);
  }
  else {
    return 0.0;
  } 
}

G4double G4AblaFissionSimfis18::dint(G4double a)
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

G4int G4AblaFissionSimfis18::idint(G4double a)
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

G4int G4AblaFissionSimfis18::idnint(G4double value)
{
  G4double valueCeil = int(std::ceil(value));
  G4double valueFloor = int(std::floor(value));

  if(std::fabs(value - valueCeil) < std::fabs(value - valueFloor)) {
    return int(valueCeil);
  }
  else {
    return int(valueFloor);
  }
}

G4double G4AblaFissionSimfis18::dmin1(G4double a, G4double b, G4double c)
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

G4double G4AblaFissionSimfis18::utilabs(G4double a)
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
