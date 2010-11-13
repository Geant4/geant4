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
// $Id: G4Incl.hh,v 1.19 2010-11-13 00:08:36 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#ifndef G4Incl_hh
#define G4Incl_hh 1

#include "globals.hh"
#include "G4InclRandomNumbers.hh"
#include "G4InclDataDefs.hh"
#include "G4Abla.hh"
#include <fstream>
#include "G4VInclLogger.hh"
#include "G4InclInput.hh"

using namespace std;
/**
 *  Class containing INCL4 hadronic cascade algorithm.
 */

class G4Incl {

public:
  /**
   *   
   * Support for Doxygen JavaDoc style.
   *
   * \author{pekka.kaitaniemi@helsinki.fi}
   */

  /**
   * Constructor to be used with Geant4.
   * @param hazard a pointer to G4Hazard structure.
   * @param calincl a pointer to G4Calincl structure.
   * @param ws a pointer to G4Ws structure.   
   * @param mat a pointer to G4Mat structure.
   * @param varntp a pointer to G4VarNtp structure.   
   */
  G4Incl(G4Hazard *hazard, G4InclInput *calincl, G4Ws *ws, G4Mat *mat, G4VarNtp *varntp);

  /**
   * Constructor for private unit testing purposes.
   * @param hazard a pointer to G4Hazard structure.
   * @param dton a pointer to G4Dton structure.
   * @param saxw a pointer to G4Saxw structure.
   * @param ws a pointer to G4Ws structure.   
   */
  G4Incl(G4Hazard *hazard, G4Dton *dton, G4Saxw *saxw, G4Ws *ws);

  G4Incl(); 
  
  ~G4Incl(); // Destructor

  void dumpParticles();
  G4double energyTest(G4int i); // Test for NaN energy of particle i.
  void dumpBl5(std::ofstream& dumpOut); // Dump the contents of G4Bl5.
  void dumpSaxw(std::ofstream& dumpOut); // Dump the contents of G4Saxw.
  void dumpBl1(std::ofstream& dumpOut); // Dump the contents of G4Bl1.
  void dumpBl2(std::ofstream& dumpOut); // Dump the contents of G4Bl2.
  void dumpBl3(std::ofstream& dumpOut); // Dump the contents of G4Bl3.

  /**
   * Set Fermi break-up use flag.
   */
  void setUseFermiBreakUp(G4bool useIt);

  /**
   * Set projectile spectator use flag.
   *
   * Select whether or not to produce so-called spectator nucleus from
   * the projectile nucleons that do not hit the target.
   */
  void setUseProjectileSpectators(G4bool useIt);

  /**
   * Set verbosity level.
   */
  void setVerboseLevel(G4int level);
  /**
   * Get verbosity level.
   */
  G4int getVerboseLevel();
  
  void setDtonData(G4Dton *newDton); // Set internal data.
  void setWsData(G4Ws *newWs); 
  void setHazardData(G4Hazard *newHazard); 
  void setSaxwData(G4Saxw *newSaxw);
  void setSpl2Data(G4Spl2 *newSpl2);
  void setMatData(G4Mat *newMat);
  void setInput(G4InclInput *newCalincl);
  void setLightNucData(G4LightNuc *newLightNuc);
  void setLightGausNucData(G4LightGausNuc *newLightGausNuc);
  void setBl1Data(G4Bl1 *newBl1);
  void setBl2Data(G4Bl2 *newBl2);
  void setBl3Data(G4Bl3 *newBl3);
  void setBl4Data(G4Bl4 *newBl4);
  void setBl5Data(G4Bl5 *newBl5);
  void setBl6Data(G4Bl6 *newBl6);
  void setBl8Data(G4Bl8 *newBl8);
  void setBl9Data(G4Bl9 *newBl9);
  void setBl10Data(G4Bl10 *newBl10);
  void setKindData(G4Kind *newKind);

  /**
   * Register the logger.
   */
  void registerLogger(G4VInclLogger *aLogger) {
    if(aLogger != 0) {
      theLogger = aLogger;
      abla->registerLogger(aLogger);
   }
  }

public:
  /**
   * Process one event with INCL4 only.
   */ 
  void processEventIncl(G4InclInput *input);

  /**
   * Process one event with INCL4 and built-in ABLA evaporation and fission.
   */ 
  void processEventInclAbla(G4InclInput *input, G4int eventnumber);

public: // Methods used to initialize INCL
  /**
   * Initialize INCL.
   *
   * @param initRandomSeed choose whether INCL should initialize random seeds.
   */
  void initIncl(G4bool initRandomSeed);

  /**
   * Initialize target materials.
   *
   * @param izmat charge number
   * @param iamat mass number
   * @param imat material number (array index)
   */
  void initMaterial(G4int izmat, G4int iamat, G4int imat);

  /**
   * A normal member taking two arguments and returning an integer value.
   *
   * @param l an integer argument.
   * @param q a constant character pointer.
   * @return The test results
   */
  G4double deutv(G4int l, G4double q);
  
  /**
   * Returns the values of the function:
   * \f[
   * (0.23162461 + (j - 1))^2
   * \f]
   *
   * @param j an integer parameter
   * @return a double value
   */
  G4double fm2(G4int j);

  /**
   * Interpolates function described by class G4Saxw around a point.
   *
   * @param xv interpolation point
   * @return a double value
   */
  G4double interpolateFunction(G4double xv);

  /**
   * Calculates the first derivative of the function stored in class G4Saxw.
   *
   * @param xv an integer parameter
   */
  void firstDerivative(G4int k);

  /**
   * Returns the values of the function:
   *  \f[
   *  \frac{r^2}{1 + e^{\frac{r - r_{0}}{A_{dif}}}}
   *  \f]
   *
   * @param r a G4double argument
   * @return a double value
   */
  G4double wsax(G4double r);

  /**
   * Returns the values of the function:
   * \f[
   *    \frac{1}{A_{dif}^{2}} \frac{r^3 e^{r - r_0}}{1 + e^{\frac{r - r_{0}}{A_{dif}}}}
   * \f]
   *
   * @param r a G4double argument
   * @return a double value
   */
  G4double derivWsax(G4double r);

  /**
   * Returns the values of the function:
   * \f[
   *    r^2 (1.0 + r_0 \frac{r^2}{A_{dif}}) e^{-\frac{r^2}{A_{dif}^2}}
   * \f]
   * @param r a G4double argument
   * @return a double value
   */
  G4double dmho(G4double r);

  /**
   * Returns the values of the function:
   * \f[
   *    -\frac{2r^4}{A_{dif}^2} (r_0(1.0 - \frac{r^2}{A_{dif}^2}) - 1.0) e^{-\frac{r^2}{A_{dif}^2}}
   * \f]
   *
   * @param r a G4double argument
   * @return a double value
   */
  G4double derivMho(G4double r);

  /**
   * Returns the values of the function:
   * \f[
   *    \frac{r^4}{A_{dif}^2} e^{-\frac{1}{2} \frac{r^2}{A_{dif}^2}}
   * \f]
   *
   * @param r a G4double argument
   * @return a double value
   */
  G4double derivGaus(G4double r);

  /**
   * Ce subroutine appele sur le premier tir va calculer la densite du deuton
   * dans l'espace des impulsions et preparer l'interpolation permettant ensuite
   * le tir au hasard d'un module de l'impulsion (q).
   * Ce subroutine remplit le G4Spl2:
   * xsp, ysp integrale normalisee de la densite de 0 a q.
   * a(),b(),c() coefs des nsp points pour une interpolation du second degre.
   * q est en fm-1. 
   */
  void densDeut();

  /**
   * Integrate using Alkazhov's method.
   *
   * @param ami a double parameter
   * @param ama a double parameter
   * @param dr a double parameter
   * @param functionChoice an integer parameter
   * @return a double value
   */
  G4double integrate(G4double ami, G4double ama, G4double step, G4int functionChoice);
    
  /**
   * Deuteron density
   *
   * @param q a double parameter
   * @return a double value
   */
  G4double dens(G4double q);

  /**
   * 
   */
  void spl2ab();

  /**
   * 
   * @param xv a double parameter
   * @return a double value
   */
  G4double splineab(G4double xv);

public: // Main INCL routines
  /**
   * INCL model as a function.
   */
  void pnu(G4int *ibert_p, G4int *nopart_p, G4int *izrem_p, G4int *iarem_p, G4double *esrem_p,
	   G4double *erecrem_p, G4double *alrem_p, G4double *berem_p, G4double *garem_p,
	   G4double *bimpact_p, G4int *l_p, G4double *xjrem, G4double *yjrem, G4double *zjrem);
     
  //C projection de JREM sur Z:
  //        MREM=ANINT(ZJREM/197.328)
  //	IF (MREM.GT.JREM) MREM=JREM
  //	IF (MREM.LT.-JREM) MREM=-JREM

  /**
   * Single nucleon-nucleon collision.
   */
  void collis(G4double *p1_p, G4double *p2_p, G4double *p3_p, G4double *e1_p, G4double *pout11_p, G4double *pout12_p, 
	      G4double *pout13_p, G4double *eout1_p, G4double *q1_p, G4double *q2_p, G4double *q3_p,
	      G4double *q4_p, G4int *np_p, G4int *ip_p, G4int *k2_p, G4int *k3_p, G4int *k4_p, 
	      G4int *k5_p, G4int *m1_p, G4int *m2_p, G4int *is1_p, G4int *is2_p);

  /**
   * This routine describes the anisotropic decay of a particle of
   * mass xi into 2 particles of masses x1,x2.
   * The anisotropy is supposed to follow a 1+3*hel*(std::cos(theta))**2 law
   * with respect to the direction of the incoming particle.
   *
   * In the input, p1,p2,p3 is the momentum of particle xi.
   *
   * In the output, p1,p2,p3 is the momentum of particle x1,
   * while q1,q2,q3 is the momentum of particle x2.
   *
   * @param p1_p pointer to momentum component 1
   * @param p2_p pointer to momentum component 2
   * @param p3_p pointer to momentum component 3
   * @param wp_p pointer to a double parameter
   * @param q1_p pointer to momentum component 1
   * @param q2_p pointer to momentum component 2
   * @param q3_p pointer to momentum component 3
   * @param wq_p pointer to a double parameter
   * @param xi_p pointer to a double parameter
   * @param x1_p pointer to momentum component 1
   * @param x2_p pointer to momentum component 2
   * @param x3_p pointer to momentum compone n3
   * @param hel_p pointer to a double parameter
   */
  void decay2(G4double *p1_p, G4double *p2_p, G4double *p3_p, G4double *wp_p, G4double *q1_p, 
	      G4double *q2_p, G4double *q3_p, G4double *wq_p, G4double *xi_p, G4double *x1_p,
	      G4double *x2_p, G4double *hel_p);

  /**
   * Time calculation.
   *
   * @param i an index of particle 1
   * @param j an index of particle 2
   */
  void time(G4int i, G4int j);

  /**
   * New time.
   *
   * @param l1 an integer parameter
   * @param l2 an integer parameter
   */
  void newt(G4int l1, G4int l2);

  /**
   * 
   *
   * @param l1 an integer parameter
   */
  void new1(G4int l1);

  /**
   * 
   *
   * @param y1 a double parameter
   * @param y2 a double parameter
   * @param y3 a double parameter
   * @param q1 a double parameter
   * @param q2 a double parameter
   * @param q3 a double parameter
   * @param q4 a double parameter
   * @param npion number of pions
   * @param l1 an integer parameter
   */
  void new2(G4double y1, G4double y2, G4double y3, G4double q1, G4double q2, G4double q3, 
	    G4double q4, G4int npion, G4int l1);

  /**
   * 
   *
   * @param y1 a double parameter
   * @param y2 a double parameter
   * @param y3 a double parameter
   * @param q1 a double parameter
   * @param q2 a double parameter
   * @param q3 a double parameter
   * @param q4 a double parameter
   * @param npion number of pions
   * @param l1 an integer parameter
   */
  void new3(G4double y1, G4double y2, G4double y3, G4double q1, G4double q2, G4double q3, 
	    G4double q4, G4int npion, G4int l1);

  /**
   * Lorentz transformation.
   *
   * @param q1 a double parameter
   * @param q2 a double parameter
   * @param q3 a double parameter
   * @param b1 a double parameter
   * @param b2 a double parameter
   * @param b3 a double parameter
   * @param E energy
   */
  void loren(G4double *q1, G4double *q2, G4double *q3, G4double *b1, G4double *b2, G4double *b3, G4double *E);

  /**
   * Pauli blocking.
   *
   * @param l an integer parameter
   * @param xr a double parameter
   * @param pr a double parameter
   * @return a double value
   */
  G4double pauliBlocking(G4int l, G4double xr, G4double pr);

  /**
   * Fit by J. Vandermeulen.
   * Low energy fit from reference J.Cugnon, D. L'hote and J. Vandermeulen, NIM B111 (1996) 215.
   *
   * @param E energy
   * @param m a double parameter m = 0, 1, 2 for nucleon-nucleon, nucleon-delta, delta-delta
   * @param i a double parameter i = 2, 0, -2  for pp, pn, nn
   * @return a double value
   */
  G4double lowEnergy(G4double E, G4double m, G4double i);

  /**
   * Total cross-sections.
   * 
   * @param E energy
   * @param m an integer parameter m=0,1,2 for nucleon-nucleon, nucleon-delta, delta-delta
   * @param i an integer parameter i = 2, 0, -2 for pp, pn, nn
   * @return a double value
   */
  G4double totalCrossSection(G4double E, G4int m, G4int i);

  /**
   *
   *
   * @param Ein energy
   * @param d a double parameter
   * @param i an integer parameter
   * @param isa an integer parameter
   * @return a double value
   */
  G4double srec(G4double Ein, G4double d, G4int i, G4int isa);

  /**
   * Delta production cross section.
   *
   * @param E energy
   * @param i an integer parameter
   * @return a double value
   */
  G4double deltaProductionCrossSection(G4double E, G4int i);

  /**
   * Sigma(pi+ + p) in the (3,3) region.
   * New fit by J. Vandermeulen and constant value above the (3,3) resonance.
   *
   * @param x a double parameter
   * @return a double value
   */
  G4double pionNucleonCrossSection(G4double x);

  /**
   * Transmission probability for a nucleon of kinetic energy
   * E on the edge of the well of depth v0 (nr approximation).
   * ,
   * of the nucleus and r is the target radius
   *
   * @param E kinetic energy
   * @param iz the isospin of the nucleon
   * @param izn the instanteneous charge
   * @param r 
   * @param v0 
   * @return a double value
   */
  //  G4double transmissionProb(G4double E, G4double iz, G4double izn, G4double r, G4double v0);
  //  G4double transmissionProb(G4double E, G4int iz, G4int izn, G4double r, G4double v0)
  G4double transmissionProb(G4double E, G4int iz, G4int ia, G4int izn, G4double r, G4double v0);
  //  G4double transmissionProb(G4double E, G4int iz, G4int ia, G4int izn,G4double R, G4double v0);
  //  G4double transmissionProb(G4double E, G4double iz, G4double izn, G4double r, G4double v0);

  void projo_spec(G4int ia1, G4int ips,
			G4double fmpinc, G4double pinc, G4double tlab);

  void ordered(G4double t, G4int nb);

  /**
   *
   *
   * @param x1 a double parameter
   * @param x2 a double parameter
   * @param x3 a double parameter
   * @param p1 a double parameter
   * @param p2 a double parameter
   * @param p3 a double parameter
   * @param E a double parameter
   * @param r2 a double parameter
   * @return a double value
   */
  G4double ref(G4double &x1, G4double &x2, G4double &x3, G4double p1,
		 G4double p2, G4double p3, G4double E, G4double r2);

  /**
   * ForceAbsor
   */
  void forceAbsor(G4int *nopart, G4int *iarem, G4int *izrem, G4double *esrem, G4double *erecrem,
		  G4double *alrem, G4double *berem, G4double *garem, G4int *jrem);

  /**
   * ForceAbs
   * @param iprojo projectile
   * @param at target mass number
   * @param zt target charge number 
   * @param ep projectile energy
   * @param bmax a double parameter
   * @param pt a double parameter
   * @return absorption probability
   */
  G4double forceAbs(G4double iprojo, G4double at, G4double zt, G4double ep, G4double bmax, G4double pt);

  /**
   * absoprption xsec revised version rkt-97/5                             
   * neutron data from barashenkov                                         
   * this gives absorption xsec for given zp,ap,zt,at,e (mev/nucleon)      
   * arguement changed to mev; then e=ep/ap mev/nucleon                    
   * can be used for neutrons also.                                        
   * this has coulomb as ours                                          
   * @param zp projectile charge number
   * @param zp projectile mass number
   * @param zt a double parameter
   * @param zt target charge number
   * @param at target mass number
   * @param ep projectile energy
   */
  G4double xabs2(G4double zp, G4double ap, G4double zt, G4double at, G4double ep);

  /**
   * Standard random number generator.
   * @param *rndm pointer to the variable reserved for random number
   * @param *seed pointer to the random seed
   */
  void standardRandom(G4double *rndm, G4long *seed);

  /**
   * First derivative of a gaussian potential.
   * @param *rndm pointer to the variable reserved for random number
   */
  void gaussianRandom(G4double *rndm);

  /**
   * Safe exponential function which eliminates the CPU under and overflows.
   * @param x a double parameter
   * @return a double value
   */
  G4double safeExp(G4double x);

  /**
   * Nuclear radius
   * @param A mass number (double parameter)
   */
  G4double radius(G4int A);
    
  /** Parametrisation de la section efficace de réaction calculée par incl4.1
   * iprojo=1 proton incident, iprojo=2, neutron incident).
   * entre al et u, entre 10 et 100 mev protons, 20 et 100 mev neutrons.
   * bon ordre de grandeur pour les noyaux légers (c, o ...), trés faux
   * a energie sup a 100 mev.
   * (Comment needs to be translated)
   * @param projectile an integer parameter (1 = proton, 2 = neutron)
   * @param E energy of the projectile (double parameter)
   * @param A target mass number (double parameter)
   * @return cross section (double value)
   */
  G4double crossSection(G4int projectile, G4double E, G4double A); 

  /**
   * coulombTransm
   * subroutine coulomb_transm(e,fm1,z1,fm2,z2,proba)
   * calcul du coulombien dans lahet (proba de transmission ou
   * d'entree dans le potentiel nucleaire).
   * @param E energy (a double parameter)
   * @param fm1 a double parameter
   * @param z1 a double parameter
   * @param fm2 a double parameter
   * @param z2 a double parameter
   * @return a double value
   */
  G4double coulombTransm(G4double E, G4double fm1, G4double z1, G4double fm2, G4double z2);

  /**
   * Clmb1
   * @param eta a double parameter \f$\eta = c_2*z_1*z_2*\sqrt{m/E}\f$                                          
   * @param rho a double parameter \f$\rho = c_3*(r_1+r_2)*\sqrt{mE}\f$                                        
   * @return a double value
   */
  G4double clmb1(G4double rho, G4double eta, G4double *ml);

  /**
   * First derivative of a gaussian potential.
   * @param eta a double parameter \f$\eta = c_2*z_1*z_2*\sqrt{m/E}\f$                                          
   * @param rho a double parameter \f$\rho = c_3*(r_1+r_2)*\sqrt{mE}\f$                                        
   * @return a double value
   */
  G4double clmb2(G4double rho, G4double eta, G4double *t1);

public: // Utilities
  /**
   * Returns the smaller of two numbers.
   * @param a a double value
   * @param b a double value
   * @return a double value
   */
  G4double min(G4double a, G4double b);

  /**
   * Returns the smaller of two numbers.
   * @param a an integer value
   * @param b an integer value
   * @return an integer value
   */
  G4int min(G4int a, G4int b);

  /**
   * Returns the greater of two numbers.
   * @param a a double value
   * @param b a double value
   * @return a double value
   */
  G4double max(G4double a, G4double b);

  /**
   * Returns the greater of two numbers.
   * @param a an integer value
   * @param b an integer value
   * @return an integer value
   */
  G4int max(G4int a, G4int b);

  /**
   * Rounds a double to the nearest int
   * @param a double parameter
   * @return an integer value
   */
  G4int nint(G4double number);

  /**
   * Calls a function
   * @param functionChoice an integer value representing the choice of
   * function (0 = wsax, 1 = derivWsax, 2 = dmho, 3 = derivMho, 4 = derivGaus)
   * @param r a double parameter
   * @return a double value
   */
  G4double callFunction(G4int functionChoice, G4double r);

  G4double am(G4double a, G4double b, G4double c, G4double d);
  G4double pcm(G4double e, G4double a, G4double c);
  G4double sign(G4double a, G4double b);
  G4double utilabs(G4double a);
  G4double amax1(G4double a, G4double b);
  G4double w(G4double a, G4double b, G4double c, G4double d);
  G4int idnint(G4double a);

  void print_log_start_step();
  void print_log_end_step();
  void print_log_entry(G4int iavatars, G4int iselected, G4int iparticles, G4int imin);
  void print_avatars();
  void print_one_avatar(G4int index);
  void print_one_particle(G4int index);
  void print_three_vector(G4double x, G4double y, G4double z);
  void print_map();

   private:

  /*
   * (Re)Initialize INCL internal variables
   */
  void clearState();

  /**
   * Random seeds for INCL4 internal random number generators.
   */
  G4Hazard *hazard;

  /** 
   * Data structure for INCL4.
   */
  G4Dton *dton;

  /** 
   * Data structure for INCL4. Contains the Woods-Saxon potential
   * functions for target nuclei.
   */
  G4Saxw *saxw;

  /** 
   * Data structure for INCL4.
   */
  G4Ws *ws;

  /**
   * G4Spl2
   */
  G4Spl2 *spl2;

  /**
   * G4LightGausNuc
   */
  G4LightGausNuc *light_gaus_nuc;

  /**
   * G4LightNuc
   */
  G4LightNuc *light_nuc;

  /**
   * INCL input data structure
   */
  G4InclInput *calincl;

  /**
   * G4Mat
   */
  G4Mat *mat;

  /**
   *
   */
  G4Bl1 *bl1;

    /**
   *
   */
  G4Bl2 *bl2;

    /**
   *
   */
  G4Bl3 *bl3;

    /**
   *
   */
  G4Bl4 *bl4;

    /**
   *
   */
  G4Bl5 *bl5;

    /**
   *
   */
  G4Bl6 *bl6;

  /**
   *
   */
  G4Bl8 *bl8;

    /**
   *
   */
  G4Bl9 *bl9;

    /**
   *
   */
  G4Bl10 *bl10;

  /**
   *
   */
  G4Kind *kindstruct;

  /**
   * Projectile properties.
   */
  G4Bev *bev;

  /**
   *
   */
  G4Paul *paul;
  
  /**
   * Detailed information of the cascade
   */
  G4VarAvat *varavat;

  /**
   * Cascade output.
   */
  G4VarNtp *varntp;

  /**
   * For storing the results of the evaporation.
   */
  G4VarNtp *evaporationResult;

  /** 
   * Defines the verbosity of console output. Values can be between 0
   * and 4 where 0 means silent and 4 the most verbose possible
   * output.
   */
  G4int verboseLevel;

  /**
   * Function ID for wsax.
   * @see wsax
   * @see integrate
   * @see callFunction
   */
  G4int wsaxFunction;

  /**
   * Function ID for derivWsax.
   * @see derivWsax
   * @see integrate
   * @see callFunction
   */
  G4int derivWsaxFunction;

  /**
   * Function ID for dmho.
   * @see derivWsax
   * @see integrate
   * @see callFunction
   */
  G4int dmhoFunction;

  /**
   * Function ID for derivMho.
   * @see derivMho
   * @see integrate
   * @see callFunction
   */
  G4int derivMhoFunction;

  /**
   * Function ID for derivGaus.
   * @see derivGaus
   * @see integrate
   * @see callFunction
   */
  G4int derivGausFunction;

  /**
   * Function ID for dens.
   * @see dens
   * @see integrate
   * @see callFunction
   */
  G4int densFunction;

  /**
   * Use fermi break-up?
   */
  G4bool useFermiBreakup;

  /**
   * Projectile spectator nucleus support?
   */
  G4bool useProjSpect;

  /**
   * Type of the particle at index i.
   *
   * The extension to large composite projectiles the value is
   * negative for projectile spectators. Otherwise the particle types
   * are given in exactly the same way as in G4Calincl::f[6].
   * @see G4Calincl::f
   */
  G4int kind[300]; //= (*kind_p);
  G4double ep[300]; // = (*ep_p);
  G4double alpha[300]; // = (*alpha_p); 
  G4double beta[300]; // = (*beta_p);
  G4double gam[300]; // = (*gam_p);

  G4VBe *be;
  G4InclProjSpect *ps;
  G4InclFermi *fermi;
  G4QuadvectProjo *qvp;
  G4Volant *volant;
  G4Abla *abla;
  G4InclRandomInterface *randomGenerator;

  G4VInclLogger *theLogger;
  G4int inside_step; // Flag to determine whether we are inside or outside a simulation step
};

#endif
