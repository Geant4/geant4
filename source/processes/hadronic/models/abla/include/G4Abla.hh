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
// $Id: G4Abla.hh,v 1.1 2008-02-27 18:31:11 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "globals.hh"

#include "G4AblaDataDefs.hh"
#include "G4InclDataDefs.hh"

/**
 *  Class containing ABLA de-excitation code.
 */

class G4Abla {

public:
  /**
   * Basic constructor.
   */
  G4Abla();

  /**
   * This constructor is used by standalone test driver and the Geant4 interface.
   *
   * @param aHazard random seeds
   * @param aVolant data structure for ABLA output
   * @param aVarNtp data structure for transfering ABLA output to Geant4 interface
   */
  G4Abla(G4Hazard *aHazard, G4Volant *aVolant, G4VarNtp *aVarntp);

  /**
   * Constructor that is to be used only for testing purposes.
   * @param aHazard random seeds
   * @param aVolant data structure for ABLA output   
   */
  G4Abla(G4Hazard *hazard, G4Volant *volant);

  /**
   * Basic destructor.
   */
  ~G4Abla();

  /**
   * Set verbosity level.
   */
  void setVerboseLevel(G4int level);

  /**
   * Main interface to the de-excitation code.
   *
   * @param nucleusA mass number of the nucleus
   * @param nucleusZ charge number of the nucleus
   * @param nucleusMass mass of the nucleus
   * @param excitationEnergy excitation energy of the nucleus
   * @param angularMomentum angular momentum of the nucleus (produced as output by INCL4)
   * @param recoilEnergy recoil energy of the nucleus
   * @param momX momentum x-component
   * @param momY momentum y-component
   * @param momZ momentum z-component
   * @param eventnumber number of the event
   */
  void breakItUp(G4double nucleusA, G4double nucleusZ, G4double nucleusMass, G4double excitationEnergy,
		 G4double angularMomentum, G4double recoilEnergy, G4double momX, G4double momY, G4double momZ,
		 G4int eventnumber);

  // Evaporation
public:
  /**
   * Initialize ABLA evaporation code.
   *
   */
  void initEvapora();

  /**
   * Coefficient of collective enhancement including damping                         
   * Input: z,a,bet,sig,u                                                  
   * Output: qr - collective enhancement factor                            
   * See  junghans et al., nucl. phys. a 629 (1998) 635                    
   * @param z charge number
   * @param a mass number
   * @param bet beta deformation
   * @param sig perpendicular spin cut-off factor
   * @param u Energy
   * @return Coefficient of collective enhancement   
   */
  void qrot(G4double z, G4double a, G4double bet, G4double sig, G4double u, G4double *qr);

  /**
   * Model de la goutte liquide de c. f. weizsacker.
   * usually an obsolete option
   */
  void mglw(G4double a, G4double z, G4double *el);

  /**
   * Mglms
   */
  void mglms(G4double a, G4double z, G4int refopt4, G4double *el);

  /**
   *
   */
  //  G4double spdef(G4int a, G4int z, G4int optxfis);

  /**
   * Calculation of fissility parameter
   */
  //  G4double fissility(int a,int z, int optxfis);

  /**
   * Main evaporation routine.
   */
  void evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf, 
	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
	       G4double *pleva_par, G4double *pxeva_par, G4double *pyeva_par,
	       G4int *ff_par, G4int *inttype_par, G4int *inum_par);

  /**
   * Calculation of particle emission probabilities.
   */
  void direct(G4double zprf,G4double a, G4double ee, G4double jprf, 
	      G4double *probp_par, G4double *probn_par, G4double *proba_par, 
	      G4double *probf_par, G4double *ptotl_par, G4double *sn_par, G4double *sbp_par, G4double *sba_par, G4double *ecn_par, 
	      G4double *ecp_par,G4double *eca_par, G4double *bp_par, G4double *ba_par, G4int inttype, G4int inum, G4int itest);

  /**
   * Level density parameters.
   */
  void densniv(G4double a, G4double z, G4double ee, G4double esous, G4double *dens, G4double bshell, G4double bs, G4double bk, 
	       G4double *temp, G4int optshp, G4int optcol, G4double defbet);

  /**
   * This subroutine calculates the fission barriers                                                                  
   * of the liquid-drop model of Myers and Swiatecki (1967).                                                                 
   * Analytic parameterization of Dahlinger 1982 
   * replaces tables. Barrier heights from Myers and Swiatecki                                                                
   */
  G4double bfms67(G4double zms, G4double ams);

  /**
   * This subroutine calculates the ordinary legendre polynomials of   
   * order 0 to n-1 of argument x and stores them in the vector pl.    
   * They are calculated by recursion relation from the first two      
   * polynomials.                                                      
   * Written by A.J.Sierk  LANL  t-9  February, 1984                   
   */
  void lpoly(G4double x, G4int n, G4double pl[]);

  /**
   * This function will calculate the liquid-drop nuclear mass for spheri
   * configuration according to the preprint NUCLEAR GROUND-STATE        
   * MASSES and DEFORMATIONS by P. Mo"ller et al. from August 16, 1993 p.
   * All constants are taken from this publication for consistency.      
   */
  G4double eflmac(G4int ia, G4int iz, G4int flag, G4int optshp);

  /**
   * Procedure for calculating the pairing correction to the binding   
   * energy of a specific nucleus.
   */
  void appariem(G4double a, G4double z, G4double *del);

  /**
   * PROCEDURE FOR CALCULATING THE PARITY OF THE NUMBER N.             
   * RETURNS -1 IF N IS ODD AND +1 IF N IS EVEN                        
   */
  void parite(G4double n, G4double *par);

  /**
   * RISE TIME IN WHICH THE FISSION WIDTH HAS REACHED      
   * 90 PERCENT OF ITS FINAL VALUE
   */
  G4double tau(G4double bet, G4double homega, G4double ef, G4double t);

  /**
   * KRAMERS FAKTOR  - REDUCTION OF THE FISSION PROBABILITY       
   * INDEPENDENT OF EXCITATION ENERGY 
   */
  G4double cram(G4double bet, G4double homega);

  /**
   * CALCULATION OF THE SURFACE BS OR CURVATURE BK OF A NUCLEUS        
   * RELATIVE TO THE SPHERICAL CONFIGURATION                           
   * BASED ON  MYERS, DROPLET MODEL FOR ARBITRARY SHAPES               
   */
  G4double bipol(int iflag, G4double y);

  /**
   * THIS SUBROUTINE RETURNS THE BARRIER HEIGHT BFIS, THE              
   * GROUND-STATE ENERGY SEGS, IN MEV, AND THE ANGULAR MOMENTUM        
   * AT WHICH THE FISSION BARRIER DISAPPEARS, LMAX, IN UNITS OF        
   * H-BAR, WHEN CALLED WITH INTEGER AGUMENTS IZ, THE ATOMIC           
   * NUMBER, IA, THE ATOMIC MASS NUMBER, AND IL, THE ANGULAR           
   * MOMENTUM IN UNITS OF H-BAR. (PLANCK'S CONSTANT DIVIDED BY         
   * 2*PI).                                                            
   */
  void barfit(G4int iz, G4int ia, G4int il, G4double *sbfis, G4double *segs, G4double *selmax);

  /**
   * TIRAGE ALEATOIRE DANS UNE EXPONENTIELLLE : Y=EXP(-X/T)
   */ 
  G4double expohaz(G4int k, G4double T);

  /**
   * DISTRIBUTION DE MAXWELL
   */
  G4double fd(G4double E);

  /**
   *FONCTION INTEGRALE DE FD(E)
   */
  G4double f(G4double E);

  /**
   * tirage aleatoire dans une maxwellienne
   */
  G4double fmaxhaz(G4double T);

  /**
   *
   */
  G4double pace2(G4double a, G4double z);

  /**
   *
   */
  void guet(G4double *x_par, G4double *z_par, G4double *find_par);

  // Fission
public:
  /**
   *
   */
  G4double spdef(G4int a, G4int z, G4int optxfis);

  /**
   *
   */
  G4double fissility(G4int a, G4int z, G4int optxfis);

//   void evapora(G4double zprf, G4double aprf, G4double ee, G4double jprf,
// 	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
// 	       G4double *pleva_par, G4double *pxeva_par);
//  G4double bfms67(G4double zms, G4double ams);
  //  void lpoly(G4double x, G4int n, G4double pl[]);
  //  G4double expohaz(G4int k, G4double T);
  //  G4double fd(G4double E);
  //  G4double f(G4double E);
  //  G4double fmaxhaz(G4double k, G4double T);
  void even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out);
  G4double umass(G4double z,G4double n,G4double beta);
  G4double ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d);
  void fissionDistri(G4double &a,G4double &z,G4double &e,
		     G4double &a1,G4double &z1,G4double &e1,G4double &v1,
		     G4double &a2,G4double &z2,G4double &e2,G4double &v2);
  void standardRandom(G4double *rndm, G4long *seed);
  G4double haz(G4int k);
  G4double gausshaz(int k, double xmoy, double sig);

    
public:
  // Coordinate system transformations:
  void lorab(G4double gam, G4double eta, G4double ein, G4double pin[],
	     G4double *eout, G4double pout[]);

  void translab(G4double gamrem, G4double etrem, G4double csrem[4], G4int nopart, G4int ndec);
  void translabpf(G4double masse1, G4double t1, G4double p1, G4double ctet1,
		  G4double phi1, G4double gamrem, G4double etrem, G4double R[][4],
		  G4double *plab1, G4double *gam1, G4double *eta1, G4double csdir[]);

  void rotab(G4double R[4][4], G4double pin[4], G4double pout[4]);

  // Utils
  G4int min(G4int a, G4int b);
  G4double min(G4double a, G4double b);
  G4int max(G4int a, G4int b);
  G4double max(G4double a, G4double b);

  G4int nint(G4double number);
  G4int secnds(G4int x);
  G4int mod(G4int a, G4int b);
  G4double dmod(G4double a, G4double b);
  G4double dint(G4double a);
  G4int idint(G4double a);
  G4int idnint(G4double value);
  G4double utilabs(G4double a);
  G4double dmin1(G4double a, G4double b, G4double c);

private:
  G4int verboseLevel;
  G4int ilast;

  G4Pace *pace;
  G4Hazard *hazard;
  G4Ald *ald;
  G4Ablamain *ablamain;
  G4Emdpar *emdpar;
  G4Eenuc *eenuc;
  G4Ec2sub *ec2sub;
  G4Ecld *ecld; 
  G4Fb *fb;
  G4Fiss *fiss;
  G4Opt *opt;
  G4Volant *volant;
  G4VarNtp *varntp;  
};
