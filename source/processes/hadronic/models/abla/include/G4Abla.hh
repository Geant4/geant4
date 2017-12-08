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
// Jose Luis Rodriguez, CEA (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4Abla_hh
#define G4Abla_hh 1

#ifdef ABLAXX_IN_GEANT4_MODE
#include "globals.hh"
#else
#include "G4INCLGeant4Compat.hh"
#include "G4INCLConfig.hh"
#endif

#include "G4AblaRandom.hh"
#include "G4AblaDataDefs.hh"

/**
 *  Class containing ABLA++ de-excitation code.
 */

class G4Abla {

public:
  /**
   * This constructor is used by standalone test driver and the Geant4 interface.
   *
   * @param aHazard random seeds
   * @param aVolant data structure for ABLA output
   * @param aVarNtp data structure for transfering ABLA output to Geant4 interface
   */
#ifdef ABLAXX_IN_GEANT4_MODE
  G4Abla(G4Volant *aVolant, G4VarNtp *aVarntp);
#else
  G4Abla(G4INCL::Config *config, G4Volant *aVolant, G4VarNtp *aVarntp);
#endif

  /**
   * Basic destructor.
   */
  ~G4Abla();

  /// \brief Dummy copy constructor
  G4Abla(G4Abla const &other);

  /// \brief Dummy assignment operator
  G4Abla &operator=(G4Abla const &other);

  /**
   * Set verbosity level.
   */
  void setVerboseLevel(G4int level);

  /**
   * Get the internal output data structure pointer.
   */
  G4Volant* getVolant() {
    return volant;
  }

  /**
   * Main interface to the de-excitation code.
   *
   * @param nucleusA mass number of the nucleus
   * @param nucleusZ charge number of the nucleus
   * @param excitationEnergy excitation energy of the nucleus
   * @param angularMomentum angular momentum of the nucleus (produced as output by INCL4)
   * @param momX momentum x-component
   * @param momY momentum y-component
   * @param momZ momentum z-component
   * @param eventnumber number of the event
   */
  void DeexcitationAblaxx(G4int nucleusA, G4int nucleusZ, G4double excitationEnergy, G4double angularMomentum, G4double momX, G4double momY, G4double momZ, G4int eventnumber);

  // Evaporation
public:
  /**
   * Initialize ABLA evaporation code.
   *
   */
  void initEvapora();

  /**
   * Initialize ABLA parameters.
   *
   */
  void SetParameters();
  void SetParametersG4(G4int z, G4int a);

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
  G4double spdef(G4int a, G4int z, G4int optxfis);

  /**
   * Calculation of fissility parameter
   */
  G4double fissility(int a,int z, int optxfis);

  /**
   * Main evaporation routine.
   */
  void evapora(G4double zprf, G4double aprf, G4double *ee_par, G4double jprf, 
	       G4double *zf_par, G4double *af_par, G4double *mtota_par,
	       G4double *vleva_par, G4double *vxeva_par, G4double *vyeva_par,
	       G4int *ff_par, G4int *fimf_par, G4double *fzimf, G4double *faimf, G4double *tkeimf_par,G4double *jprfout,G4int *inttype_par, G4int *inum_par,G4double EV_TEMP[200][5],G4int *iev_tab_temp_par);

  /**
   * Calculation of particle emission probabilities.
   */
  void direct(G4double zprf,G4double a, G4double ee, G4double jprf, G4double *probp_par, G4double *probd_par, G4double *probt_par, G4double *probn_par, G4double *probhe_par, G4double *proba_par, G4double *probg_par,G4double *probimf_par, G4double *probf_par, G4double *ptotl_par, G4double *sn_par, G4double *sbp_par,  G4double *sbd_par,  G4double *sbt_par, G4double *sbhe_par, G4double *sba_par, G4double *ecn_par, G4double *ecp_par,G4double *ecd_par,G4double *ect_par,G4double *eche_par, G4double *eca_par, G4double *ecg_par, G4double *bp_par, G4double *bd_par, G4double *bt_par, G4double *bhe_par, G4double *ba_par,G4double *sp,G4double *sd,G4double *st,G4double *she,G4double *sa, G4double * ef, G4double *ts1, G4int inttype, G4int inum, G4int itest, G4int *sortie, G4double *tcn,
G4double *jprfn, G4double *jprfp, G4double *jprfd, G4double *jprft, G4double *jprfhe, G4double *jprfa, G4double *tsum);

  /**
   * Calculation of fission and the particle emission probabilities after fission.
   */
void fission(G4double AF,G4double ZF,G4double EE,G4double JPRF,
        G4double *VX1_FISSION,G4double *VY1_FISSION,G4double *VZ1_FISSION,
        G4double *VX2_FISSION,G4double *VY2_FISSION,G4double *VZ2_FISSION,
        G4int *ZFP1,G4int *AFP1,G4int *ZFP2,G4int *AFP2,G4int *imode, 
        G4double *VX_EVA_SC, G4double *VY_EVA_SC, G4double *VZ_EVA_SC,
        G4double EV_TEMP[200][5],G4int *IEV_TAB_FIS);

  /**
   * Calculation of lorentz's boost
   */
void lorentz_boost(G4double VXRIN,G4double VYRIN,G4double VZRIN,G4double VXIN,G4double VYIN,G4double VZIN,G4double *VXOUT,G4double *VYOUT,G4double *VZOUT);

  /**
   * Calculation of unstable nuclei
   */
void unstable_nuclei(G4int AFP,G4int ZFP,G4int *AFPNEW,G4int *ZFPNEW,G4int &IOUNSTABLE,G4double VX,G4double VY,G4double VZ,G4double *VP1X,G4double *VP1Y,G4double *VP1Z,G4double BU_TAB_TEMP[200][5],G4int *ILOOP);

  /**
   * Calculation of unstable nuclei tke
   */
void unstable_tke(G4double AIN,G4double ZIN,G4double ANEW,G4double ZNEW,G4double VXIN,G4double VYIN,G4double VZIN,G4double *V1X,G4double *V1Y,G4double *V1Z,G4double *V2X,G4double *V2Y,G4double *V2Z);

  /**
   * Calculation of tke for breakup fragments
   */
void tke_bu(G4double Z,G4double A,G4double ZALL,G4double AAL,G4double *VX,G4double *VY,G4double *VZ);

  /**
   * Calculation of the angular momentum of breakup fragments
   * according to Goldhaber model
   */
void AMOMENT(G4double AABRA,G4double APRF,G4int IMULTIFR,G4double *PX,G4double *PY,G4double *PZ);

  /**
   * Calculation of particle emission barriers.
   */
  void barrs(G4int Z1,G4int A1,G4int Z2,G4int A2,G4double *sBARR,G4double *sOMEGA);

  /**
   * Calculation of particle emission between the saddle and scission point.
   */
  void evap_postsaddle(G4double A, G4double Z, G4double E_scission_pre, G4double *E_scission_post, G4double *A_scission, G4double *Z_scission,
G4double &vx_eva,G4double &vy_eva,G4double &vz_eva);

  /**
   * Calculation of imfs.
   */
  void imf(G4double ACN,G4double ZCN,G4double TEMP,G4double EE,G4double *ZIMF,G4double *AIMF,G4double *BIMF,G4double *SBIMF,G4double *TIMF,G4double JPRF);

  /**
   * Calculation of omega at saddle point.
   */
  void fomega_sp(G4double AF,G4double Y,G4double *MFCD,G4double *sOMEGA,G4double *sHOMEGA);

  /**
   * Calculation of omega at ground state.
   */
  void fomega_gs(G4double AF,G4double ZF,G4double *K1,G4double *sOMEGA,G4double *sHOMEGA);

  /**
   * Calculation of tunnelling effect in fission.
   */
  G4double tunnelling(G4double A,G4double ZPRF,G4double Y,G4double EE,G4double EF,G4double TEMP,G4double DENSG,G4double DENSF,G4double ENH_FACT);

  /**
   * Calculation of fission width at the saddle point according to B&W.
   */
  void fission_width(G4double ZPRF,G4double A,G4double EE,G4double BS,G4double BK,G4double EF,G4double Y,G4double *GF,G4double *TEMP,G4double JPR,G4int IEROT,G4int FF_ALLOWED,G4int OPTCOL,G4int OPTSHP,G4double DENSG);

  /**
   * Calculation of unbound nuclei.
   */
void unbound(G4double SN,G4double SP,G4double  SD,G4double ST,G4double SHE,G4double SA,G4double BP,G4double BD,G4double BT,G4double BHE,G4double BA,G4double *PROBF,G4double *PROBN,G4double *PROBP,G4double *PROBD,G4double *PROBT,G4double *PROBHE,G4double *PROBA,G4double *PROBIMF,G4double *PROBG,G4double *ECN,G4double *ECP,G4double *ECD,G4double *ECT,G4double *ECHE,G4double *ECA);

  /**
   * Calculation of the fission distribution.
   */
  void fissionDistri(G4double &a,G4double &z,G4double &e,
		   G4double &a1,G4double &z1,G4double &e1,G4double &v1,
		   G4double &a2,G4double &z2,G4double &e2,G4double &v2,
                   G4double &vx_eva_sc,G4double &vy_eva_sc,G4double &vz_eva_sc);

  /**
   * Calculation of even-odd effects in fission.
   */
  void even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out);

  /**
   * Functions for the fission model.
   */
  G4double umass(G4double z,G4double n,G4double beta);
  G4double ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d);
  G4double Uwash(double E, double Ecrit,double Freduction,double gamma);
  G4double frldm(double z,double n,double beta);
  G4double eflmac_profi(double a,double z);
  G4double gausshaz(int k, double xmoy, double sig);
  G4double haz(G4int k);

  /**
   * Level density parameters.
   */
  void densniv(G4double a, G4double z, G4double ee, G4double ef, G4double *dens, G4double bshell, G4double bs, G4double bk, 
	       G4double *temp, G4int optshp, G4int optcol, G4double defbet, G4double *ecor, G4double jprf, G4int ifis,G4double *qr);

  /**
   * Calculation of the fission probability modified by transient time effects.
   */
  void part_fiss(G4double BET,G4double GP,G4double GF,G4double Y,G4double TAUF,G4double TS1,G4double TSUM,G4int *CHOICE,G4double ZF,G4double AF,G4double FT,G4double *T_LAPSE,G4double *GF_LOC);

  G4double func_trans(G4double TIME,G4double ZF,G4double AF,G4double BET,G4double Y,G4double FT,G4double T_0);

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
   * Calculation of decay widths for light particles.
   */
  G4double width(G4double AMOTHER,G4double ZMOTHER,G4double APART,G4double ZPART,G4double TEMP,G4double B1,G4double SB1,G4double EXC);

  /**
   * Calculation of penetration factors for light charged particles.
   */
  G4double pen(G4double A, G4double ap, G4double omega, G4double T);

  /**
   * Calculation of mean value of orbital angular momentum.
   */
  void lorb(G4double AMOTHER,G4double ADAUGHTER,G4double LMOTHER,G4double EEFINAL,G4double *LORBITAL,G4double *SIGMA_LORBITAL);

  /**
   * Calculation of BS and BK for the nuclear-level density.
   */
  void bsbkbc(G4double A,G4double Z,G4double *BS,G4double *BK,G4double *BC);

  /**
   * Special functions used for the emission of particles.
   */
  G4double erf(G4double x);

  G4double gammp(G4double a, G4double x);

  void gcf(G4double *gammcf,G4double a,G4double x,G4double gln);

  void gser(G4double *gamser,G4double a,G4double x,G4double gln);

  G4double fvmaxhaz(G4double T);

  G4double fvmaxhaz_neut(G4double x);

  /**
   * Random numbers.
   */
  void standardRandom(G4double *rndm, G4long *seed);

  /**
   * LOGARITHM OF THE GAMM FUNCTION
   */ 
  G4double gammln(G4double xx);

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
   * tirage aleatoire dans une maxwellienne
   */
  G4double fmaxhaz_old(G4double T);

  /**
   * Random generator according to the
     powerfunction y = x**(lambda) in the range from xmin to xmax
   */
  G4int IPOWERLIMHAZ(G4double lambda,G4int xmin,G4int xmax);

  /**
   *
   */
  G4double pace2(G4double a, G4double z);

  /**
   *
   */
  void guet(G4double *x_par, G4double *z_par, G4double *find_par);

  /**
   * Limits of existing nuclei
   */
  void isostab_lim(G4int z, G4int *nmin, G4int *nmax);

  /**
   * Fill the data array for INCL
   */
  void FillData(G4int IMULTBU,G4int IEV_TAB);
    
public:
  // Utils
  G4int min(G4int a, G4int b);
  G4double min(G4double a, G4double b);
  G4int max(G4int a, G4int b);
  G4double max(G4double a, G4double b);
  G4double DSIGN(G4double a, G4double b);
  G4int ISIGN(G4int a, G4int b);
  G4int nint(G4double number);
  G4int secnds(G4int x);
  G4int mod(G4int a, G4int b);
  G4double dmod(G4double a, G4double b);
  G4double dint(G4double a);
  G4int idint(G4double a);
  G4int idnint(G4double value);
  G4double utilabs(G4double a);
  G4double dmin1(G4double a, G4double b, G4double c);
  G4Ec2sub* getFrldmTable() {
    return ec2sub;
  }

private:
  G4int verboseLevel;
  G4int ilast;
  G4double T_freeze_out_in;
  G4int IEV_TAB_SSC;
  G4double BU_TAB[200][11],EV_TAB[200][5],EV_TAB_SSC[200][5];
  G4int gammaemission;
  G4double T_freeze_out;
  G4Pace *pace;
  G4Ald *ald;
  G4Eenuc *eenuc;
  G4Ec2sub *ec2sub;
  G4Ecld *ecld; 
  G4Mexp *masses;
  G4Fb *fb;
  G4Fiss *fiss;
  G4Opt *opt;
  G4Volant *volant;
  G4VarNtp *varntp;  
#ifndef ABLAXX_IN_GEANT4_MODE
  G4INCL::Config *theConfig;
#endif
};

#endif
