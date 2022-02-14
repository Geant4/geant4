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
// *                                                                  *
// * Parts of this code which have been  developed by Abdel-Waged     *
// * et al under contract (31-465) to the King Abdul-Aziz City for    *
// * Science and Technology (KACST), the National Centre of           *
// * Mathematics and Physics (NCMP), Saudi Arabia.                    *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4UrQMD1_3Interface.hh
/// \brief Definition of the G4UrQMD1_3Interface class
//
//

#ifndef G4UrQMD1_3Interface_hh
#define G4UrQMD1_3Interface_hh

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:           G4UrQMD1_3Model.hh
//
// Version:          0.B
// Date:             20/12/12
// Author:           Kh. Abdel-Waged and Nuha Felemban
// Revised by:       V.V. Uzhinskii
//                   SPONSERED BY
// Customer:         KAUST/NCMP
// Contract:         31-465
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include "G4SystemOfUnits.hh"

//  coms
//
const G4int  nmax  = 500;
const G4int  nspl  = 500;
const G4int  smax  = 500;
//  comres
const G4int  minnuc=1;
const G4int  minmes=100;
const G4int  maxmes=132;
const G4int  numnuc=16;
const G4int  numdel=10;
const G4int  maxnuc=minnuc+numnuc-1;
const G4int  mindel=minnuc+maxnuc;
const G4int  maxdel=mindel+numdel-1;
const G4int  minres=minnuc+1;
const G4int  maxres=maxdel;
const G4int  numlam=13;
const G4int  numsig=9;
const G4int  numcas=6;
const G4int  numome=1;
const G4int  minlam=mindel+numdel;
const G4int  maxlam=minlam+numlam-1;
const G4int  minsig=minlam+numlam;
const G4int  maxsig=minsig+numsig-1;
const G4int  mincas=minsig+numsig;
const G4int  maxcas=mincas+numcas-1;
const G4int  minome=mincas+numcas;
const G4int  maxome=minome+numome-1;
const G4int  minbar=minnuc;
const G4int  maxbar=maxome;
const G4int  offmeson=minmes;
const G4int  maxmeson=maxmes;
const G4int  maxbra=11;
const G4int  maxbrm=25;
const G4int  maxbrs1=10;
const G4int  maxbrs2=3;
const G4int  nsigs = 10;
const G4int  itblsz= 100;
const G4int  maxreac = 13;
const G4int  maxpsig = 12;
//
//comwid
//
const G4int    widnsp=120;
const G4double mintab=0.10;
const G4double maxtab1=5.0;
const G4double maxtab2=50.0;
const G4int    tabver=9;
//
// options
//
const G4int numcto=400;
const G4int numctp=400;
const G4int maxstables=20;
//
// colltab (collision tables)
//
const G4int ncollmax = 100;
//  
// inputs
//
const G4int aamax=300;
//
// newpart (new created particles)
//
const G4int  mprt=200;
const G4int  oprt=2;
//
// boxinc
//
const G4int bptmax=20;
//

// This next line is required as the default version of FORTRAN LOGICAL is
// four bytes long, whereas storage for G4bool is one byte.
//
// comnorm
const G4int n = 400;
//
// comstr
const G4int njspin=8;
//
//iso
const G4int jmax=7;

// This next line is required as the default version of FORTRAN LOGICAL is
// four bytes long, whereas storage for G4bool is one byte.
//

typedef G4int ftnlogical;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Standard common block for UrQMD
// Common options for coms.f
//  20 commons
//
//
struct ccurqmd13urqmdparams
{
G4int u_at,u_zt,u_ap,u_zp;
G4double u_elab,u_imp;
G4int u_sptar,u_spproj;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13sys
{  
G4int npart, nbar, nmes, ctag,nsteps,uid_cnt,
  ranseed,event,ap,at,zp,zt,eos,dectag,
  nhardres, nsoftres, ndecres, nelcoll, nblcoll;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13rsys
{
G4double time,acttime,bdist,bimp,bmin,ebeam,ecm;
};

struct ccurqmd13comseed
{
  ftnlogical firstseed;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13logic
{
  ftnlogical lsct[nmax], logSky, logYuk, logCb, logPau;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13mdprop
{
 G4double r0_t[nmax], rx_t[nmax], ry_t[nmax], rz_t[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13cuts
{
  G4double cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13spdata
{
  G4double spx[nspl], spPauy[nspl], outPau[nspl], 
    spCby[nspl],  outCb[nspl],
    spYuky[nspl], outYuk[nspl],
    spSkyy[nspl], outSky[nspl],
    spdwwy[nspl], outdww[nspl];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13isys
{
  
G4int spin[nmax],ncoll[nmax],charge[nmax],ityp[nmax],
  lstcoll[nmax],
  iso3[nmax],origin[nmax],strid[nmax],uid[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13coor
{
G4double r0[nmax], rx[nmax], ry[nmax], rz[nmax],
  p0[nmax], px[nmax], py[nmax], pz[nmax],
  fmass[nmax], rww[nmax],dectime[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13frag
{
G4double tform[nmax], xtotfac[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13aios
{
G4double airx[nmax], airy[nmax], airz[nmax],
  aipx[nmax], aipy[nmax], aipz[nmax],
  aorx [4][nmax], aory[4][nmax], aorz[4][nmax],
  aopx[4][nmax], aopy[4][nmax], aopz[4][nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13pots
{
G4double Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
  gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
  dt,da, db,dtimestep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13scoor
{
G4double r0s[smax], rxs[smax], rys[smax], rzs[smax],
  p0s[smax], pxs[smax], pys[smax], pzs[smax],
  sfmass[smax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13sisys
{
  G4int sspin[smax], scharge[smax], sityp[smax], siso3[smax],
    suid[smax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13ssys
{
  G4int  nspec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13rtdelay
{
G4double p0td[nmax][2],pxtd[nmax][2],pytd[nmax][2],pztd[nmax][2],
  fmasstd[nmax][2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13itdelay
{
G4int ityptd[nmax][2],iso3td[nmax][2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13svinfo
{
G4int itypt[2],uidt[2],origint[2],iso3t[2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13ffermi
{
G4double ffermpx[nmax], ffermpy[nmax], ffermpz[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13peq
{
G4double peq1, peq2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Definition for Collision Term
// Commons  comres
// 4 commons
//

struct ccurqmd13versioning
{
char versiontag[45];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13resonances
{

  G4double massres[maxbar-minbar+1],widres[maxbar-minbar+1];
  G4double massmes[maxmes-minmes+1];
  G4double widmes[maxmes-minmes+1];
  G4double mmesmn[maxmes-minmes+1];
  G4double branres[maxdel-minnuc][maxbra+1];
  G4double branmes[maxmes-minmes][maxbrm+1];
      
  G4double branbs1[maxsig-minlam][maxbrs1+1];
  G4double branbs2[maxcas-mincas][maxbrs2+1];

  G4int  bs1type[maxbrs1+1][4],bs2type[maxbrs2+1][4];
  G4int lbs1[maxsig-minlam][maxbrs1+1];
  G4int lbs2[maxcas-mincas][maxbrs2+1];
  G4int lbm[maxmes-minmes][maxbrm+1];

  G4int  jres[maxbar-minbar+1];
  G4int  jmes[maxmes-minmes+1];
  G4int lbr[maxdel-minnuc][maxbra+1];
  G4int  brtype[maxbra+1][4];
  G4int  pares[maxbar-minbar+1],pames[maxmes-minmes+1];
  G4int  bmtype[maxbrm+1][4];
  G4int  isores[maxbar-minbar+1], isomes[maxmes-minmes+1];
  G4int  strres[maxbar-minbar+1],strmes[maxmes-minmes+1];
  G4int mlt2it[maxmes-minmes];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13sigtabi
{
G4int sigmaln[maxreac][2][maxpsig];
G4int sigmainf[20][nsigs];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct  ccurqmd13sigtabr
{
G4double  sigmas[itblsz][nsigs],sigmascal[5][nsigs];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//comwid
struct ccurqmd13decaywidth
{
G4double tabx [widnsp];
G4double fbtaby [2][maxbar-minbar+1][widnsp];
G4double  pbtaby[maxbra+1][maxbar-minbar+1][2][widnsp];
G4double  fmtaby [2][maxmes-minmes+1][widnsp];
G4double  pmtaby [maxbrm+1][maxmes-minmes+1][2][widnsp];
G4int     wtabflg;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13brwignorm
{
G4double bwbarnorm[maxbar-minbar+1];
G4double bwmesnorm[maxmes-minmes+1];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13xsections
{
G4double tabxnd [widnsp];
G4double frrtaby[maxdel-1][2][2][widnsp];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13tabnames
{
char tabname[77];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// options
//
struct ccurqmd13options
{
G4int    CTOption[numcto];
G4double CTParam[numctp];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13optstrings
{
char ctodc[numcto][2];
char ctpdc[numctp][2];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13loptions
{
ftnlogical fixedseed,bf13,bf14,bf15,bf16,bf17,bf18,bf19,
  bf20;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13stables
{
G4int nstable;
G4int stabvec[maxstables];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//colltab
//
struct ccurqmd13colltab
{
G4double cttime[ncollmax+1],ctsqrts[ncollmax],
  ctsigtot[ncollmax],tmin;
G4int    cti1[ncollmax],cti2[ncollmax];
G4int    nct,actcol;
ftnlogical ctvalid[ncollmax];
G4int    ctsav[ncollmax];
G4int    nsav,apt;
G4double ctcolfluc[ncollmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// inputs
//
struct ccurqmd13inputs
{
G4int  nevents,spityp[2],prspflg;
G4int  trspflg,spiso3[2],outsteps,bflag,srtflag,efuncflag;
G4int  nsrt,firstev,npb;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13input2
{
G4double srtmin,srtmax,pbeam,betann,betatar,betapro;
G4double pbmin,pbmax;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13protarints
{
G4int pt_iso3[2][aamax],pt_ityp[2][aamax],pt_spin[2][aamax];
G4int pt_charge[2][aamax],pt_aa[2],pt_uid[2][aamax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13protarreals
{
G4double pt_r0[2][aamax],pt_rx[2][aamax],pt_ry[2][aamax],
  pt_rz[2][aamax],pt_fmass[2][aamax],pt_dectime[2][aamax];
G4double pt_p0[2][aamax],pt_px[2][aamax],pt_py[2][aamax],
  pt_pz[2][aamax];
G4double pt_rho[2][aamax];
G4double pt_pmax[2][aamax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// newpart
struct ccurqmd13inewpart
{
G4int itypnew[mprt],i3new[mprt],itot[mprt],inew[mprt],nexit;
G4int iline,strcount,pslot[oprt],nstring1, nstring2,
  sidnew[mprt],itypold[oprt],iso3old[oprt];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13rnewpart
{
G4double pnew[mprt][5],xnew[mprt][4],betax,betay,betaz, 
  pold[oprt][5],p0nn,pxnn,pynn,pznn,pnn, mstring[2],
  pnnout,xtotfacold[oprt];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13fnewpart
{
G4double leadfac[mprt];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// boxinc
//
struct ccurqmd13boxic
{
G4int cbox;
G4int boxflag;
G4int mbox;
G4int bptityp[bptmax],bptiso3[bptmax],bptpart[bptmax];
G4int edensflag,para,solid, mbflag,mtest;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13boxrc
{
G4double bptpmax[bptmax];
G4double edens;
G4double lbox;
G4double lboxhalbe;
G4double lboxd;
G4double mbp0, mbpx, mbpy, mbpz;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// comnorm
struct ccurqmd13normsplin
{
G4double x_norm[n][4],y_norm[n][4];
G4double y2a[n][4],y2b[n][4], dx;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// comstr
struct ccurqmd13FRGSPA
{
G4double pjspns, pmix1s[njspin][3], pmix2s[njspin][3], pbars, parqls, parrs;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13FRGCPA
{
G4double pjspnc, pmix1c[njspin][3], pmix2c[njspin][3], pbarc;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13coparm
{
G4double parm[njspin];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct ccurqmd13const
{
G4double pi;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//// freezeout
//
struct ccurqmd13frcoor
{
G4double frr0[nmax], frrx[nmax], frry[nmax], frrz[nmax],
  frp0[nmax], frpx[nmax], frpy[nmax], frpz[nmax];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  input
struct ccurqmd13values
{
G4double valint[1];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// cascinit
struct ccurqmd13ini
{
ftnlogical bcorr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// iso
struct ccurqmd13factorials
{
G4double logfak[101];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct ccurqmd13cgks
{
G4double cgktab[jmax+1][2*jmax+1][2*jmax+1][jmax+1][jmax+1];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// UrQMD
//
struct ccurqmd13energies
{
G4double ekinbar, ekinmes, esky2, esky3, eyuk, ecb, epau;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// urqmd
extern "C"
{
extern int time_ ();
extern void loginit_();
extern void sseed_ (int*);
extern void uinit_ (int*);
extern void urqmd_ ();
extern int pdgid_ (int*, int*); //ityp

extern void g4urqmdblockdata_ ();

// urqmdparams
extern struct ccurqmd13urqmdparams  urqmdparams_;
//coms
extern struct ccurqmd13sys     sys_;
extern struct ccurqmd13rsys    rsys_;
extern struct ccurqmd13comseed comseed_;
extern struct ccurqmd13logic   logic_;
extern struct ccurqmd13mdprop  mdprop_;
extern struct ccurqmd13cuts    cuts_;
extern struct ccurqmd13spdata  spdata_;
extern struct ccurqmd13isys    isys_;
extern struct ccurqmd13coor    coor_;
extern struct ccurqmd13frag    frag_;
extern struct ccurqmd13aios    aios_;
extern struct ccurqmd13pots    pots_;
extern struct ccurqmd13scoor   scoor_;
extern struct ccurqmd13sisys   sisys_;
extern struct ccurqmd13ssys    ssys_;
extern struct ccurqmd13rtdelay rtdelay_;
extern struct ccurqmd13itdelay itdelay_;
extern struct ccurqmd13svinfo  svinfo_;
extern struct ccurqmd13ffermi  ffermi_;
extern struct ccurqmd13peq     peq_;
//comres
extern struct ccurqmd13versioning  versioning_;
extern struct ccurqmd13resonances  resonances_;
extern struct ccurqmd13sigtabi  sigtabi_;
extern struct ccurqmd13sigtabr sigtabr_;

//comwid
extern struct ccurqmd13decaywidth decaywidth_;
extern struct ccurqmd13brwignorm  brwignorm_;
extern struct ccurqmd13xsections  xsections_;
extern struct ccurqmd13tabnames   tabnames_;
//options
extern struct ccurqmd13options    options_;
extern struct ccurqmd13optstrings optstrings_;
extern struct ccurqmd13loptions    loptions_;
extern struct ccurqmd13stables     stables_;
//colltab
extern struct ccurqmd13colltab     colltab_;
//inputs
extern struct ccurqmd13inputs      inputs_;
extern struct ccurqmd13input2      input2_;
extern struct ccurqmd13protarints  protarints_;
extern struct ccurqmd13protarreals protarreals_;
//newpart
extern struct ccurqmd13inewpart    inewpart_;
extern struct ccurqmd13rnewpart    rnewpart_;
extern struct ccurqmd13fnewpart    fnewpart_;
//bocinc
extern struct ccurqmd13boxic       boxic_;
extern struct ccurqmd13boxrc       boxrc_;
// comnorm
struct ccurqmd13normsplin  normsplin_;
//comstr
struct ccurqmd13FRGSPA  FRGSPA_;
struct ccurqmd13FRGCPA  FRGCPA_;
struct ccurqmd13coparm  coparm_;
struct ccurqmd13const   const_;
// freezeout
struct ccurqmd13frcoor  frcoor_;
//urqmd
extern struct ccurqmd13energies  energies_;
//input
extern struct ccurqmd13values values_;
// cascinit
extern struct ccurqmd13ini ini_;
//iso
extern struct ccurqmd13factorials factorials_;
extern struct ccurqmd13cgks  cgks_;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
#endif
