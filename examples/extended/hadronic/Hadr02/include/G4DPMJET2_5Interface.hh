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
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4DPMJET2_5Interface.hh
/// \brief Definition of the G4DPMJET2_5Interface class
//
// $Id: G4DPMJET2_5Interface.hh 77519 2013-11-25 10:54:57Z gcosmo $
//

#ifndef G4DPMJET2_5Interface_hh
#define G4DPMJET2_5Interface_hh 1

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5Interface.hh
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
const int intmd  =   252;
const int maxpro =     8;
const int intmx  =  2488;
const int mxafbk =    16;
const int mxffbk =     6;
const int mxnfbk =    10;
const int mxzfbk =     9;
const int mxpsst =   300;
const int mxpsfb = 41000;
const int nxafbk =    17;                        // NXAFBK = MXAFBK + 1 
const int nxnfbk =    15;                        // NXNFBK = MXNFBK + MXFFBK / 3
const int nxzfbk =    14;                        // NXZFBK = MXZFBK + MXFFBK / 3
const int nmxhkk = 89998;
const int neb    =    50;
//
//
// This next line is required as the default version of FORTRAN LOGICAL is
// four bytes long, whereas storage for G4bool is one byte.
//
typedef G4int ftnlogical;

struct ccdpm25casadi
{
  G4double casaxx;
  G4int    icasad;
};
struct ccdpm25cmhico
{
  G4double cmhis;
};
struct ccdpm25cronin
{
  G4double cronco;
  G4int    mkcron;
};
struct ccdpm25colle
{
  G4int    nevhad, nvers, ihadrz, nfile;
};
struct ccdpm25collis
{
  G4double s;                       //sometimes called ss
  G4int    ijprox, ijtar;
  G4double ptthr, ptthr2;
  G4int    iophrd, ijprlu, ijtalu;
};
struct ccdpm25coulo
{
  G4int    icoul;
};
struct ccdpm25diffra
{
  G4int    isingd, idiftp, ioudif, iflagd;
};
struct ccdpm25diqsum
{
  G4int ndvuu, ndvus, ndvss, nvduu, nvdus, nvdss,
        ndsuu, ndsus, ndsss, nsduu, nsdus, nsdss,
        ndzuu, ndzus, ndzss, nzduu, nzdus, nzdss,
        nadvuu, nadvus, nadvss, navduu, navdus, navdss,
        nadsuu, nadsus, nadsss, nasduu, nasdus, nasdss,
        nadzuu, nadzus, nadzss, nazduu, nazdus, nazdss;
};
struct ccdpm25diquax
{
  G4double amedd;
  G4int    idiqua, idiquu;
};
struct ccdpm25diqrej
{
  G4int idiqre[7], idvre[3], ivdre[3], idsre[3], isdre[3],
    idzre[3], izdre[3], idiqrz[7];
};
struct ccdpm25dprin
{
  G4int    ipri, ipev, ippa, ipco, init, iphkk, itopd, ipaupr;
};
struct ccdpm25dropjj
{
  G4double dropjt, dropva;
};
struct ccdpm25droppt
{
  ftnlogical
           intpt, fermp, ihadss, ihadsv, ihadvs, ihadvv, ihada,
           ipadis, ishmal, lpauli;
};
struct ccdpm25edens
{
  G4int    ieden;
};
struct ccdpm25evappp
{
  G4int    ievap;
};
struct ccdpm25ferfor
{
  G4int    iferfo;
};
struct ccdpm25final
{
  G4int    ifinal;
};
struct ccdpm25fluctu
{
  G4int    ifluct;
};
struct ccdpm25secint
{
  G4int    isecin;
};
struct ccdpm25frbkcm
{
  G4double amufbk, eexfbk[mxpsst], amfrbk[mxpsst],
           exfrbk[mxpsfb], sdmfbk[mxpsfb], coufbk[mxpsfb],
           exmxfb, r0frbk, r0cfbk, c1cfbk, c2cfbk;
  G4int    ifrbkn[mxpsst], ifrbkz[mxpsst],
           ifbksp[mxpsst], ifbkpr[mxpsst], ifbkst[mxpsst],
           ipsind[2][mxzfbk+1][mxnfbk+1], jpsind[mxafbk+1],
           ifbind[2][mxzfbk+1][nxnfbk+1], jfbind[nxafbk+1],
           ifbcha[mxpsfb][5], iposst, iposfb, ifbstf,
           ifbfrb, nbufbk;
  ftnlogical
           lfrmbk, lncmss;
};
struct ccdpm25gluspl
{
  G4int    nugluu, nsgluu;
};
struct ccdpm25hadthr
{
  G4double ehadth;
  G4int    inthad;
};
struct ccdpm25hdjase
{
  G4int    nhse1, nhse2, nhse3, nhase1, nhase2, nhase3;
};
struct ccdpm25hettp
{
  G4int    nhstp, nbertp, iosub, insrs;
};
struct ccdpm25ifragm
{
  G4int    ifrag;
};
//struct ccdpm25infore
//{
//  G4int    ifrej;
//};
struct ccdpm25inpflg
{
  G4int    iang, ifiss, ib0, igeom, istrag, keydk;
};
struct ccdpm25kglaub
{
  G4int    jglaub;
};
struct ccdpm25nstari
{
  G4int    nstart;
};
struct ccdpm25ncshxx
{
  G4int    ncouxh, ncouxt;
};
struct ccdpm25nncms
{
  G4double gamcm, bgcm, umo, pcm, eproj, pproj;
};
struct ccdpm25nucc
{
  G4int    it, itz, ip, ipz, ijproj, ibproj, ijtarg, ibtarg;
};
struct ccdpm25nuccc
{
  G4int    jt, jtz, jp, jpz, jjproj, jbproj, jjtarg, jbtarg;
};
struct ccdpm25nucimp
{
  G4double prmom[248][5], tamom[248][5], prmfep, prmfen, tamfep,
           tamfen, prefep, prefen, taefep, taefen, prepot[210],
           taepot[210], prebin, taebin, fermod, etacou;
};
struct ccdpm25nuclea
{
  G4double pfermp[2], pfermn[2], fermdd,
           ebindp[2], ebindn[2], epot[210][2],
           etacoo[2];
  G4int    icoull;
};
struct ccdpm25parevt
{
  G4double dpower, fsprd0, fshpfn, rn1gsc, rn2gsc;
  ftnlogical
           ldiffr[39], lpower, linctv, levprt, lheavy,
           ldeexg, lgdhpr, lpreex, lhlfix, lprfix, lparwv;
  G4int    ilvmod, jlvmod;
  ftnlogical
           llvmod, lsngch, lschdf;
};
struct ccdpm25pomtab
{
  G4int    ipomta;
};
struct ccdpm25pomtyp
{
  G4int    ipim, icon, isig, lmax, mmax, nmax;
  G4double difel, difnu;
};
struct ccdpm25popcor
{
  G4double pdb, ajsdef;
};
struct ccdpm25popcck
{
  G4double pdbck, pdbse, pdbseu;
  G4int    ijpock, irejck, ick4, ihad4, ick6, ihad6,
           irejse, ise4, ise6, irejs3, ise43, ise63, irejs0,
           ihada4, ihada6, irejsa, isea4, isea6, ireja3,
           isea43, isea63, ireja0;
};
//struct ccdpm25inxdpm
//{
//  G4int    intdpm;
//};
struct ccdpm25projk
{
  G4int    iprojk;
};
struct ccdpm25promu
{
  G4int    ipromu;
};
struct ccdpm25pshow
{
  G4int    ipshow;
};
struct ccdpm25ptlarg
{
  G4double xsmax;
};
struct ccdpm25ptsamp
{
  G4int    isampt;
};
struct ccdpm25pydat1
{
  G4int    mstu[200];
  G4double paru[200];
  G4int    mstj[200];
  G4double parj[200];
};
struct ccdpm25recom
{
  G4int    irecom;
};
struct ccdpm25seadiq
{
  ftnlogical
           lseadi;
};
struct ccdpm25seaqxx
{
  G4double seaqx, seaqxn;
};
struct ccdpm25seasu3
{
  G4double seasq;
};
struct ccdpm25sincha
{
  G4int isicha;
};
struct ccdpm25stars
{
  G4int istar2, istar3;
};
struct ccdpm25strufu
{
  G4int istrum, istrut;
};
struct ccdpm25taufo
{
  G4double taufor;
  G4int    ktauge, itauve, incmod;
};
struct ccdpm25user1
{
  char     titled[80], projty[8], targty[8];
};
struct ccdpm25user2
{
  G4double cmener, sdfrac, ptlar;
  G4int    istruf, isingx, idubld;
};
struct ccdpm25vxsvd
{
  G4double vxsp[50], vxst[50], vxsap[50], vxsat[50],
           vxvp[50], vxvt[50], vxdp[50], vxdt[50];
  G4int    nxsp, nxst, nxsap, nxsat, nxvp, nxvt, nxdp, nxdt;
};
struct ccdpm25xseadi
{
  G4double xseacu, unon, unom, unosea, cvq, cdq, csea, ssmima,
           ssmimq, vvmthr;
};
struct ccdpm25zentra
{
  G4int    icentr;
};
struct ccdpm25bufueh
{
  G4double annvv, annss, annsv, annvs, anncc,
           anndv, annvd, annds, annsd,
           annhh, annzz,
           ptvv, ptss, ptsv, ptvs, ptcc, ptdv, ptvd, ptds, ptsd,
           pthh, ptzz,
           eevv, eess, eesv, eevs, eecc, eedv, eevd, eeds, eesd,
           eehh, eezz,
           anndi, ptdi, eedi,
           annzd, anndz, ptzd, ptdz, eezd, eedz;
};

struct ccdpm25bufues
{
  G4double bnnvv, bnnss, bnnsv, bnnvs, bnncc,
           bnndv, bnnvd, bnnds, bnnsd,
           bnnhh, bnnzz,
           bptvv, bptss, bptsv, bptvs, bptcc, bptdv,
           bptvd, bptds, bptsd,
           bpthh, bptzz,
           beevv, beess, beesv, beevs, beecc, beedv,
           beevd, beeds, beesd,
           beehh, beezz,
           bnndi, bptdi, beedi,
           bnnzd, bnndz, bptzd, bptdz, beezd, beedz;
};

struct ccdpm25ncouch
{
  G4double acouvv, acouss, acousv, acouvs,
           acouzz, acouhh, acouds, acousd,
           acoudz, acouzd, acoudi,
           acoudv, acouvd, acoucc;
};

struct ccdpm25ncoucs
{
  G4double bcouvv, bcouss, bcousv, bcouvs,
           bcouzz, bcouhh, bcouds, bcousd,
           bcoudz, bcouzd, bcoudi,
           bcoudv, bcouvd, bcoucc;
};
struct ccdpm25dshm
{
  G4double rash, rbsh, bmax, bstep, sigsh, rosh, gsh,
           bsite[200][2];
  G4int    nstatb, nsiteb;
};
struct ccdpm25rptshm
{
  G4double rproj,rtarg,bimpac;
};
struct ccdpm25dtumat
{
  G4double bsiten[50][24][200], bsitem[50][24][200],
           rprojj[50], rtagg[50], bstepp[50], bmaxx[50],
           ntaxx[50], nztaxx[50], nprxx[50], nzprxx[50];
};
struct ccdpm25collap
{
  G4double s3;
  G4int    ijproj1, ijtar1;
  G4double ptthr1, ptthr3;
  G4int    iophrd1, ijprlu1,
           ijtalu1;
};
struct ccdpm25diqi
{
  G4int    ipvq[248], ippv1[248], ippv2[248], itvq[248],
           ittv1[248], ittv2[248], ipsq[intmx], ipsq2[intmx],
           ipsaq[intmx], ipsaq2[intmx], itsq[intmx],
           itsq2[intmx], itsaq[intmx], itsaq2[intmx],
           kkproj[248], kktarg[248];
};
struct ccdpm25dpar
{
  char     aname[210];
  G4double aam[210], ga[210], tau[210];
  G4int    iich[210],
           iibar[210], k1[210], k2[210];
};
struct ccdpm25extevt
{
  G4int    idres[nmxhkk], idxres[nmxhkk], nobam[nmxhkk],
           idbam[nmxhkk], idch[nmxhkk], npoint[10];
};
struct ccdpm25hkkevt
{
  G4int    nhkk, nevhkk, isthkk[nmxhkk], idhkk[nmxhkk],
           jmohkk[nmxhkk][2], jdahkk[nmxhkk][2];
  G4double phkk[nmxhkk][5],
           vhkk[nmxhkk][4], whkk[nmxhkk][4];
};
struct ccdpm25ifroto
{
  G4int    ifrovp[248], itovp[248], ifrosp[intmx],
           ifrovt[248], itovt[248], ifrost[intmx],
           jsshs[intmx], jtshs[intmx], jhkknp[248],
           jhkknt[248],
           jhkkpv[intmx], jhkkps[intmx],
           jhkktv[intmx], jhkkts[intmx],
           mhkkvv[intmx], mhkkss[intmx],
           mhkkvs[intmx], mhkksv[intmx],
           mhkkhh[intmx],
           mhkkdv[248], mhkkvd[248],
           mhkkds[intmd], mhkksd[intmd];
};
struct ccdpm25paname
{
  char     btype[30][8];  
};
struct ccdpm25shmakl
{
  G4int    jssh[intmx], jtsh[intmx], inter1[intmx],
           inter2[intmx];
};
struct ccdpm25sigma
{
  G4double sigsof, bs, zsof, sighar, fill[7];
};
struct ccdpm25xsecpt
{
  G4double ptcut, sigs, dsigh;
};
struct ccdpm25nucros
{
  G4double dsigsu, dsigmc;
  G4int    ndsig;
};
struct ccdpm25hboo
{
  G4int ihbook;
};
struct ccdpm25xsecnu
{
  G4double ecmuu, ecmoo;
  G4int    ngritt, nevtt;
};
struct ccdpm25glaber
{
  G4double ecmnn[neb], ecmnow,
           xstot[neb], xsela[neb],
           xsqep[neb], xsqet[neb],
           xsqe2[neb], xspro[neb],
           xetot[neb], xeela[neb],
           xeqep[neb], xeqet[neb],
           xeqe2[neb], xepro[neb],
           bslope,     elabb[neb];
};

extern "C"
{
 
 extern void parpt_ (int*, double*, double*, int*, int*);
 extern void csj1mi_ (double*, double*);
 extern void ddatar_ ();
 extern void dhadde_ ();
 extern void dchant_ ();
 extern void dchanh_ ();
 extern void defaul_ (double*, double*);
 extern void defaux_ (double*, double*);
 extern void lundin_ ();
 extern void rndmst_ (int*, int*, int*, int*);
 extern void berttp_ ();
 extern void incini_ ();
 extern void distr_ (int*, int*, double*, int*);
 extern void shmakf_ (int*, int*, int*, int*);
 extern void shmaki_ (int*, int*, int*, int*, double*, double*, double*);
 extern void prblm2_ (double*);
 extern void jtdtu_ (int*);
 extern void samppt_ (int*, double*);
 extern void dpmevt_ (double*, int*, int*, int*, int*, int*, int*, int*);
 extern void kkinc_ (double*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
 extern double rd2in_ (int*, int*);
 extern double rd2out_ (int*, int*);
 extern void xsglau_ (int*, int*, int*, int*);

 extern void g4dpmjet_initialise_block_data_ ();
 extern void g4dpmjet_open_nuclear_bin_ (int*, int*, int*, char*);
 extern void g4dpmjet_close_nuclear_bin_ (int*);
 extern void g4dpmjet_open_fort6_ (int*, int*, char*);
 extern void g4dpmjet_close_fort6_ ();

 extern struct ccdpm25casadi casadi_;
 extern struct ccdpm25cmhico cmhico_;
 extern struct ccdpm25cronin cronin_;
 extern struct ccdpm25colle  colle_;
 extern struct ccdpm25collis collis_;
 extern struct ccdpm25coulo  coulo_;
 extern struct ccdpm25diffra diffra_;
 extern struct ccdpm25diqsum diqsum_;
 extern struct ccdpm25diquax diquax_;
 extern struct ccdpm25diqrej diqrej_;
 extern struct ccdpm25dprin  dprin_;
 extern struct ccdpm25dropjj dropjj_;
 extern struct ccdpm25droppt droppt_;
 extern struct ccdpm25edens  edens_;
 extern struct ccdpm25evappp evappp_;
 extern struct ccdpm25ferfor ferfor_;
 extern struct ccdpm25final  final_;
 extern struct ccdpm25fluctu fluctu_;
 extern struct ccdpm25secint secint_;
 extern struct ccdpm25frbkcm frbkcm_;
 extern struct ccdpm25gluspl gluspl_;
 extern struct ccdpm25hadthr hadthr_;
 extern struct ccdpm25hdjase hdjase_;
 extern struct ccdpm25hettp  hettp_;
 extern struct ccdpm25ifragm ifragm_;
 //extern struct ccdpm25infore infore_;
 extern struct ccdpm25inpflg inpflg_;
 extern struct ccdpm25kglaub kglaub_;
 extern struct ccdpm25nstari nstari_;
 extern struct ccdpm25ncshxx ncshxx_;
 extern struct ccdpm25nncms  nncms_;
 extern struct ccdpm25nucc   nucc_;
 extern struct ccdpm25nuccc  nuccc_;
 extern struct ccdpm25nucimp nucimp_;
 extern struct ccdpm25nuclea nuclea_;
 extern struct ccdpm25parevt parevt_;
 extern struct ccdpm25pomtab pomtab_;
 extern struct ccdpm25pomtyp pomtyp_;
 extern struct ccdpm25popcor popcor_;
 extern struct ccdpm25popcck popcck_;
 //extern struct ccdpm25inxdpm inxdpm_;
 extern struct ccdpm25projk  projk_;
 extern struct ccdpm25promu  promu_;
 extern struct ccdpm25pshow  pshow_;
 extern struct ccdpm25ptlarg ptlarg_;
 extern struct ccdpm25ptsamp ptsamp_;
 extern struct ccdpm25pydat1 pydat1_;
 extern struct ccdpm25recom  recom_;
 extern struct ccdpm25seadiq seadiq_;
 extern struct ccdpm25seaqxx seaqxx_;
 extern struct ccdpm25seasu3 seasu3_;
 extern struct ccdpm25sincha sincha_;
 extern struct ccdpm25stars  stars_;
 extern struct ccdpm25strufu strufu_;
 extern struct ccdpm25taufo  taufo_;
 extern struct ccdpm25user1  user1_;
 extern struct ccdpm25user2  user2_;
 extern struct ccdpm25vxsvd  vxsvd_;
 extern struct ccdpm25xseadi xseadi_;
 extern struct ccdpm25zentra zentra_;
 extern struct ccdpm25bufueh bufueh_;
 extern struct ccdpm25bufues bufues_;
 extern struct ccdpm25ncouch ncouch_;
 extern struct ccdpm25ncoucs ncoucs_;
 extern struct ccdpm25dshm   dshm_;
 extern struct ccdpm25rptshm rptshm_;
 extern struct ccdpm25dtumat dtumat_;
 extern struct ccdpm25collap collap_;
 extern struct ccdpm25diqi   diqi_;
 extern struct ccdpm25dpar   dpar_;
 extern struct ccdpm25extevt extevt_;
 extern struct ccdpm25hkkevt hkkevt_;
 extern struct ccdpm25ifroto ifroto_;
 extern struct ccdpm25paname paname_;
 extern struct ccdpm25shmakl shmakl_;
 extern struct ccdpm25sigma  sigma_;
 extern struct ccdpm25xsecpt xsecpt_;
 extern struct ccdpm25nucros nucros_;
 extern struct ccdpm25hboo   hboo_;
 extern struct ccdpm25xsecnu xsecnu_;
 extern struct ccdpm25glaber glaber_;
}
// NOTE Should there be a semicolon after the close-curly bracket ??  Some
// compilers think yes, others no!
#endif
