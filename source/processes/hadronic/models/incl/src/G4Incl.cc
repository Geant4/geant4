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
// $Id: G4Incl.cc,v 1.1 2007-05-23 10:25:37 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4Incl.hh"

#include <iostream>

G4Incl::G4Incl()
  :verboseLevel(0)
{
  // ID codes for each function to be used for integration routine.
  wsaxFunction = 0;
  derivWsaxFunction = 1;
  dmhoFunction = 2;
  derivMhoFunction = 3;
  derivGausFunction = 4;
  densFunction = 5;
}

G4Incl::G4Incl(G4Hazard *hazard, G4Dton *dton, G4Saxw *saxw, G4Ws *ws)
  :verboseLevel(0),
   hazard(hazard), dton(dton), saxw(saxw), ws(ws)
{

}

G4Incl::G4Incl(G4Hazard *hazard, G4Calincl *calincl, G4Ws *ws, G4Mat *mat, G4VarNtp *varntp)
  :verboseLevel(0),
   hazard(hazard), calincl(calincl), ws(ws), mat(mat), varntp(varntp)
{
  light_gaus_nuc = (G4LightGausNuc*) malloc(sizeof(G4LightGausNuc));
  light_nuc = (G4LightNuc*) malloc(sizeof(G4LightNuc));
  spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));
  saxw = (G4Saxw*) malloc(sizeof(G4Saxw));
  dton = (G4Dton*) malloc(sizeof(G4Dton));
  bl1 = (G4Bl1*) malloc(sizeof(G4Bl1));
  bl2 = (G4Bl2*) malloc(sizeof(G4Bl2));
  bl3 = (G4Bl3*) malloc(sizeof(G4Bl3));
  bl4 = (G4Bl4*) malloc(sizeof(G4Bl4));
  bl5 = (G4Bl5*) malloc(sizeof(G4Bl5));
  bl6 = (G4Bl6*) malloc(sizeof(G4Bl6));
  bl8 = (G4Bl8*) malloc(sizeof(G4Bl8));
  bl9 = (G4Bl9*) malloc(sizeof(G4Bl9));
  bl10 = (G4Bl10*) malloc(sizeof(G4Bl10));
  kindstruct = (G4Kind*) malloc(sizeof(G4Kind));
  paul = (G4Paul*) malloc(sizeof(G4Paul));
  varavat = (G4VarAvat*) malloc(sizeof(G4VarAvat));

  volant = (G4Volant*) malloc(sizeof(G4Volant));
  abla = new G4Abla(hazard, volant);
  abla->initEvapora();
}

G4Incl::~G4Incl()
{

}

void G4Incl::setVerboseLevel(G4int level)
{
  verboseLevel = level;
}

G4int G4Incl::getVerboseLevel()
{
  return verboseLevel;
}

void G4Incl::setDtonData(G4Dton *newDton)
{
  dton = newDton;
}

void G4Incl::setWsData(G4Ws *newWs)
{
  ws = newWs;
}

void G4Incl::setHazardData(G4Hazard *newHazard)
{
  hazard = newHazard;
}

void G4Incl::setSaxwData(G4Saxw *newSaxw)
{
  saxw = newSaxw;
}

void G4Incl::setSpl2Data(G4Spl2 *newSpl2)
{
  spl2 = newSpl2;
}

void G4Incl::setCalinclData(G4Calincl *newCalincl)
{
  calincl = newCalincl;
}

void G4Incl::setMatData(G4Mat *newMat)
{
  mat = newMat;
}

void G4Incl::setLightNucData(G4LightNuc *newLightNuc)
{
  light_nuc = newLightNuc;
}

void G4Incl::setLightGausNucData(G4LightGausNuc *newLightGausNuc)
{
  light_gaus_nuc = newLightGausNuc;
}

void G4Incl::setBl1Data(G4Bl1 *newBl1)
{
  bl1 = newBl1;
}

void G4Incl::setBl2Data(G4Bl2 *newBl2)
{
  bl2 = newBl2;
}

void G4Incl::setBl3Data(G4Bl3 *newBl3)
{
  bl3 = newBl3;
}

void G4Incl::setBl4Data(G4Bl4 *newBl4)
{
  bl4 = newBl4;
}

void G4Incl::setBl5Data(G4Bl5 *newBl5)
{
  bl5 = newBl5;
}

void G4Incl::setBl6Data(G4Bl6 *newBl6)
{
  bl6 = newBl6;
}

void G4Incl::setBl8Data(G4Bl8 *newBl8)
{
  bl8 = newBl8;
}

void G4Incl::setBl9Data(G4Bl9 *newBl9)
{
  bl9 = newBl9;
}

void G4Incl::setBl10Data(G4Bl10 *newBl10)
{
  bl10 = newBl10;
}

void G4Incl::setKindData(G4Kind *newKind)
{
  kindstruct = newKind;
}

// INCL main routine for event processing

void G4Incl::processEventWithEvaporation()
{
  const G4double uma = 931.4942;
  const G4double melec = 0.511;
  const G4double fmp = 938.2723;
  const G4double fmn = 939.56563;
  
  int ntrans = 0;
  int nabs = 0;
  int nretir = 0;
  int nfis = 0;
  double muln_t = 0.0;
  double muln_2 = 0.0;
  double ener_2 = 0.0;
  double muln_20 = 0.0;
  double ener_20 = 0.0;
  double muln_max = 0.0;
  double ener_max = 0.0;
  int nbevhbk = 0;

  int impulse1 = 0; // conservation d'impulsion fausse sans fission
  int impulse2 = 0;
  int impulse3 = 0;
  int imp_f1 = 0; //      ! conservation d'impulsion fausse avec fission
  int imp_f2 = 0;
  int imp_f3 = 0;
        
  double S1x = 0.0;
  double S2x = 0.0;
  double S1y = 0.0;
  double S2y = 0.0;
  double S1z = 0.0;
  double S2z = 0.0;
  
  double ap, zp, mprojo, pbeam;

  // pi+
  if(calincl->f[6] == 3.0) { 
    mprojo = 139.56995;
    ap = 0.0;
    zp = 1.0;
  }

  // pi0
  if(calincl->f[6] == 4.0) {
    mprojo = 134.9764;
    ap = 0.0;
    zp = 0.0;
  }

  // pi-
  if(calincl->f[6] == 5.0) {
    mprojo = 139.56995;
    ap = 0.0;
    zp = -1.0;
  }

  // coulomb en entree seulement pour les particules ci-dessous

  // proton
  if(calincl->f[6] == 1.0) {
    mprojo = 938.27231;
    ap = 1.0;
    zp = 1.0;
  }

  // neutron  
  if(calincl->f[6] == 2.0) {
    mprojo = 939.56563;
    ap = 1.0;
    zp = 0.0;
  }

  // deuteron
  if(calincl->f[6] == 6.0) {
    mprojo = 1875.61276;
    ap = 2.0;
    zp = 1.0;
  }

  // triton
  if(calincl->f[6] == 7.0) {
    mprojo = 2808.95;
    ap = 3.0;
    zp = 1.0;
  }

  // He3
  if(calincl->f[6] == 8.0) {
    mprojo = 2808.42;
    ap = 3.0;
    zp = 2.0;
  }

  // Alpha
  if(calincl->f[6] == 9.0) {
    mprojo = 3727.42;
    ap = 4.0;
    zp = 2.0;
  }

  pbeam = sqrt(calincl->f[2]*(calincl->f[2] + 2.0*mprojo));         

  double at = calincl->f[0];
  double zt = calincl->f[1];
  double eap = calincl->f[2];
        
  calincl->f[3] = 0.0;    //     !seuil sortie proton
  calincl->f[7] = 0.0;  //       !seuil sortie neutron

  double proba_trans;
  
  if((calincl->f[6] == 1.0) || (calincl->f[6] >= 6.0)) {
    proba_trans = coulombTransm(calincl->f[2],ap,zp,calincl->f[0],calincl->f[1]);
  }

  double ntrans_coul = calincl->icoup*(1.0 - proba_trans);
  calincl->icoup = calincl->icoup - ntrans_coul;
  std::cout <<"icoup, proba_trans " << calincl->icoup << " " << proba_trans << std::endl;

  // old main loop
//   if(i == ievtest) then
//             write(6,*) 'le tir selectionne commence: i=',i
//         endif
        
// 990     continue        !retirage si nopart= -100 (void event)
            
// c stockage des generateurs pour reproduire les bugs...
//            ialview = ial
//                 do ihaz=1,19
//                    iyv(ihaz)=iy(ihaz)
//                 enddo
                
// c          write(6,*)'event',i

  double mtota = 0.0; // counters of a masses emitted as alpha                             
  double malpha1 = 0.0;                             
  double malpha2 = 0.0;
                                  
  int ap1 = idnint(ap);                                                  
  int zp1 = idnint(zp);                                                  
//C*********************** Appel de la CASCADE ****************************
  int ibert = 1;
//   if(i == 1) {
//     ibert = 0;
//     ibert = i;
//   }           

  int nopart;
  int izrem;
  int iarem;
  double esrem;
  double erecrem;
  double berem;
  double garem;
  double bimpac;
  int jrem;
  double alrem;

  std::cout <<"Calling pnu..." << std::endl;
  pnu(&ibert, &nopart,&izrem,
            &iarem,&esrem,&erecrem,&alrem,&berem,&garem,&bimpac,&jrem);
  std::cout <<"returned from pnu..." << std::endl;

  //  C Absorption forc%Gï¿½%@ pour p (10-100 MeV) et n (20-100MeV)
  forceAbsor(nopart, iarem, izrem, esrem, erecrem, alrem, berem, garem, jrem);
//         IF(nopart.eq.-1)THEN
//                 ntrans=ntrans+1
//                 go to 99        ! tirage non exploite, au suivant!
//         ENDIF
        
//         IF(nopart.eq.0)THEN
//                 nabs=nabs+1


// C               go to 99        ! tirage non exploite, au suivant!
// C Et si! Il faut evaporer une absorption... (Important a basse energie 
// C                                           en dessous de 200 MeV)
//         ENDIF
        
//         IF(nopart.eq.-100)THEN
//                 nretir=nretir+1 
//                 go to 990       ! Retirage
//         ENDIF



// C ----------------- OK, valuable cascade, we can continue --------------


// Generate output data:
//  C ----------------- OK, valuable cascade, we can continue --------------

  G4double zprf = double(izrem);    // nuclear charge of the prefragment
  G4double aprf = double(iarem);    // mass number of the prefragment
  G4double ee = double(esrem);        // excitation energy of the prefragment
  G4double jprf = 0.0;                // angular momentum of the prefragment

  // mean angular momentum of prefragment                                  

  jprf = 0.165 * pow(at,(2.0/3.0)) * aprf*(at - aprf)/(at - 1.0);                               
  if (jprf < 0) {
    jprf = 0.0;
  }
  // check m.de jong, ignatyuk, schmidt nuc.phys a 613, pg442, 7th line
  jprf = sqrt(2*jprf);

  jprf = jrem;
  varntp->jremn = jrem;      // jrem copie dans le ntuple

  G4double numpi = 0;  // compteurs de pions, neutrons protons
  G4double multn = 0; 
  G4double multp = 0;

  // ecriture dans le ntuple des particules de cascade (sauf remnant)      
  varntp->ntrack = nopart;          // nombre de particules pour ce tir
  varntp->massini = iarem;
  varntp->mzini = izrem;
  varntp->exini = esrem;
  varntp->bimpact = bimpac;
  std::cout <<"Bimpact: " << varntp->bimpact << std::endl;
  
  //  three ways to compute the mass of the remnant: 
  //                -from the output of the cascade and the canonic mass
  //                -from energy balance (input - all emitted energies)
  //                -following the approximations of the cugnon code (esrem...)
  G4double mcorem = mprojo + calincl->f[2] + abla->pace2(double(calincl->f[0]),double(calincl->f[1]))
    + calincl->f[0]*uma - calincl->f[1]*melec;

  G4double pxbil = 0.0;
  G4double pybil = 0.0;
  G4double pzbil = 0.0;         

  std::cout <<"Number of particles: " << nopart << std::endl;
  
  for(G4int j = 0; j < nopart; j++) {
    varntp->itypcasc[j] = 1;
    // kind(): 1=proton, 2=neutron, 3=pi+, 4=pi0, 5=pi -      
    if(kind[j] == 1) { 
      varntp->avv[j] = 1;
      varntp->zvv[j] = 1;
      varntp->plab[j] = sqrt(ep[j]*(ep[j]+1876.5592)); // cugnon
      multp = multp + 1;
      mcorem = mcorem - ep[j] - 938.27231;
      std::cout <<"Proton produced! " << std::endl;
      std::cout <<"Momentum: "<< varntp->plab[j] << std::endl;
    }
//         if(kind(j).eq.2) then
//                 avv(j)=1
//                 zvv(j)=0
// c               plab(j)=sqrt(ep(j)*(ep(j)+1879.13126))     !ok  
//                     plab(j)=sqrt(ep(j)*(ep(j)+1876.5592))  !cugnon
//                 multn=multn+1
//                 mcorem=mcorem -ep(j) -939.56563
//         endif   
//         if(kind(j).eq.3) then
//                 avv(j)=-1
//                 zvv(j)=1
// c               plab(j)=sqrt(ep(j)*(ep(j)+279.1399))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+276.0))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -139.56995
//         endif
//         if(kind(j).eq.4) then
//                 avv(j)=-1
//                 zvv(j)=0
// c               plab(j)=sqrt(ep(j)*(ep(j)+269.9528))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+276.0))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -134.9764
//         endif
//         if(kind(j).eq.5) then
//                 avv(j)=-1
//                 zvv(j)=-1
// c               plab(j)=sqrt(ep(j)*(ep(j)+279.1399))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+276.0))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -139.56995
//         endif
//         if(kind(j).eq.6) then
//                 avv(j)=2
//                 zvv(j)=1
// c               plab(j)=sqrt(ep(j)*(ep(j)+279.1399))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+2.*1874.34))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -1874.34
//         endif
//         if(kind(j).eq.7) then
//                 avv(j)=3
//                 zvv(j)=1
// c               plab(j)=sqrt(ep(j)*(ep(j)+279.1399))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+2.*2806.359))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -2806.359
//         endif
//         if(kind(j).eq.8) then
//                 avv(j)=3
//                 zvv(j)=2
// c               plab(j)=sqrt(ep(j)*(ep(j)+279.1399))       !ok
//                     plab(j)=sqrt(ep(j)*(ep(j)+2.*2807.119))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -2807.119
//         endif
//         if(kind(j).eq.9) then
//                 avv(j)=4
//                 zvv(j)=2
// 			   plab(j)=sqrt(ep(j)*(ep(j)+2.*3724.818))      !cugnon
//                 numpi=numpi+1
//                 mcorem=mcorem -ep(j) -3724.818
//         endif
//         enerj(j)=ep(j)
//         tetlab(j)=180.*acos(gam(j))/3.141592654
//         philab(j)=180.*atan2(beta2(j),alpha(j))/3.141592654
        
//         pxbil = pxbil+plab(j)*alpha(j)
//         pybil = pybil+plab(j)*beta2(j)
//         pzbil = pzbil+plab(j)*gam(j)
        
// 900   continue
  }
}

// Init routines

void G4Incl::initIncl(G4bool initRandomSeed)
{
  // Subroutine for initialisation of intranuclear cascade incl
  //
  // this will  read some specific parameters for incl,
  // prepare the saxon-wood density for each nucleus
  // compute the deuteron momentum space density from paris pot.
  // print some global informations 
  //
  // input: should contain z and a for nbmat nucleus considered in this problem 
  //
  // input: should contain a seed (ial, odd and of 5 digits) to start the work.     

  G4double xrand;
  G4double ialdep;
  G4int imat;
  G4int iamat, izmat;

  // for the 19 secondary seeds of hazard:
  G4int nbtirhaz[IGRAINESIZE] = {38,82,76,18,39,31,41,59,26,54,
				 14,84,13,15,91,89,10,6,52};

  std::cout << " " << std::endl;	
  std::cout <<"*************** version incl 4.2 (C++)*********************" << std::endl;
  std::cout << "* stopping time and potential can be changed     ****" << std::endl;
  std::cout << "* input of first random numbers                  ****" << std::endl;
  std::cout << "* bimpact is output instead of sepa as in incl 3.0 **" << std::endl;
  std::cout << "* implementation of surface w.s.          4/2000  ***" << std::endl;
  std::cout << "* interaction only with \"participants\"    4/2000  ***" << std::endl;
  std::cout << "* cdpp:coherent dynamical pauli principle 5/2001  ***" << std::endl;
  std::cout <<"* paris momentum density for the deuteron 4/2001  ***" << std::endl;
  std::cout << "* nd-nn (*3) and phase space delta width  4/2001  ***" << std::endl;
  std::cout << "************* incl 4.0 -> incl4.1 *******************" << std::endl;
  std::cout <<"* no lower cut on pi-n interaction	   2/2002  ****" << std::endl;       
  std::cout <<"* init of first avatars for participants  2/2002  ***" << std::endl;     
  std::cout <<"* output of excit energy for absorption   2/2002  ***" << std::endl;
  std::cout <<"************* incl 4.1 -> incl4.2 *******************" << std::endl;
  std::cout <<"* increased absorption below 100 mev 	   5/2002  ****" << std::endl;       
  std::cout <<"* coulomb transmis. on projo from lahet  5/2002  ****" << std::endl;       
  std::cout <<"* light targets(gaussian a<5 mho 4<a<19) 6/2002  ****" << std::endl;       
  std::cout <<"*****************************************************" << std::endl;
  std::cout <<"* corrections 11/2003 numerical from j hendricks    *" << std::endl;
  std::cout <<"*           l of the remnant (from actual c.m.)     *" << std::endl;
  std::cout <<"*****************************************************" << std::endl;
  std::cout <<" " << std::endl;

  // specific parameters for incl:	 
  // espace de phases test (r et p) pour pauli: 
  // valeur recommandee par j.c. v-test=0.592 h**3:
  G4double rbl = 2.;
  G4double pbl = 200.;

  // valeur pour avoir v-test=2 h**3 (avec pbl=200)
  rbl = 3.1848;

  // constants and derived data
  G4double hc = 197.328;
  G4double fmp = 938.2796;
  G4double pf = 1.37*hc;  
  G4double tf = sqrt(pf*pf+fmp*fmp) - fmp;
  G4double binding = calincl->f[4] - tf; // Orig: f(5)
  std::cout <<"v0 nuclear potential (a>5) " << calincl->f[4] << "(mev) and scaling" << " stopping time factor: " << calincl->f[5] << std::endl; //Orig f(6)
  std::cout <<"fermi momentum: " << pf << " binding energy: " << binding << std::endl;
  std::cout <<"r*p cell for pauli stat: " << rbl << " " << pbl << std::endl;
  std::cout <<"pauli strict (1) or statistic (0): your choice: " << ws->npaulstr << std::endl;
  std::cout <<" " << std::endl;
  std::cout <<"your choice, nosurf=  " << ws->nosurf << "means" << std::endl;
  std::cout <<"nosurf=-2, with w.s. density and incl4 stopping time" << std::endl;
  std::cout <<"nosurf=1, sharp surface" << std::endl;
  std::cout <<"nosurf=0, with w.s. density, time without b dep" << std::endl;
  std::cout <<"nosurf=-1, with w.s. density, time with b dep." << std::endl;
  std::cout <<"rmaxws=ro+xfoisa*a, xfoisa= " << ws->xfoisa << " above a=19" << std::endl;

  G4int k1 = 5;
  G4int k2 = 0;
  G4int k3 = 0;
  G4int k4 = 2;
  G4int k5 = 1;
  G4int k6 = 0;

  std::cout << " " << std::endl;
  std::cout <<"k1,k2,k3,k4,k5,k6 " << k1 << " " << k2 << " " << k3 << " " << k4 << " " << k5 << " " << k6 << std::endl;
  std::cout <<"meaning: k1=5, reference frame = lab system" << std::endl;
  std::cout <<"k2=0, relativistic kinematics" << std::endl;
  std::cout <<"k3=0, deltas are produced" << std::endl;
  std::cout <<"k4=2, delta has a exponentially random lifetime" << std::endl;
  std::cout <<"k5=1, delta-nucleon=delta-delta = nucleon-nucleon elastic x-section" << std::endl;
  std::cout <<" " << std::endl;

  // preparation of 19 other seeds (can also be initialized from outside):
  if(initRandomSeed) {
    ialdep=hazard->ial;
    for(int i = 0; i < IGRAINESIZE; i++) {
      for(int j = 0; j < nbtirhaz[i]; j++) {
	standardRandom(&xrand,&(hazard->ial));
      }

      // Zero is not accepted as random seed!
      do {
	standardRandom(&xrand,&(hazard->ial));
      } while(xrand == 0);

      xrand = xrand*100000;

      while(xrand < 10000) {
	xrand = xrand*10;
      }
      hazard->igraine[i] = (int) xrand;
      if(hazard->igraine[i] == ((hazard->igraine[i]/2)*2)) {
	hazard->igraine[i] = hazard->igraine[i] + 1;
      }
    }

    hazard->ial= int(ialdep);
  }

  // calculation with realistic nuclear density (saxon-wood)
  if (ws->nosurf <= 0) {
    // prepare nucleus density for nbmat nucleus defined in struct mat
    if(mat->nbmat >= 500) {
      std::cout <<"You need " << mat->nbmat << " nuclei in your problem. The maximum number of nuclei is 500 " << std::endl;
      return;
    }

    for(G4int i = 0; i < mat->nbmat; i++) {
      imat=i;
      izmat = int(mat->zmat[i]);
      iamat = int(mat->amat[i]);
      std::cout <<"Initializing material with values:" << std::endl;
      std::cout <<"imat = " << i << std::endl;
      std::cout <<"izmat = " << izmat << std::endl;
      std::cout <<"iamat = " << iamat << std::endl;

      initMaterial(izmat, iamat, imat);
    }
  }

  // deuteron density in momentum space:
  densDeut();
}


void G4Incl::initMaterial(G4int izmat, G4int iamat, G4int imat)
{
  G4double res_dws;
  G4double fnor;

  G4double rcour, geom;
  G4int nbr;

  G4double step, f_r;

  // rms espace r, espace p, fermi momentum and energy for light gauss nuc.      
  const G4double datarms1t[LGNSIZE] = {0.0, 0.0, 0.0, 0.0, 0.0, 2.10, 1.80, 1.80, 1.63};
  const G4double datapf1t[LGNSIZE] = {0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 110.0, 110.0, 153.0};

  for(G4int i = 0; i < LGNSIZE; i++) {
    light_gaus_nuc->rms1t[i] = datarms1t[i];
    light_gaus_nuc->pf1t[i] = datapf1t[i];
  }

  // fermi 2 param from a=19 to 28, modified harm oscil a=6 to 18
  // (h. de vries et al. at. data and nuc. data tab. 36 (1987) 495)
  const G4double datarln[LNSIZE] = {0.,0.,0.,0.,0.,0.334,0.327,0.479,0.631,0.838,
				    0.811,1.07,1.403,1.335,1.25,1.544,1.498,1.513,
				    2.58,2.77, 2.775,2.78,2.88,2.98,3.22,3.03,2.84,
				    3.14,0.,0.};

  const G4double dataaln[LNSIZE] = {0.,0.,0.,0.,0.,1.78,1.77,1.77,1.77,1.71,
				    1.69,1.69,1.635,1.730,1.81,1.833,1.798,
				    1.841,0.567,0.571, 0.560,0.549,0.550,0.551,
				    0.580,0.575,0.569,0.537,0.,0.};

  for(G4int i = 0; i < LNSIZE; i++) {
    light_nuc->r[i] = datarln[i];
    light_nuc->a[i] = dataaln[i];
  }

  std::cout <<"***************************************************" << std::endl;
  std::cout <<"Nuclear density for nucleus (z, a): " << izmat << " " << iamat << " " << imat << std::endl;

  const G4double fmp = 938.2796;

  // parametres moyens de densite de la cible (fermi 2 parametres)
  if (iamat >= 28) {
    ws->r0 = (2.745e-4*iamat+1.063)*pow(iamat, 0.33333333);
    ws->adif = 1.63e-4*iamat+0.510;
    ws->rmaxws = ws->r0 + (ws->xfoisa)*(ws->adif);
  }
  else if(iamat >= 19) {
    ws->r0 = light_nuc->r[iamat];
    ws->adif = light_nuc->a[iamat];
    ws->rmaxws = ws->r0 + (ws->xfoisa)*(ws->adif);
  }
  else if(iamat >= 6) {
    ws->r0 = light_nuc->r[iamat];
    ws->adif = light_nuc->a[iamat];
    ws->rmaxws = 5.5 + 0.3*(iamat-6.)/12.;
  }
  else if(iamat >= 2) {
    if(iamat == 2) {
      ws->r0=light_gaus_nuc->rms1t[5]; // Orig: rms1t(6)
      light_gaus_nuc->pfln[5] = light_gaus_nuc->pf1t[5]*1.291;  // Orig [6], sqrt(5/3)=1.291
      light_gaus_nuc->tfln[5] = sqrt(pow(light_gaus_nuc->pfln[5],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[5] = light_gaus_nuc->tfln[5] + 2.22;
      std::cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[5] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[5] << " " << light_gaus_nuc->tfln[5] << std::endl;
    }
    if(iamat == 3 && izmat == 1) {
      ws->r0=light_gaus_nuc->rms1t[6]; // Orig: rms1t(7)
      light_gaus_nuc->pfln[6] = light_gaus_nuc->pf1t[6]*1.291;  // Orig [7], sqrt(5/3)=1.291
      light_gaus_nuc->tfln[6] = sqrt(pow(light_gaus_nuc->pfln[6],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[6] = light_gaus_nuc->tfln[6] + 4.24;
      std::cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[6] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[6] << " " << light_gaus_nuc->tfln[6] << std::endl;
    }
    if(iamat == 3 && izmat == 2) {
      ws->r0 = light_gaus_nuc->rms1t[7]; // Orig: rms1t(8)
      light_gaus_nuc->pfln[7] = light_gaus_nuc->pf1t[7]*1.291;   //!sqrt(5/3)=1.291
      light_gaus_nuc->tfln[7] = sqrt(pow(light_gaus_nuc->pfln[7],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[7] = light_gaus_nuc->tfln[7] + 3.86;
      std::cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[7] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[7] << " " << light_gaus_nuc->tfln[7] << std::endl;
    }
    if(iamat == 4) {
      ws->r0 = light_gaus_nuc->rms1t[8]; // Orig: rms1t(9)
      light_gaus_nuc->pfln[8] = light_gaus_nuc->pf1t[8]*1.291;   // !sqrt(5/3)=1.291
      light_gaus_nuc->tfln[8] = sqrt(pow(light_gaus_nuc->pfln[8],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[8] = light_gaus_nuc->tfln[8] + 9.43;
      std::cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[8] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[8] << " " << light_gaus_nuc->tfln[8] << std::endl;
    }
    ws->adif = 0.57735*ws->r0;
    ws->rmaxws = ws->r0 + 2.5;
  }
  ws->drws = (ws->rmaxws)/29.0;

  // bmax for sigma geom and various projectiles (p,n,pion/d/t/he3/he4/)
  G4int j;
  for(G4int i = 0; i < MATGEOSIZE; i++) { // Orig: do i=1,6
    j = i;
    if(i >= 2) {
      j = i + 3;
    }
    mat->bmax_geo[i][imat] = (ws->rmaxws) + (light_gaus_nuc->rms1t[j]);
  }

  // preparation de la distribution w.s.:
  if (iamat >= 19) {
    G4double step = 0.2;
    res_dws = integrate(0.0, 13.5, step, derivWsaxFunction);
  }
  else { 
    // preparation de la distribution m.h.o.:
    if(iamat >= 6) {
      step=0.1;
      res_dws = integrate(0.0, 10.0, step, derivMhoFunction);
    }
    else {
      // preparation de la distribution gaussienne:
      //	 G4double cte = pow(ws->adif,3)*sqrt(2.*3.141592654);        
      res_dws = 3.0*(pow(ws->adif, 3)*sqrt(2.0*3.141592654))/2.0;
    }
  }
  fnor = res_dws;

  // calcul de q/pf=f(r)      
  nbr = int(floor((ws->rmaxws)/(ws->drws) + 1.5));
  rcour = -1*(ws->drws);

  j = 0;
  for(G4int i = 0; i < nbr; i++) { // do i=1,nbr
    rcour = rcour + (ws->drws);
    if(i == 0) { // 1->0
      f_r = 0.0;
      saxw->x[j][imat] = f_r;
      saxw->y[j][imat] = 0.0;		//!on impose x(1)=0., y(1)=0.
      res_dws = 0.0;
    }
    else {
      step = rcour/20.;
      if(step >= 0.05) {
	step = 0.05;
      }
      if (iamat >= 19) {
	//integ(ws, dton, 0.,rcour,step,&derivwsax,&res_dws);
	res_dws = integrate(0.0, rcour, step, derivWsaxFunction);
	f_r = res_dws/fnor;
      }
      else { 
	if(iamat >= 6) {
	  //integ(ws, dton, 0.,rcour,step,&derivmho,&res_dws);
	  res_dws = integrate(0.0, rcour, step, derivMhoFunction);
	  f_r = res_dws/fnor;
	}
	else { 
	  //integ(ws, dton, 0.,rcour,step,&derivgaus,&res_dws);
	  res_dws = integrate(0.0, rcour, step, derivGausFunction);
	  f_r = res_dws/fnor;
	}
      }
      // modif le 20/10/2003; éviter les valeurs négatives avant **1/3 !
      //       }
      if(f_r >= 0.0)  {
	f_r = pow(f_r,(1./3.));
	std::cout <<"j = " << j << std::endl;
	saxw->x[j][imat] = f_r;
	saxw->y[j][imat] = rcour;
      }
    }
    j = j + 1;
  }
  saxw->n = j;
  saxw->x[j-1][imat] = 1.; // !on impose saxw->x[nbpinter-1]=1. (y=rmax)

  // interpolation de f_inv(r) (fonction inverse de f(r))       
  // flin2(imat, saxw, ws);
  firstDerivative(imat);
     
  if(iamat >= 19) {
    std::cout <<"Wood-Saxon density, r0 = " << ws->r0 << " a = " << ws->adif << std::endl;
  }
  if(iamat >= 6 && iamat <= 19) {
    std::cout <<"Modif. harm. oscil. density, alpha = " << ws->r0 << " a = " << ws->adif << std::endl;
  }
  if(iamat >= 2 && iamat <= 6) {
    std::cout <<"Gaussian density, r.m.s = " << ws->r0 << " sigma = " << ws->adif << std::endl;
  }
  geom = 31.41592653*pow(ws->rmaxws,2);
  std::cout <<"For incident nucleons or pions rmax = " << ws->rmaxws << " and geometrical (pi*rmaxws*rmaxws) reaction cross section (mb) is " << geom << std::endl;
  for(G4int k = 2; k < MATGEOSIZE; k++) {
    std::cout << "Rmaxws for d/t/3he/4he = " << mat->bmax_geo[k][imat] << std::endl;
  }

  std::cout <<"Exact calculation of the r(q) function for the target nucleus density q/pf  r(q/pf)" << std::endl;
  std::cout <<"*********************************************" << std::endl;
}

G4double G4Incl::deutv(G4int l, G4double q)
{
  //
  
  G4double res = 0.0;

  if (l == 0) {
    for(G4int i = 0; i < DTONSIZE; i++) {
      res = res + dton->c[i]/(pow(q,2) + fm2(i+1));
    }
  }
  if(l != 0) {
    for(G4int i = 0; i < DTONSIZE; i++) { 
      res = res + dton->d[i]/(pow(q,2) + fm2(i+1));
    }
  }

  return res*sqrt(2./M_PI)*dton->fn; // See G4InclDataDefs.hh
}

G4double G4Incl::fm2(G4int j)
{
  // Fm2
  
  return pow((0.23162461 + (j - 1)),2);
}

G4double G4Incl::interpolateFunction(G4double xv)
{
  // fonction d'interpolation au point xv ( meme hors bornes )             
  // de la fn x->y dont les derivees premieres (s) ont ete                 
  // evaluees par l'appel prealable de flin2                              
  // les indices vont de 1 a n

  G4double tz = xv - saxw->x[0][saxw->k];

  if(tz < 0) {
    return (saxw->y[0][saxw->k] + saxw->s[0][saxw->k]*tz);
  }
  else if(tz == 0) {
    return (saxw->y[0][saxw->k]);
  }
  else {
    G4int j;
    for(G4int i = 1; i < saxw->n; i++) {
      j = i - 1;
      tz = xv - saxw->x[j][saxw->k];
      if(tz < 0) {
	G4double dgx = xv - saxw->x[j][saxw->k];
	return(saxw->y[j][saxw->k] + saxw->s[j][saxw->k]*dgx);
      }
      else if(tz == 0) {
	return saxw->y[j][saxw->k];
      }
    }
  }

  return 0.0;
}

void G4Incl::firstDerivative(G4int k)
{
  for(G4int i=0; i < saxw->n-1; i++) {
    saxw->s[i][k] = (saxw->y[i+1][k] - saxw->y[i][k]) / (saxw->x[i+1][k] - saxw->x[i][k]);
  }
  saxw->s[saxw->n-1][k] = saxw->s[saxw->n-2][k];
}

G4double G4Incl::wsax(G4double r) {
  return pow(r,2) / (1.0+exp(r-(ws->r0)/(ws->adif)));
}

G4double G4Incl::derivWsax(G4double r)
{
  G4double derivwsax = pow(r,3)*exp((r-(ws->r0))/(ws->adif))/pow((1.0+exp((r-(ws->r0))/(ws->adif))),2);
  return derivwsax/(ws->adif);
}

G4double G4Incl::dmho(G4double r)
{
  G4double arg=pow((r/(ws->adif)),2);
  return r*r*(1.+(ws->r0)*arg)*exp(-arg);
}

G4double G4Incl::derivMho(G4double r)
{
  G4double arg=pow((r/(ws->adif)),2);
  return -2.*r*r*arg*((ws->r0) -1.-(ws->r0)*arg)*exp(-arg);
}

G4double G4Incl::derivGaus(G4double r)
{
  G4double arg=pow((r/(ws->adif)),2);
  return r*r*arg*exp(-arg/2.);      
}

void G4Incl::densDeut()
{
  // ce subroutine appele sur le premier tir va calculer la densite du deuton
  // dans l'espace des impulsions et preparer l'interpolation permettant ensuite
  // le tir au hasard d'un module de l'impulsion (q).
  // ce subroutine remplit le common /spl2/:
  // xsp(0:1), ysp integrale normalisee de la densite de 0 a q.
  // a(),b(),c() coefs des nsp points pour une interpolation du second degre.
  // q est en fm-1. 

  //    495	      dimension q(100),f(100)
  //    496	      common/spl2/ xsp(100),ysp(100),a(100),b(100),cc(100),nbp
  G4double cData[DTONSIZE] = {0.88688076e+00,-0.34717093e+00,-.30502380e+01,
			      .56207766e+02,-.74957334e+03,.53365279e+04,-.22706863e+05,
			      .60434469e+05,-.10292058e+06,.11223357e+06,-.75925226e+05,
			      .29059715e+05,-.48157368e+04};

  G4double dData[DTONSIZE] = {.23135193e-01,-.85604572e+00,.56068193e+01,
			      -.69462922e+02,.41631118e+03,-.12546621e+04,.12387830e+04,
			      .33739172e+04,-.13041151e+05,.19512524e+05,-.15634324e+05,
			      .66231089e+04,-.11698185e+04};

  G4double fnData = 0.28212e+00;

  for(G4int i = 0; i < DTONSIZE; i++) {
    dton->c[i] = cData[i];
    dton->d[i] = dData[i];
  }
  dton->fn = fnData;

  //    509	c avec fn=.28212 les fo radiales suivantes sont normalisees a:          deu00470
  //    510	c somme(0,infini)(deut0(q)**2 + deut2(q)**2))*q*q*dq = 1./4*pi          deu00480
  //    511	c et ceci dans l'espace r et dans l'espace q. pd=5.74%                  deu00490
  //    512	cjcd
  //    513	      common /inout/ in, io, itty, iscrt                                
  //    514	cjcd

  const G4int qsize = 100;
  G4double q[qsize];
  G4double f[qsize];
  G4double dq=0.01;
  q[0]=0.0;
  for(G4int i = 1; i < 50; i++) {
    q[i] = q[i-1] + dq;
    f[i] = 0.0;
  }

  spl2->n = 77; // nombre de points de calcul

  dq=0.1;
  for(G4int i = 50; i < spl2->n; i++) {
    q[i] = q[i-1] + dq;
  }

  f[0]=0.0;

  G4double sumint=0.0;

  // the id if the function we wish to integrate (in this case: G4Incl::dens
  for(G4int i = 1; i < spl2->n; i++) {
    dq = (q[i]-q[i-1])/10.0;
    sumint = sumint + integrate(q[i-1], q[i], dq, densFunction);
    f[i] = sumint;
  }

  for(G4int i = 0; i < spl2->n; i++) {
    spl2->x[i] = f[i]/f[spl2->n-1];
    spl2->y[i] = q[i];
  }

  spl2ab();

  std::cout << "deuteron density in q space from Paris potential: " << spl2->n << " Exact values from 0 to " 
	    << q[spl2->n-1] << " fm-1 " << std::endl;
}

G4double G4Incl::integrate(G4double ami, G4double ama, G4double step, G4int functionChoice)
{
  G4double res;
  G4double x1[5];
  G4double ri = ami;
  G4double ra = ama;
  G4double nb;
  G4double acont = 1.0;

  if(ama <= ami) {
    acont = -1.0;
    ri = ama;
    ra = ami;
  }
  
  x1[0] = 95.0/288.0;
  x1[1] = 317.0/240.0;
  x1[2] = 23.0/30.0;
  x1[3] = 793.0/720.0;
  x1[4] = 157.0/160.0;
  nb = ((ra - ri)/step) + 1.0000000001;
  G4double dr = (ra - ri)/(nb - 1);
  res = 0.0;

  if(nb >= 10) {
    for(G4int i=0; i < 5; i++) {
      res = res + (callFunction(functionChoice, ri) + callFunction(functionChoice, ra))*x1[i];
      ri = ri + dr;
      ra = ra - dr;
    }
  }
  else {
    std::cout <<"pas assez de points d|integration" << std::endl;
    return (res*step*acont);
  }
  nb = nb - 10;
  if(nb != 0) {
    for(G4int i = 1; i < nb; i++) {
      res = res + callFunction(functionChoice, ri);
      ri = ri + dr;
    }
  }

  return (res*step*acont);
}


G4double G4Incl::dens(G4double q)
{
  return q*q*(pow(deutv(0,q),2)+pow(deutv(2,q),2));
}

void G4Incl::spl2ab() 
{
  G4int i, j, k;

  for(i=0; i <= spl2->n-3; i++) {
    j = i + 1;
    k = i + 2;

    spl2->c[i] = ((spl2->y[k]-spl2->y[i])*(spl2->x[j]-spl2->x[i])-(spl2->x[k]-spl2->x[i])*(spl2->y[j]-spl2->y[i]))
      /((spl2->x[j]-spl2->x[i])*(spl2->x[k]-spl2->x[i])*(spl2->x[k]-spl2->x[j]));

    spl2->b[i] = (spl2->y[j]-spl2->y[i])/(spl2->x[j]-spl2->x[i]);

    spl2->a[i] = spl2->y[i];
  }

  for(i = spl2->n-2; i < spl2->n; i++) {
    spl2->c[i] = spl2->c[spl2->n-3];
    spl2->b[i] = spl2->b[spl2->n-3];
    spl2->a[i] = spl2->a[spl2->n-3]; 
  }
}

G4double G4Incl::splineab(G4double xv)
{
  G4double tz;
  G4int j;
  
  tz = xv-spl2->x[0];

  if(tz < 0) {
    return spl2->a[0] + spl2->b[0] * tz + spl2->c[0] * tz * (xv - spl2->x[1]);
  }
  if(tz == 0) {
    return spl2->y[0];                                                       
  }
  if(tz > 0) {
    for(G4int i = 1; i <= spl2->n-1; i++) {
      j = i;
      tz = xv - spl2->x[i];                                  

      if(tz < 0) { 
	j = j - 1;
	tz = xv - spl2->x[j];
	return spl2->a[j] + spl2->b[j] * tz + spl2->c[j] * tz * (xv - spl2->x[j+1]);
      }
      if(tz == 0) {
	return spl2->y[j];                                                    
      }
    }
  }

  // Returns 0.0 if the point xv is outside the defined region (xv > spl2->x[spl2->n-1])
  if(verboseLevel > 0) {
    // G4cout <<"G4Incl::splineab : requested point outside defined region! Returning 0.0." << G4endl;
  }
  return 0.0; 
}

// Actual calculation

void G4Incl::pnu(G4int *ibert_p, G4int *nopart_p, G4int *izrem_p, G4int *iarem_p, G4double *esrem_p,
		 G4double *erecrem_p, G4double *alrem_p, G4double *berem_p, G4double *garem_p,
		 G4double *bimpact_p, G4int *l_p)
{
  G4int ibert = (*ibert_p);
  //  float f[15]; // = (*f_p);
  G4int nopart = (*nopart_p);
//   G4int kind[300]; //= (*kind_p);
//   G4double ep[300]; // = (*ep_p);
//   G4double alpha[300]; // = (*alpha_p); 
//   G4double beta[300]; // = (*beta_p);
//   G4double gam[300]; // = (*gam_p);
  G4int izrem = (*izrem_p);
  G4int iarem = (*iarem_p);
  G4double esrem = (*esrem_p); 
  G4double erecrem = (*erecrem_p);
  G4double alrem = (*alrem_p);
  G4double berem = (*berem_p);
  G4double garem = (*garem_p);
  G4double bimpact = (*bimpact_p);
  G4int l = (*l_p);

  //alog 
  G4double aml1;
  G4double aml2; 
  G4double amlnew; 
  G4double arg; 
  G4double b1;
  G4double b2; 
  G4double b3;
  G4double bb2; 
  G4double be;
  G4double bmass[2000]; 
  G4double bmax2; 
  G4double c1; 
  G4double c2; 
  G4double cb0; 
  G4double cchi; 
  G4double ccr; 
  G4double cg; 
  G4double cif; 
  G4double cmultn; 
  G4double cobe; 
  G4double coeffb0; 
  G4double comom;
  G4double cstet; 
  G4double dis1; 
  G4double dis2; 
  G4double dis3; 
  G4double dist; 
  G4double eb0; 
  G4double ecoreh5; 
  G4double efer; 
  G4double egs; 
  G4double eh5; 
  G4double eh6; 
  G4double eij; 
  G4double ekout; 
  G4double elead; 
  G4double energie_in; 
  G4double ener_max; 
  G4double eout; 
  G4double eps_c[BL1SIZE]; 
  G4double epsv; 
  G4double erecg; 
  G4double erem; 
  G4double exi; 
  G4double expob0; 
  G4double factemp; 
  G4double fffc; 
  G4double fm; 
  G4double g1; 
  G4double g2; 
  G4double ge; 
  G4double geff; 
  G4double gg; 
  G4double gl1; 
  G4double gl2; 
  G4double gpsg; 
  G4int i1; 
  G4int i20; 
  G4int ic33; 
  G4int ich1; 
  G4int ich2; 
  G4int ich3; 
  G4int ich4; 
  G4int ichd; 
  G4int ichpion; 
  G4int idecf; 
  G4int idep; 
  G4int iej; 
  G4int iejn; 
  G4int iejp; 
  G4int i_emax; 
  G4int iflag; 
  G4int iflag20; 
  G4int iflag40; 
  G4int iflag60; 
  G4int ilm; 
  G4int imin; 
  G4int indic[2000]; 
  G4int inrem; 
  G4int ip; 
  G4int ipi[2000]; 
  G4int iqe; 
  G4int irem; 
  G4int irst_avatar; 
  G4int isos; 
  G4int itch; 
  G4int iteste; 
  G4int itt; 
  G4int ixr1; 
  G4int ixr2; 
  G4int ixr3; 
  G4int k; 
  G4int kcol; 
  G4int kd; 
  G4int klm; 
  G4int l1; 
  G4int l2; 
  G4int ldel; 
  G4int lead; 
  G4int led; 
  G4int lnew; 
  G4int lp; 
  G4int lp1; 
  G4double mcdd; 
  //G4double mg; 
  G4int mg;
  G4double mpaul1; 
  G4double mpaul2; 
  G4double mrdd; 
  G4double mrdn; 
  G4double mrdp; 
  G4double mrnd; 
  G4double mrnn; 
  G4double mrpd; 
  G4int n20; 
  G4int nbalttf; 
  G4int nbquit; 
  G4int nbtest; 
  G4int nc[300]; 
  G4int ncol; 
  G4int ncol_2c; 
  G4int next; 
  G4int nmiss; 
  G4int np; 
  G4int npidir; 
  G4int npion; 
  G4int npproj[300]; 
  G4int npx; 
  G4int nsum_col; 
  G4double p1v; 
  G4double p2v; 
  G4double p3_c[BL1SIZE]; 
  G4double p3v; 
  G4double pfrem1; 
  G4double pfrem2; 
  G4double pfrem3; 
  G4double pfreml; 
  G4double pfreml2; 
  G4double phi; 
  G4double p_mod; 
  G4double pot; 
  G4double pout1; 
  G4double pout2; 
  G4double pout3; 
  G4double pppp; 
  G4double prem1; 
  G4double prem2; 
  G4double prem3; 
  G4double psf; 
  G4double pspr; 
  G4double ptotl; 
  G4double q1[BL1SIZE]; 
  G4double q2[BL1SIZE]; 
  G4double q3[BL1SIZE]; 
  G4double q4[BL1SIZE]; 
  G4double qdeut; 
  G4double qqq; 
  G4double r22; 
  G4double rcm1; 
  G4double rcm2; 
  G4double rcm3; 
  G4double rcorr; 
  G4double rhopi; 
  G4double rndm; 
  G4double rr; 
  G4double rrrr; 
  G4double s; 
  G4double s1t1; 
  G4double s2t1; 
  G4double s3t1; 
  G4double schi; 
  G4double sepa; 
  G4double sif; 
  G4double sitet; 
  G4double sp1t1; 
  G4double sp2t1; 
  G4double sp3t1; 
  G4double sq; 
  G4double sueps; 
  G4double t[50]; 
  G4double t0; 
  G4double t1; 
  G4double t2; 
  G4double t3; 
  G4double t33; 
  G4double t4; 
  G4double t5; 
  G4double t6; 
  G4double t7; 
  G4double t8; 
  G4double tau; 
  G4double tbid; 
  G4double tdel; 
  G4double temfin; 
  G4double tim; 
  G4double timi; 
  G4double tlabu; 
  G4double tp; 
  G4double tref; 
  G4double tri; 
  G4double tt31; 
  G4double tt32; 
  G4double tt33; 
  G4double tt34; 
  G4double tt35; 
  G4double tt36; 
  G4double tte; 
  G4double u; 
  G4double v; 
  G4double var_ab; 
  G4double x; 
  G4double x1l1; 
  G4double x1l2; 
  G4double x1_target; 
  G4double x2cour; 
  G4double x2l1; 
  G4double x2l2; 
  G4double x2_target; 
  //G4double x3; 
  G4double x3l1; 
  G4double x3l2; 
  G4double x3_target; 
  G4double xapres; 
  G4double xavant; 
  G4double xbl1; 
  G4double xbl2; 
  G4double xc; 
  G4double xe; 
  G4double xga; 
  G4double xl1; 
  G4double xl2; 
  G4double xl3; 
  G4double xlab; 
  G4double xleng; 
  G4double xlengm; 
  G4double xmodp; 
  G4double xpb; 
  G4double xq; 
  G4double xr1; 
  G4double xr2; 
  G4double xr3; 
  G4double xr4; 
  G4double xr5; 
  G4double xr6; 
  G4double xr7; 
  G4double xr8; 
  G4double xv; 
  G4double xxx; 
  G4double xy1; 
  G4double xy2; 
  G4double xy3; 
  G4double xye; 
  G4double y; 
  G4double y1[BL3SIZE];
  G4double y2[BL3SIZE]; 
  G4double y3[BL3SIZE]; 
  G4double ym[2000]; 
  G4double z; 
  G4double za_i; 
  G4double zai2; 
  G4double zshif; 
  G4double ztouch; 
  G4double ztu; 

  // LIEGE INC-model as a subroutine 

  // The Liege INC model has been applied to several systems and in
  // several conditions. Refinements are still in progress

  // PLEASE refer to this version as INCL4.1 in order to avoid
  // confusion when comparing to the results of your colleagues.

  // DIFFERENT from INCL2.0 in the sense that the cascade is stopped
  // when the excitation energy vanishes if this occurs before the
  // predetermined stopping time (herein denoted as temfin) Special
  // are taken to avoid emission of slow particles due to the
  // imperfect Pauli blocking

  // PLEASE notice: There are basically only two parameters in the
  // model: the average potential depth, denoted as V0, and the time
  // at which the cascade is stopped denoted as temfin. In this
  // program, the "standard" values (those G4introduced in the ref
  // NPA620(1997)475) are V0=40MeV and temfin=1.25*some function of
  // incident energy and impact parameter. You may, of course, change
  // these parameters V0 and the numerical coefficient in temfin,
  // within reasonable limits (i.e. V0 cannot be lower than 38MeV;
  // V0=45MeV is recommended for heavy nuclei). If you do, PLEASE
  // indicate your choice, once again, for the same reason as above.

  // The description of the cascade model(incl2.0) can be found in:
  // J.C., C.VOLANT & S.VUILLIER, NPA620(1997)475 It is basically the
  // same model as described in J.C. NPA462(1987)751 (version 7 (in
  // our jargon), sketched below) + a refinement of the
  // parametrization of the cross-sections, based on J.C., D. L'HOTE,
  // J.VANDERMEULEN, NIM B111(1996)215 and J.C., S.LERAY, E.MARTINEZ,
  // Y.PATIN & S.VUILLIER PRC56(1998)2431

  // technical notes: 
  // 1.for the parametrizations of cross  sections, see
  //   notes of 4/10/96, 9/10/96, 31/12/97 and 13/10/98
  // 2.temfin=1.25*... 18/6/98
  // 3.sepa in concordance with v0-tf 2/7/98
  // 4.special care for stopping the cascade before t=temfin 27/04/99

  //     84	c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //     85	C                                                                       P-N00030
  //     86	C     VERSION 7:  2/2/93                                                P-N00040
  //     87	C
  //     88	C++++++++++++ DESCRIPTION OF INPUT AND OUTPUT+++++++++++++++++++++++++++++
  //     89	C
  //     90	C **INPUT DATA**
  //     91	C
  //     92	C  IBERT=O IN THE FIRST CALL
  //     93	C        1 IN THE SUBSEQUENT CALLS
  //     94	C
  //     95	C  F= (REAL) ARRAY OF DIMENSION 8
  //     96	C
  //     97	C    F(1)= A (TARGET)
  //     98	C    F(2)= Z (TARGET)
  //     99	C    F(3)= KINETIC ENERGY (IN MEV) OF THE INCIDENT PARTICLE
  //    100	C    F(4)= SUPPOSED TO BE THE MINIMUM  PROTON ENERGY REQUIRED TO LEAVE
  //    101	C          THE TARGET. IN THIS CASCADE MODEL, IT IS ZERO
  //    102	C    F(5)= Nuclear potential V0 (standard value=45 MeV for heavy nuclei)
  //    103	C    F(6)= Rescale the cascade duration (the standard value t0 is MULTIPLIED
  //    104	c 	   by this value. F(6)=1. is the standard)
  //    105	C    F(7)= TYPE OF INCIDENT PARTICLE
  //    106	C           1.00 FOR PROTON
  //    107	C           2.00 FOR NEUTRON
  //    108	C           3.00 FOR PI+
  //    109	C           4.00 FOR PI0
  //    110	C           5.00 FOR PI-
  //    111	C           6.00 FOR DEUTERON
  //    112	C           7.00 FOR TRITON
  //    113	C           8.00 FOR HE3
  //    114	C           9.00 FOR HE4
  //    115	C    F(8)= SUPPOSED TO BE THE MINIMUM NEUTRON ENERGY REQUIRED TO LEAVE
  //    116	C          THE TARGET. IN THIS CASCADE MODEL, IT IS ZERO
  //    117	C
  //    118	C                                 NOSURF=1 Sharp density (hard sphere), 
  //    119	C                                 NOSURF=0 Wood-Saxon density, stopping time "70" 
  //    120	C                                                without B (impact) dependence.
  //    121	C                                 NOSURF=-1 Wood-Saxon density, stopping time "70" 
  //    122	C                                                with B (impact) dependence
  //    123	C                                 (on peut toujours nenormaliser ces fonctions 
  //    124	C                                  de temps avec le facteur F(6): t=t0*F(6) ) 
  //    125	C                                XFOISA      Rmaxws = R0 + XFOISA*A
  //    126	C					     Bmax = Rmaxws for pions and nucleons
  //    127	C					     Bmax = Rmaxws + rms1t (data) for composits
  //    128	C         Pauli strict (1) or statistic (0) or without pauli (2):     NPAULSTR
  //    129	C
  //    130	C    F(9)= imat, target material identifier for the right choose of Sax.-Wood density
  //    131	c
  //    132	c common/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
  //    133	c     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
  //    134	c ......20 numbers in the main routine to initialize the random numbers
  //    135	c
  //    136	C **OUTPUT DATA**
  //    137	C
  //    138	C  NOPART=-1 PSEUDO REACTION (VOID EVENT)
  //    139	C          0 ABSORPTION
  //    140	C         >0 TRUE EVENT, = NUMBER OF PARTICLES EMITTED (EXCLUDING THE REMNANT)
  //    141	C
  //    142	C  FOR N=1,NOPART:
  //    143	C  KIND(N)= TYPE OF PARTICLES (SAME CONVENTION AS FOR F(7), BUT IN G4INTEGERS)
  //    144	C
  //    145	C  EP(N)=  KINETIC ENERGY
  //    146	C
  //    147	C  ALPHA(N),BETA(N),GAM(N)= DIRECTION COSINES
  //    148	C
  //    149	C  IZREM= Z (REMNANT)
  //    150	C
  //    151	C  IAREM= A (REMNANT)
  //    152	C
  //    153	C  ESREM= EXCITATION ENERGY OF THE REMNANT
  //    154	C
  //    155	C  ERECREM= RECOIL ENERGY OF THE REMNANT
  //    156	C
  //    157	C  ALREM,BEREM,GAREM=DIRECTION COSINES OF THE REMNANT
  //    158	C
  //    159	C  BIMPACT impact parameter
  //    160	C
  //    161	C  L G4intrinsic momentum of the remnant in units of h/2pi=hbar=197.328
  //    162	C+++++++++ DESCRIPTION OF THE INC MODEL ++++++++++++++++++++++++++++++++++++
  //    163	C
  //    164	C     MODEL DESCRIBED IN J.CUGNON (NP A462(1987)751)                    P-N00050
  //    165	C     =MODEL (DR) OF J.CUGNON,D.KINET,J.VANDERMEULEN(NP A379(1982)567)  P-N00060
  //    166	C                                                                       P-N00110
  //    167	C        +REFLECTION OR TRANSMISSION ON THE POTENTIAL WALL              P-N00120
  //    168	C              (THE POTENTIAL DEPTH IS THE SAME FOR NUCLEONS & DELTA'S) P-N00130
  //    169	C              (CONTAINS A COULOMB BARRIER)                             P-N00140
  //    170	C                                                                       P-N00150
  //    171	C        +ABSORPTION OF THE PION ABOVE THE (3,3) RESONANCE (NOT IN      P-N00160
  //    172	C                VERSION 2)                                             P-N00170
  //    173	C                                                                       P-N00180
  //    174	C        +POSSIBLE PAULI BLOCKING OF TWO BODY COLLISIONS                P-N00190
  //    175	C        +POSSIBLE PAULI BLOCKING OF DELTA DECAY                        P-N00200
  //    176	C                 THE PAULI BLOCKING IS APPLIED TO THE NUCLEONS         P-N00210
  //    177	C                 ONLY.THE PAULI BLOCKING FACTORS ARE EVALUATED         P-N00220
  //    178	C                 BY COUNTING THE NUCLEONS INSIDE A VOLUME IN           P-N00230
  //    179	C                 PHASE SPACE.THE EXTENSION OF THIS VOLUME IS           P-N00240
  //    180	C                 OF THE ORDER OF H**3                                  P-N00250
  //    181	C                                                                       P-N00260
  //    182	C    ADDITIONAL FEATURES:                                               P-N00270
  //    183	C                                                                       P-N00280
  //    184	C        +ISOSPIN (WITH NEW PN BACKWARD-FORWARD ASYMMETRY)              P-N00290
  //    185	C                                                                       P-N00300
  //    186	C        +"LONGITUDINAL GROWTH" OF THE BARYONS (NOT ACTIVATED HERE)
  //    187	C
  //    188	C        + PARTICLE #1 IS ALWAYS THE FASTEST PARTICLE IN THE Z-DIRECTIONP-N00100
  //    189	C                                 (NOT ACTIVATED HERE)
  //    190	C        +SIMPLIFIED NEUTRON EVAPORATION AT THE END OF THE CASCADE      P-N00310
  //    191	C                                 (NOT PRESENT HERE)
  //    192	C                 
  //    193	C        +POSSIBLE CONSERVATION OF ANGULAR MOMENTUM (NOT ACTIVATED
  //    194	C                    HERE, COPIED FROM P_NUCJ)
  //    195	C        P-NU7=SAME AS P-NU6 + EVAPORATION                              P-N00330
  //    196	C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++P-N00340
  //    197	 
  //    198	      DIMENSION f(15),kind(300),ep(300),alpha(300),beta(300),gam(300)
  //    199	      DIMENSION bmass(300)
  //    200	      COMMON/hazard/ial,IY1,IY2,IY3,IY4,IY5,IY6,IY7,IY8,IY9,IY10,
  //    201	     s               IY11,IY12,IY13,IY14,IY15,IY16,IY17,IY18,IY19
  //    202	      COMMON/kind/kindf7
  //    203	      DIMENSION IND(20000),JND(20000)                                   P-N00350
  //    204	      DIMENSION INDIC(3000)                                             P-N00360
  //    205	      DIMENSION NPAR(625,15),NIMP(600,15),NNCO(15),NIMPP(600)           P-N00370
  //    206	      DIMENSION NENTR(10,8,20,10),NOUT1(15),NOUT2(15)                   P-N00380
  //    207	 
  //    208	c     DIMENSION TEM(15),NSR(40),NSP(40),NSR1(40),NSP1(40)               P-N00400
  //    209	      DIMENSION TEM(15),NSR1(40),NSP1(40)                               P-N00400
  //    210	      DIMENSION T(200),LINE(132),Q1(200),Q2(200),Q3(200),Q4(200),NC(300)P-N00410
  //    211	      DIMENSION Y1(200),Y2(200),Y3(200),YM(200),IPI(200)                P-N00420
  //    212	      DIMENSION NRNN(15),NRND(15),NRDD(15),NRDN(15),NRDP(15),NRPD(15),NCP-N00430
  //    213	     -DD(15),NPAUL1(15),NPAUL2(15)                                      P-N00440
  //    214	      DIMENSION NPDIR(600)                                              P-N00450
  //    215	      DIMENSION NEJ(6,15),NRES(6,15),NPIA(6,15),NCHPRO(15),NDEL(15)     P-N00460
  //    216	      DIMENSION EDEP1(15),EDEP2(15),EDEP3(15),EDEP4(15),NG4INT(15)        P-N00470
  //    217	     -,EPAR1(15),EPAR2(15),EPAR3(15),EPAR4(15),ENPI(15),E1(15),EZ3(15)  P-N00480
  //    218	      DIMENSION IHF1(50),IHF2(50),IHF3(50),IHF4(50),IHP(2,100),IHC(50), P-N00490
  //    219	     -IHE(2,100),IHF5(100),IHREM(100,100)
  //    220	     
  //    221	      DIMENSION JPARTICIP(300),eps_c(4),p3_c(4)
  //    222	     
  //    223	C Dialogue with INCL: function R(q/pf) for each nucleus
  //    224	      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPG4INTER,IMAT
  //    225	      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX      
  //    226	
  //    227	C RMS espace R, espace P, Fermi momentum and energy for light gauss nuc.      
  //    228	      COMMON/light_gaus_nuc/rms1t(9),pf1t(9),pfln(9),tfln(9),vnuc(9)
  //    229	      
  //    230	C Fermi 2 param from A=19 to 28, modified harm oscil A=6 to 18
  //    231	C (H. De Vries et al. At. Data and Nuc. Data Tab. 36 (1987) 495)
  //    232	      COMMON/light_nuc/R_light_nuc(30),a_light_nuc(30)      
  //    233	
  //    234	c common for study of avatars through an ntuple optionally produced
  //    235	      real*4 bavat,timeavat,energyavat
  //    236	      G4integer bloc_paul,bloc_cdpp,go_out,avm,del1avat,del2avat
  //    237	      parameter (avm=1000)
  //    238	      common/var_avat/kveux,bavat,nopartavat,ncolavat,
  //    239	     s r1_in(3),r1_first_avat(3),
  //    240	     s epsd(250),eps2(250),eps4(250),eps6(250),epsf(250),
  //    241	     s nb_avat,
  //    242	     s timeavat(avm),l1avat(avm),l2avat(avm),jpartl1(avm),jpartl2(avm),
  //    243	     s del1avat(avm),del2avat(avm),energyavat(avm),
  //    244	     s bloc_paul(avm),bloc_cdpp(avm),go_out(avm)
  //    245	
  //    246	      dimension npproj(300)
  //    247	      common/spl2/ xsp(100),ysp(100),ad(100),bd(100),cd(100),ndeut
  //    248	c deutons      
  //    249	c     dimension ltt(15)                                                 p-n00520
  //    250	      common/bl1/p1(300),p2(300),p3(300),eps(300),ind1(300),ind2(300),tap-n00530
  //    251	      common/bl2/crois(19900),k,ind,jnd                                 p-n00540
  //    252	      common/bl3/r1,r2,x1(300),x2(300),x3(300),ia1,ia2,rab2             p-n00550
  //    253	      common/bl4/tmax5                                                  p-n00560
  //    254	      common/bl5/tlg(300),nesc(300)                                     p-n00570
  //    255	      common/bl6/xx10,isa                                               p-n00580
  //    256	      common/bl8/rathr,ramass
  //    257	      common/bl9/hel(300),l1,l2
  //    258	      common/bl10/ri4,rs4,r2i,r2s,pf
  //    259	      common/paul/ct0,ct1,ct2,ct3,ct4,ct5,ct6,pr,pr2,xrr,xrr2,
  //    260	     s            cp0,cp1,cp2,cp3,cp4,cp5,cp6
  //    261	     
  //    262	      dimension ia1t(9),iz1t(9),fmpinct(9)      
  //    264	      data line/132*1h*/,hc,fmp,fmd,fmpi/197.328,938.2796,1232.,138.00/ p-n00590
  G4double hc = 197.328;
  G4double fmp = 938.2796;
  G4double fmpi = 138.00;
  //    265	      data /ia1t,iz1t,fmpinct/1,1,0,0,0,2,3,3,4,1,0,1,0,-1,1,1,2,2,
  //    266	     -938.2796,938.2796,138.0,138.0,138.0,1874.35,2806.8,2806.8,3727.
  G4int ia1t[9] = {1,1,0,0,0,2,3,3,4};
  G4int iz1t[9] = {1,0,1,0,-1,1,1,2,2};
  G4double fmpinct[9] = {-938.2796,938.2796,138.0,138.0,138.0,1874.35,2806.8,2806.8,3727.};
  //    266	     
  //    267	c      data rms1t,pf1t/0.,0.,0.,0.,0.,2.10,1.80,1.80,1.63,
  //    268	c     -0.,0.,0.,0.,0.,77.,110.,110.,153./
  //    269	
  //    270	
  //    271	c deutons
  //    272	      data nentr/16000*0/                                               p-n00620
  //    389	c                                                                       p-n01700
  //    390	ccc   reading of the data                                               p-n01710
  //    391	c                                                                       p-n01720
  //    392	c-------explanation of some physical quantities------------------------ p-n01730
  //    393	c           (basically input data, when run as a program)
  //    394	c                                                                       p-n01740
  //    395	c     ia1=mass number of the incident ion                               p-n01750
  //    396	c                                                                       p-n01760
  //    397	c     ia2=mass number of the target                                     p-n01770
  //    398	c     iz1=atomic number of the projectile                               p-n01780
  //    399	c                                                                       p-n01790
  //    400	c     iz2=atomic number of the target                                   p-n01800
  //    401	c     r01=radius parameter of the projectile                            p-n01810
  //    402	c           r01 should be put to 1.000                                  p-n01820
  //    403	c     r02=radius parameter of the target                                p-n01830
  //    404	c     adif1=diffuseness of the projectile                               p-n01840
  //    405	c           adif1 should be put to 1.0000                               p-n01850
  //    406	c     adif2=diffuseness of the target                                   p-n01860
  //    407	c                                                                       p-n01870
  //    408	c     tlab = incident energy  (in mev)                                  p-n01880
  //    409	c                                                                       p-n01890
  //    410	c     k1=1 reference frame=c.m. of the covering matter                  p-n01900
  //    411	c     k1=2 reference frame=c.m. frame of the incident ion and its       p-n01910
  //    412	c          G4intercept                                                    p-n01920
  //    413	c     k1=3 reference frame=c.m. frame of a n-n system with the same     p-n01930
  //    414	c          kinematics as the two ions                                   p-n01940
  //    415	c     k1=4 reference frame=c.m. frame for the total  system             p-n01950
  //    416	c     k1=5 reference frame=lab system                                   p-n01960
  //    417	c              k1 should be put to 5 in this version                    p-n01970
  //    418	c     k2=0 relativistic kinematics                                      p-n01980
  //    419	c     k2=1 non-relativistic kinematics(abandonned in this version)      p-n01990
  //    420	c     k3=0  deltas are produced                                         p-n02000
  //    421	c     k3=1  no delta production                                         p-n02010
  //    422	c     k4=0 the delta is given a vanishing lifetime                      p-n02020
  //    423	c     k4=1 the delta has a very large lifetime                          p-n02030
  //    424	c     k4=2 the delta has a exponentially random lifetime                p-n02040
  //    425	c     k5=0 no delta-nucleon,delta-delta G4interactions                    p-n02050
  //    426	c     k5=1 delta-nucleon=delta-delta=nucleon-nucleon elastic x-section  p-n02060
  //    427	c     k6=0 no angular momentum conservation                             p-n02070
  //    428	c     k6=1 angular momentum conservation                                p-n02080
  //    429	c                                                                       p-n02090
  //    430	c     b=impact parameter                                                p-n02100
  //    431	c                                                                       p-n02110
  //    432	c     rbl(pbl) is the radius in real(momentum) space of the volume      p-n02120
  //    433	c         on which the nucleons are counted to evaluate the pauli       p-n02130
  //    434	c         blocking factors                                              p-n02140
  //    435	c     recommended values are 2 fm and 200 mev/c respectively            p-n02150
  //    436	c                                                                       p-n02160
  //    437	c     nrun= number of  runs (irrelevant here)                           p-n02170
  //    438	c                                                                       p-n02180
  //    439	c     ntm=number of G4intermediate times at which the spatial and momentump-n02190
  //    440	c         distributions are stored (not relevant here)                  p-n02200
  //    441	c                                                                       p-n02210
  //    442	c     tem(it)=values of the G4intermediate times                          p-n02220
  //    443	c                                                                       p-n02230
  //    444	c     v0=depth of the nucleon potential                                 p-n02240
  //    445	c     v1=depth of the delta potential                                   p-n02250
  //    446	c         in this version v1 should be equal to v0                      p-n02260
  //    447	c         (attention! a  v0 > 0 corresponds to an attractive potential) p-n02270
  //    448	c                                                                       p-n02280
  //    449	c     ncase=number of cells necessary to depict the spatial distributionp-n02290
  //    450	c         in the x and z directions                                     p-n02300
  //    451	c                                                                       p-n02310
  //    452	c     xp,xg,zp,zg=limits of the box in the x and z directions           p-n02320
  //    453	c     in order to have the same scale in the two directions ,put xg-xp= p-n02330
  //    454	c         0.9*(zg-zp)                                                   p-n02340
  //    455	c     dy=dimension of the box in the y direction                        p-n02350
  //    456	c                                                                       p-n02360
  //    457	c     rap1,rap2,pp1,pp2=limits of the box in the rapidity-p(perpendicu- p-n02370
  //    458	c          -lar) space                                                  p-n02380
  //    459	c     the box is divided G4into 15 cells along p(perp) and 40 cells along p-n02390
  //    460	c         the rapidity                                                  p-n02400
  //    461	c                                                                       p-n02410
  //    462	c     en,sn=average energy carried away by a nucleon,separation energy  p-n02420
  //    463	c                                                                       p-n02430
  //    464	c       (ntm,ncase,xp,...,sn are not used here)
  //    465	c-----------------------------------------------------------------------p-n02440

  G4int jparticip[300];
  G4double beproj = 0.;
  bl3->ia2 = G4int(calincl->f[0]); // f(1)->f[0] and so on..., calincl added
  G4int iz2 = G4int(calincl->f[1]);
  G4double r02 = 1.12;
  kindstruct->kindf7 = G4int(calincl->f[6] + 0.1);

  bl3->ia1 = ia1t[G4int(kindstruct->kindf7)];
  G4int iz1 = iz1t[G4int(kindstruct->kindf7)];
  G4double fmpinc = fmpinct[G4int(kindstruct->kindf7)];
  G4double rms1 = light_gaus_nuc->rms1t[G4int(kindstruct->kindf7)];
  G4double pf1 = light_gaus_nuc->pf1t[G4int(kindstruct->kindf7)];
  G4double tlab = calincl->f[2];

  G4int k1 = 5;
  G4int k2 = 0;
  G4int k3 = 0;
  G4int k4 = 2;
  G4int k5 = 1;
  G4int k6 = 0;

  // material number:      
  saxw->imat = G4int(calincl->f[8] + 0.5); // f(9) -> f[8]
  // espace de phases test (r et p) pour pauli: 
  // valeur recommandee par j.c. v-test=0.589 h**3:
  G4double rbl = 2.0;
  G4double pbl=200.;

  paul->xrr = rbl;
  paul->xrr2 = (paul->xrr) * (paul->xrr);
  paul->pr=pbl;
  paul->pr2 = paul->pr*(paul->pr);

  G4double tem[10];
  tem[0] = 100000.0;  // tem(1) -> tem[0]
  // temfin (time at which the inc is stopped), tmax5 defined after chosing b

  G4double v0 = calincl->f[4]; // f(5)->f[4]
  G4double v1 = v0;
  bl8->rathr = 0.;
  bl8->ramass = 0.;

  // constants and derived data
  bl10->pf = 1.37*hc;
  G4double tf = sqrt(bl10->pf*(bl10->pf)+fmp*fmp)-fmp;
  G4double g0 = 115.0;
  G4double th = 0.;
  G4double pm2 = fmp*fmp;
  G4int ia = bl3->ia1 + bl3->ia2;
  G4int a2 = bl3->ia2;
  bl3->r2 = r02*pow(a2,0.33333333);

  // parametres moyens de densite de la cible (fermi 2 parametres)
  if (bl3->ia2 > 28) { //then
    ws->r0 = (2.745e-4*bl3->ia2+1.063)*pow(bl3->ia2,0.33333333);
    ws->adif = 1.63e-4*bl3->ia2 + 0.510;
    ws->rmaxws = ws->r0 + ws->xfoisa*(ws->adif);
  }
  else if(bl3->ia2 >= 0.19) { //then
    ws->r0 = light_nuc->r[bl3->ia2];
    ws->adif = light_nuc->a[bl3->ia2];
    ws->rmaxws = ws->r0 + ws->xfoisa*(ws->adif);
  }
  else if(bl3->ia2>=6) { //then
    ws->r0 = 1.581*(light_nuc->a[bl3->ia2]) * (2.0 + 5.0 * (light_nuc->r[bl3->ia2])/(2.0 + 3.0*(light_nuc->r[bl3->ia2])));
    ws->adif = light_nuc->a[bl3->ia2];
    ws->rmaxws = 5.5 + 0.3*(bl3->ia2 - 6.0)/12.0;
  }
  else if(bl3->ia2 >= 2) { // then
    if(bl3->ia2 == 2) { //then
      ws->r0 = light_gaus_nuc->rms1t[5]; // rms1t(6) -> rms1t[5]
      bl10->pf = light_gaus_nuc->pfln[5]; //pfln(6)->pfln[5]
      tf = light_gaus_nuc->tfln[5]; // tfln(6) -> tfln[5]
      v0 = light_gaus_nuc->vnuc[5]; // vnuc(6) -> vnuc(5)
    } //endif
    if(bl3->ia2 == 3 && iz2 == 1) { //then
      ws->r0 = light_gaus_nuc->rms1t[6]; //rms1t(7)->rms1t[6] and so on...
      bl10->pf = light_gaus_nuc->pfln[6];
      tf = light_gaus_nuc->tfln[6];
      v0 = light_gaus_nuc->vnuc[6];
    } //endif
    if(bl3->ia2 == 3 && iz2 == 2) { //then
      ws->r0 = light_gaus_nuc->rms1t[7]; //rms1t(8)->rms1t[7] and so on...
      bl10->pf = light_gaus_nuc->pf1t[7];
      tf = light_gaus_nuc->tfln[7];
      v0 = light_gaus_nuc->vnuc[7];
    } //endif
    if(bl3->ia2 == 4) { //then
      ws->r0 = light_gaus_nuc->rms1t[8]; // rms1t(9) -> rms1t[8] and so on...
      bl10->pf = light_gaus_nuc->pf1t[8];
      tf = light_gaus_nuc->tfln[8];
      v0 = light_gaus_nuc->vnuc[8];
    } //endif
    v1 = v0;
    ws->adif = 0.57735*(ws->r0);
    ws->rmaxws = ws->r0+2.5;
  } // end if
  if(ws->nosurf > 0) { //then
    ws->adif=0.0;
    ws->rmaxws=ws->r0;
  } //end if
  ws->drws = ws->rmaxws/29.0;

  // on impose le rayon du noyau:
  // ....voir coherence function wsax(r) et derivwsax(r)
  bl3->r2 = ws->r0;

  G4double tnr = tlab;
  G4double binc = sqrt(tnr*tnr + 2.0*tlab*fmpinc)/(tnr+fmpinc);
  G4double ginc=1.0/sqrt(1.0 - binc*binc);
  G4double pinc = fmpinc*binc*ginc;

  for(G4int i = 0; i < ia; i++) {
    jparticip[i]=0;
  }

  // skip initialisations
  G4double efrun = 0.0;
  G4int iavat = 0;

  // generation of the initial distribution in the rest frame of the iop-n03870
  // surface
  //    700	c      bmax=r2+2.2*adif
  //    701	c      r2i=(r2-2.2*adif)
  //    702	c      r2s=(r2+2.2*adif)
  //    703	c **********
  //    704	      

  ws->bmax = ws->rmaxws;		// maximum extension of the nucleus ( w.s.)

  //  fin surface (que faire aux pions ?)
  //    709	c      if (kindf7.gt.2) bmax=r2+2.2
  //    710	c      if (kindf7.gt.2) bmax=bmax	! a.b. (avec w.s., idem les nucleons)
  // deutons cv 22/01/2001
  if (kindstruct->kindf7 <= 2)  { //then
    ws->bmax = ws->bmax;     // comme alain
  }
  else {
    if (kindstruct->kindf7 < 6) { // then    
      ws->bmax = ws->bmax;   // comme alain
    }
    else {
      beproj = fmpinc - bl3->ia1*fmp;
      ws->bmax = ws->rmaxws + rms1;     // maximum extension of the nucleus ( w.s.)
    }
  }

  // deutons     
  G4double al;
  standardRandom(&al, &(hazard->ial));
  G4double b = sqrt(al)*(ws->bmax);
  G4double bred = b/bl3->r2;
  //G4double bimpact=b;
  bimpact = b;
  G4double tnor;

  if(ws->nosurf != -2) { // la suite, c'est la version temps avant 2001
    if(ws->nosurf <= 0) {
      tnor = 70.0/30.6349;         
    }
    // PK endif removed
    else {
      tnor=1.;
    } //endif               
    if(ws->nosurf == 0) {
      bred = 0.;
    }
    if (kindstruct->kindf7 <= 2) {
      if (tlab < 400.0) {
	cb0 = 6.86 - 0.0035 * tlab;
	eb0 = 0.32 - 0.00005 * tlab;
      }
      else {
	if (tlab < 1000.0) { //then
	  cb0 = 5.23 + 0.000575 * tlab;
	  eb0 = 0.32 - 0.00005 * tlab;
	}
	else {
	  cb0 = 5.73 + 0.00007 * tlab;
	  eb0 = 0.283 - 0.000013 * tlab;
	}
      }
      temfin = 1.25*cb0/amax1(1.0,0.854 + 0.438*bred)*pow(bl3->ia2,(eb0/amax1(1.0,0.941+0.177*bred)));
      temfin = temfin*tnor;
    }
    else {
      if (kindstruct->kindf7 < 6) {
	// here for pions:
	temfin = 30.0*pow((float(bl3->ia2)/208.0),0.25)*(1.0 - 0.2*bred)*(1.0 - tlab/1250.0);
	// correction for pions in the case nosurf=0 or -1 (a.b., c.v. 2/2002)
	temfin = temfin*tnor;
      }
      else {
	// deutons
	tlabu = tlab/bl3->ia1;
	if (tlabu <= 400) {
	  coeffb0 = -0.0035*tlabu + 6.86;
	  expob0 = -0.00005*tlabu + 0.32;
	}
	else {
	  if (tlabu <= 1000) { //then    
	    coeffb0 = 0.000575*tlabu + 5.23;
	    expob0 = -0.00005*tlabu + 0.32;
	  }
	  else {
	    coeffb0 = 0.00007*tlabu + 5.73;
	    expob0 = -0.000013*tlabu + 0.283;
	  }
	}
	if (bred <= 0.33333) {
	  xc = 1.0;
	  xe = 1.0;
	}
	else {
	  xc = 0.438*bred + 0.854;
	  xe = 0.177*bred + 0.941;
	}
	temfin = 1.25*(coeffb0/xc)*pow((bl3->ia2),(expob0/xe));
	// same renormalisation of time for p,n and composit particles.
	temfin = temfin*tnor;
      }
    }
  }
  else { //ici,nosurf=-2 c'est la fonction temps 2001 (avec surface).
    if(kindstruct->kindf7 >= 3 && kindstruct->kindf7 <= 5) {
      // here for pions (arbitrary multiplied by 2 for surface) (a.b., c.v.) 2/2002:
      // temfin=60.*(float(ia2)/208.)**0.25*(1.-0.2*bred)*(1.-tlab/1250.)
      // modified in april 2003 (more reasonable but not yet checked!)
      temfin = 25.5*pow(bl3->ia2,0.16);  // pb208->60fm/c
    }
    else {
      // here for other hadrons
      temfin = 29.8*pow(bl3->ia2,0.16);  // pb208->70fm/c
    }
  }

  // deutons
  // a way to change stopping time f[5] not used here
  factemp = calincl->f[5]; // f(6)->f[5]
  // attention !!! 30/04/2001 scaling time is now a multiplication factor
  temfin = temfin*factemp;

  exi = 0.0;
  nbquit = 0;
  iqe = 0;
  idecf = 0;
  bl4->tmax5 = temfin+0.1;
  npion = 0;
  efer = 0.0;

  // deutons
  if (bl3->ia1 > 1) {
    goto pnu7;
  }

  // deutons
  if (bl3->ia1 == 1) { //then
    bl5->nesc[0] = 0; // nesc(1)->nesc[0] and so on...
    bl1->ind2[0] = 2*iz1 - 1;
    bl1->ind1[0] = 0; 
    bl3->x1[0] = 0.0;
    bl3->x2[0] = 0.0;
    bl3->x3[0] = 0.0;
    bl1->p1[0] = 0.0;
    bl1->p2[0] = 0.0;
    bl1->p3[0] = 0.0;
    bl9->hel[0] = 0.0;
    jparticip[0] = 1;
  }
  else {
    npion = 1;
    ipi[0] = 8 - 2*(kindstruct->kindf7);
    y1[0] = 0.0;
    y2[0] = 0.0;
    y3[0] = 0.0;
    q1[0] = 0.0;
    q2[0] = 0.0;
    q3[0] = 0.0;
    q4[0] = fmpi;
  }

  // deutons
  goto pnu9;
  //    850	
 pnu7: 
  s1t1 = 0.0;
  s2t1 = 0.0;
  s3t1 = 0.0;
  sp1t1 = 0.0;
  sp2t1 = 0.0;
  sp3t1 = 0.0;
  for(G4int i = 0; i < (bl3->ia1-1); i++) {
    bl9->hel[i] = 0;
    bl5->nesc[i] = 0;
    bl1->ind2[i] = 1;
    if (i > iz1) {
      bl1->ind2[i] = -1;
    }
    bl1->ind1[i] = 0;
    // deutons
    jparticip[i] = 1;
    std::cout <<"Particle " << i << " is now participant." << std::endl;
    // deutons 
    gaussianRandom(&xga);
    bl3->x1[i] = xga*rms1*0.57735;
    s1t1 = s1t1 + bl3->x1[i];
    gaussianRandom(&xga);
    bl3->x2[i] = xga*rms1*0.57735;
    s2t1 = s2t1 + bl3->x2[i];
    gaussianRandom(&xga);
    bl3->x3[i] = xga*rms1*0.57735; 
    s3t1 = s3t1 + bl3->x3[i];

    if(kindstruct->kindf7 == 6) { //then
      // deuteron density from paris potential in q space:      
      standardRandom(&xq, &(hazard->igraine[9]));
      qdeut = splineab(xq) * 197.3289;
      standardRandom(&u, &(hazard->igraine[10]));
      cstet = u*2 - 1;
      sitet = sqrt(1.0 - pow(cstet,2));
      standardRandom(&v, &(hazard->igraine[11]));
      phi = 2.0*3.141592654*v;

      bl1->p1[i] = qdeut*sitet*cos(phi);
      bl1->p2[i] = qdeut*sitet*sin(phi);
      bl1->p3[i] = qdeut*cstet;
    }
    else {
      // density of composite as a gaussien in q space:
      gaussianRandom(&xga);
      bl1->p1[i] = xga*pf1*0.57735;
      gaussianRandom(&xga);
      bl1->p2[i] = xga*pf1*0.57735;
      gaussianRandom(&xga);
      bl1->p3[i] = xga*pf1*0.57735;
    }

    bl1->eps[i] = w(bl1->p1[i],bl1->p2[i],bl1->p3[i],fmp);

    sp1t1 = sp1t1 + bl1->p1[i];
    sp2t1 = sp2t1 + bl1->p2[i];
    sp3t1 = sp3t1 + bl1->p3[i];
  }

  bl9->hel[bl3->ia1] = 0;
  bl5->nesc[bl3->ia1] = 0;
  bl1->ind2[bl3->ia1] = -1;
  bl1->ind1[bl3->ia1] = 0;
  bl3->x1[bl3->ia1] = -s1t1;
  bl3->x2[bl3->ia1] = -s2t1;
  bl3->x3[bl3->ia1] = -s3t1;
  bl1->p1[bl3->ia1] = -sp1t1;
  bl1->p2[bl3->ia1] = -sp2t1;
  bl1->p3[bl3->ia1] = -sp3t1;
  bl1->eps[bl3->ia1] = w(bl1->p1[bl3->ia1],bl1->p2[bl3->ia1],bl1->p3[bl3->ia1],fmp);

  std::cout <<"Reached line 1746" << std::endl;
  
  // deutons
  jparticip[bl3->ia1] = 1;
    std::cout <<"Particle " << bl3->ia1 << " is now participant." << std::endl;
 pnu9: // continue
  // deutons
  // target preparation for 1 < a < 5 (with sum of momentum =0)
  std::cout <<"Placing nucleons into the nucleus..." << std::endl;
  if(bl3->ia2 >= 2 && bl3->ia2 <= 4) {
  pnu1633:  
    s1t1 = 0.0;
    s2t1 = 0.0;
    s3t1 = 0.0;
    sp1t1 = 0.0;
    sp2t1 = 0.0;
    sp3t1 = 0.0;
    efer = 0.0;
    for(G4int i = bl3->ia1+1; i < ia-1; i++) {
      bl1->ind2[i] = 1;
      bl5->nesc[i] = 0;
      if (i > (iz2+bl3->ia1)) {
	bl1->ind2[i] = -1;
      }
      for(G4int j = 0; j < 7; j++) {
	standardRandom(&t[j], &(hazard->ial));
      }

      t[1] = -1.0 + 2.0*t[1];      // t(2)->t[1]                                        
      t[2] = 6.283185*t[2];   // t(3)->t[2]                                     
      t[4] = -1.0 + 2.0*t[4];     // t(5)->t[4]                                  
      t[5] = 6.283185*t[5];   // t(6) -> t[5]                                      
      t1 = t[1];              // t(2)->t[1]                                           
      t2 = sqrt(1.0 - t1*t1);                                       
      t3 = cos(t[2]);  //t(3)->t[2]                                            
      t4 = sin(t[2]);   //t(3)->t[2]                                                                                      
      t5 = t[4];       // t(5)->t[4]                                                 
      t6 = sqrt(1.0 - t5*t5);                                         
      t7 = cos(t[5]);   //t(6) -> t[5]                                             
      t8 = sin(t[5]);   // t(6)->t[5]                                          
      if (ws->nosurf == 1) {
	x = bl3->r2*pow(t[0],0.33333333); // t(1)->t[0]                                      
	y = (bl10->pf)*pow(t[3],0.33333333); // t(4)->t[3]  
      }
      else {
	// surface..w.s.: impulsion (sphere dure), puis r(q)
	t33 = pow(t[6],0.33333333); // t(7)->t[6]
	y = (bl10->pf)*t33;

	rr = interpolateFunction(t33);
	x = rr*pow(t[3],0.33333333); // t(4)->t[3] 
	// fin surface on a redefini x et y
      }
      bl3->x1[i] = x*t2*t3;                                                
      bl3->x2[i] = x*t2*t4; 
      bl3->x3[i] = x*t1;                                                    
      s1t1 = s1t1 + bl3->x1[i];
      s2t1 = s2t1 + bl3->x2[i];
      s3t1 = s3t1+bl3->x3[i];
      bl1->p1[i] = y*t6*t7; 
      bl1->p2[i] = y*t6*t8;                                                  
      bl1->p3[i] = y*t5;                                                   
      sp1t1 = sp1t1+bl1->p1[i];
      sp2t1 = sp2t1+bl1->p2[i];
      sp3t1 = sp3t1+bl1->p3[i];
      bl1->ind1[i] = 0;                                                      
      bl1->eps[i] = w(bl1->p1[i],bl1->p2[i],bl1->p3[i],fmp);                            
      bl9->hel[i] = 0.0;
      efer = efer + bl1->eps[i] - fmp;                                        
    }

    bl9->hel[ia] = 0; 
    bl5->nesc[ia] = 0;
    bl1->ind2[ia] = -1;
    bl1->ind1[ia] = 0;
    bl3->x1[ia] = -s1t1;
    bl3->x2[ia] = -s2t1;
    bl3->x3[ia] = -s3t1;
    bl1->p1[ia] = -sp1t1;
    bl1->p2[ia] = -sp2t1;
    bl1->p3[ia] = -sp3t1;
    p_mod = sqrt(pow(bl1->p1[ia],2) + pow(bl1->p2[ia],2) + pow(bl1->p3[ia],2));
    if(p_mod > ((bl10->pf)+0.05)) {
      goto pnu1633;
    }
    bl1->eps[ia] = w(bl1->p1[ia],bl1->p2[ia],bl1->p3[ia],fmp);
    efer = efer + bl1->eps[ia]-fmp;                                        
    //    998	c for rho(r),rho(q) checking
    //    999	c      do i=bl3->ia1+1,ia
    //   1000	c      r_dist = sqrt(x1(i)*x1(i)+x2(i)*x2(i)+x3(i)*x3(i))
    //   1001	c      call hfill(3,r_dist,0.,1.)
    //   1002	c      p_mod = sqrt(p1(i)**2+p2(i)**2+p3(i)**2)
    //   1003	c      call hfill(4,p_mod,0.,1.)
    //   1004	c      enddo


  } //end if       !(bl3->ia2 >= 2 && bl3->ia2 <= 4)

  // target preparation for a > 4
  if(bl3->ia2 > 4) {
    x1_target = 0.0;
    x2_target = 0.0;
    x3_target = 0.0;
    for(G4int i = bl3->ia1; i < ia; i++) { //do 1 i=bl3->ia1+1,ia
      bl5->nesc[i] = 0;
      bl1->ind2[i] = 1;
      if (i > (iz2+bl3->ia1)) {
	bl1->ind2[i] = -1;
      }
      // surface ajout de t(7) surface.f avait do 6 j=2,7 ici 1,6 ?
      for(G4int j = 0; j < 7; j++) {
	standardRandom(&t[j], &(hazard->ial));
      }
      //    pnu6: // continue
      t[1] = -1.0 + 2.0*t[1]; //t(2)->t[1]
      t[2] = 6.283185*t[2];
      t[4] = -1.0 + 2.0*t[4]; // t(5)->t[4]
      t[5] = 6.283185*t[5]; //t(6)->t[5]
      t1 = t[1]; // t(2)->t[1]
      t2 = sqrt(1.0 - t1*t1);
      t3 = cos(t[2]);  //t(3)->t[2]
      t4 = sin(t[2]); //t(3)->t[2]
      t5 = t[4]; //t(5)->t[4]
      t6 = sqrt(1.-t5*t5);
      t7 = cos(t[5]);  //t(6)->t[5]
      t8 = sin(t[5]);   // t(6)->t[5]
      if (ws->nosurf == 1) {
	x = bl3->r2*pow(t[0],0.33333333); // t(1)->t[0]
	y = bl10->pf*pow(t[3],0.33333333); // t(4)->t3
      }
      else {
	// surface..w.s.: impulsion (sphere dure), puis r(q)
	t33 = pow(t[6],0.33333333); // t(7)->t[6]
	std::cout <<"t33 = " << t33 << std::endl;
	y=bl10->pf*t33;
	rr=interpolateFunction(t33);
	std::cout <<"rr = " << rr << std::endl;
	x=rr*pow(t[3],0.33333333); // t(4)->t[3]       
	// fin surface on a redefini x et y
      }
      bl3->x1[i] = x*t2*t3;
      bl3->x2[i] = x*t2*t4;
      bl3->x3[i] = x*t1;
      bl1->p1[i] = y*t6*t7;
      bl1->p2[i] = y*t6*t8;
      bl1->p3[i] = y*t5;
      x1_target = x1_target + bl3->x1[i];
      x2_target = x2_target + bl3->x2[i];
      x3_target = x3_target + bl3->x3[i];
      bl1->ind1[i] = 0;
      bl1->eps[i] = w(bl1->p1[i],bl1->p2[i],bl1->p3[i],fmp);
      bl9->hel[i] = 0.0;
      efer = efer + bl1->eps[i] - fmp;
    }

    x1_target = x1_target/bl3->ia2;
    x2_target = x2_target/bl3->ia2;
    x3_target = x3_target/bl3->ia2;

  }

  efrun = efrun + efer;

  // location of incident particle at point (b,z)
  r22 = bl3->r2*(bl3->r2);
  z = ws->bmax * (ws->bmax) - b*b;       // for the wood-saxon density...

  if (z < 0.0) {
    z=0.0;
  }

  z = sqrt(z);
  // random azimuthal direction of the impact parameter (sept 99)

  if (kindstruct->kindf7 <= 2) {
    standardRandom(&tbid, &(hazard->igraine[13]));
    tbid = tbid*6.283185;
    bl3->x1[0] = bl3->x1[0] + b*cos(tbid);  //x1(1)->x1[0]                                        
    bl3->x2[0] = bl3->x2[0] + b*sin(tbid); //x2(1)->x2[0]
    bl3->x3[0] = bl3->x3[0] - z;
    // pour le ntuple des avatars:
    if(varavat->kveux == 1) {
      varavat->r1_in[0] = bl3->x1[0]; //r1_in(1)->r1_in[0] and x1(1)->x1[0]
      varavat->r1_in[1] = bl3->x2[0]; //r1_in(2)->r1_in[1] and x1(2)->x1[1]
      varavat->r1_in[2] = bl3->x3[0]; //r1_in(3)->r1_in[2] and x1(3)->x1[2]
    } //endif
  }
  else {
    if (kindstruct->kindf7 > 6) { //then ! pour les pions on laisse
      //call standardRandom(tbid,iy14)
      standardRandom(&tbid, &(hazard->igraine[13]));
      tbid = tbid*6.283185;
      y1[0] = y1[0] + b*cos(tbid); //y1(1)->y1[0]
      y2[0] = y2[0] + b*sin(tbid); //y2(1)->y2[0]
      y3[0] = y3[0] - z;
    }
    else {
      // deutons
      //nmiss=0.;
      nmiss = 0;
      xlengm=1000.0;
      
      for(G4int i = 0; i < bl3->ia1; i++) {
	bl3->x3[i] = bl3->x3[i]/ginc;
	zai2 = ws->rmaxws*(ws->rmaxws) - pow((b+bl3->x1[i]),2) - pow(bl3->x2[i],2);
	if (zai2 > 0.0) {
	  goto pnu22;
	}
	ztu = -sqrt(zai2);
	// r22 remplace par rmaxws*rmaxws et r2 par rmaxws cv correct ?
	za_i = 2.0*(ws->rmaxws) + ztu;
	xleng = za_i - bl3->x3[i];
	if (xleng > xlengm) {
	  goto pnu21;
	}
	ilm = i;
	xlengm = xleng;
	ztouch = ztu;
	goto pnu21;
      pnu22:
	nmiss = nmiss + 1; 
      }
    pnu21:
      if (nmiss == bl3->ia1) { //then
	nopart = -1;
	//	return;
	std::cout <<"nmiss == bl3->ia1" << std::endl;
	goto pnureturn;
      }
      else {
	zshif = bl3->x3[ilm] - ztouch;
	standardRandom(&tbid, &(hazard->igraine[13]));
	tbid = tbid*6.283185;
	for(G4int i = 0; i < bl3->ia1; i++) {
	  xxx = bl3->x1[i] + b;
	  bl3->x1[i] = xxx*cos(tbid) - bl3->x2[i]*sin(tbid);
	  bl3->x2[i] = xxx*sin(tbid) + bl3->x2[i]*cos(tbid);
	  bl3->x3[i] = bl3->x3[i] - zshif;
	  //   1135	c           x1(i)=x1(i)+b
	  //   1136	c           x3(i)=x3(i)-zshif
	} //enddo
	if (fabs(pow(bl3->x1[ilm],2)+pow(bl3->x2[ilm],2)+pow(bl3->x3[ilm],2)-ws->rmaxws*(ws->rmaxws)) > 0.01) {
	  std::cout <<"wrong position" << std::endl;
	}
      } //endif
    } //endif         
  } //endif  

  std::cout <<"Reached line 1986..." << std::endl;

  // initial momentum for all type of incident particles:
  xl1 = b*pinc*sin(tbid);                                  
  xl2 = -b*pinc*cos(tbid);                                          
  xl3 = 0.0;                                           

  // transcription in the general frame of reference
  // (here,=lab frame)

  be = 0.0;
  ge = 1.0;
  b1 = (binc - be)/(1.0 - be*binc);
  b2 = -be;
  g1 = 1.0/sqrt(1.0 - b1*b1);
  g2 = 1.0;
  // deutons
  // here for nucleons
  if (kindstruct->kindf7 <= 2) {
    bl1->eps[0] = g1*fmp + v0; // eps(1)->eps[0]
    bl1->p3[0] = sqrt(pow(bl1->eps[0],2) - pow(fmp,2)); //p3(1)->p3(0), eps(1)->eps(0)
  }
  else {
    // here for pions
    if (kindstruct->kindf7 > 6) { //then
      q4[0] = g1*fmpi; // q4(1)->q4[0]
      q3[0] = b1*q4[0];
    }
    else {
      // here for composite projectiles:
      // the kinetic energy is below the threshold. put all
      // fermi momentum to 0... projectile nucleons not on shell!
      energie_in = tlab + fmpinc;
      if((energie_in) <= (bl3->ia1*fmp)) {
	for(G4int i = 0; i < bl3->ia1; i++) {
	  bl1->eps[i] = energie_in/bl3->ia1;
	  bl1->p1[i] = 0.0;
	  bl1->p2[i] = 0.0;
	  bl1->p3[i] = pinc/bl3->ia1;
	}
	goto pnu1871;
      }
      // here the composit is above threshold
      for(G4int i = 0; i < bl3->ia1; i++) { //do i=1,bl3->ia1	!save e,p in the composit rest frame
	eps_c[i] = bl1->eps[i];
	p3_c[i] = bl1->p3[i];
      } //enddo

      nbtest = bl3->ia1 - 1;
      if(kindstruct->kindf7 == 6) {
	nbtest = 2;
      }
      iflag = 0;
    pnu1870:
      sueps = 0.0;

      iflag = iflag + 1;

      for(G4int i = 0; i < bl3->ia1; i++) { //do i=1,bl3->ia1
	tte = eps_c[i];
	bl1->eps[i] = g1*(eps_c[i] + b1*p3_c[i]);
	bl1->p3[i] = g1*(b1*tte + p3_c[i]);
	sueps = sueps + bl1->eps[i];
      } //enddo

      cobe = (tlab + fmpinc)/sueps;

      // off shell problem for incident clusters (a.b. 2/2002)

      if(iflag == nbtest) { // too much..all momentum to 0
	for(G4int klm = 0; klm < bl3->ia1; klm++) { //do klm=1,bl3->ia1
	  eps_c[klm] = fmp;
	  bl1->p1[klm] = 0.0;
	  bl1->p2[klm] = 0.0;
	  p3_c[klm] = 0;
	}
	goto pnu1870;
      }
      for(G4int i = 0; i < bl3->ia1; i++) { //do i=1,bl3->ia1
	arg = pow((cobe*(bl1->eps[i])),2)-pm2;
	if (arg <= 0.) { //then	! put maximum momentum to 0. 
	  i_emax = 1; //	!find maximum
	  ener_max = bl1->eps[0]; //eps(1)->eps[0]
	  for(G4int i = 1; i < bl3->ia1; i++) { //do klm=2,bl3->ia1	
	    if(bl1->eps[klm] > ener_max) {
	      ener_max = bl1->eps[klm];
	      i_emax = klm;
	    }
	  }
	  eps_c[i_emax] = fmp;
	  bl1->p1[i_emax] = 0.0;
	  bl1->p2[i_emax] = 0.0;
	  p3_c[i_emax] = 0.0;

	  if(i_emax == bl3->ia1) { //    circular permut if the last one
	    epsv = eps_c[bl3->ia1]; //	 permutation circulaire
	    p1v = bl1->p1[bl3->ia1];
	    p2v = bl1->p2[bl3->ia1];
	    p3v = p3_c[bl3->ia1];
	    for(G4int k = bl3->ia1-1; k >= 0; k--) { //do k=bl3->ia1-1,1,-1
	      eps_c[k+1] = eps_c[k];
	      bl1->p1[k+1] = bl1->p1[k];
	      bl1->p2[k+1] = bl1->p2[k];
	      p3_c[k+1] = p3_c[k];
	    }
	    eps_c[0] = epsv;  //eps_c(1)->eps_c[0]
	    bl1->p1[0] = p1v;  //p1(1)->p1[0]
	    bl1->p2[0] = p2v;  //p2(1)->p2[0]
	    p3_c[0] = p3v; //p3_c(1)->p3_c(0)		fin permut.
	  }
	  sp1t1 = 0.0;   // re-compute the last one 
	  sp2t1 = 0.0;
	  sp3t1 = 0.0;
	  for(G4int j = 0; j < bl3->ia1-1; j++) { //do j=1,bl3->ia1-1
	    sp1t1 = sp1t1 + bl1->p1[j];
	    sp2t1 = sp2t1 + bl1->p2[j];
	    sp3t1 = sp3t1 + p3_c[j];
	  }
	  bl1->p1[bl3->ia1] = -sp1t1;
	  bl1->p2[bl3->ia1] = -sp2t1;
	  p3_c[bl3->ia1] = -sp3t1;
	  eps_c[bl3->ia1] = w(bl1->p1[bl3->ia1],bl1->p2[bl3->ia1],p3_c[bl3->ia1],fmp);	

	  goto pnu1870;  // ..and boost all of them.
	}
      }

      for(G4int i = 0; i < bl3->ia1; i++) { //do i=1,bl3->ia1
	arg = pow((cobe*(bl1->eps[i])),2) - pm2;
	comom = sqrt(arg/(pow(bl1->eps[i],2) - pm2));
	bl1->p1[i] = comom*(bl1->p1[i]);
	bl1->p2[i] = comom*(bl1->p2[i]);
	bl1->p3[i] = comom*(bl1->p3[i]);
	bl1->eps[i] = bl1->eps[i]*cobe;
	if (fabs(am(bl1->p1[i],bl1->p2[i],bl1->p3[i],bl1->eps[i])-fmp) > 0.01) {
	  std::cout <<"wrong correction " << i << std::endl;                  
	}
      }
      bl1->eps[ilm] = bl1->eps[ilm] + v0;  

    }
  }

 pnu1871:
  // evaluation of the times t(a,b)
  k = 0;
  kcol = 0;
  std::cout <<"Evaluating times..." << std::endl;
  if (kindstruct->kindf7 <= 2) {
    // modif s.vuillier tient compte propagation projectile,1e collision
    // imposee pour lui (c'est une maniere de faire!!)
    for(G4int i = 0; i < ia; i++) { //do 40 i=1,ia
      //		tref=ref(x1(i),x2(i),x3(i),p1(i),p2(i),p3(i),eps(i),r22)          p-n04740
      tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
      if (tref > bl4->tmax5) {
	std::cout <<"tref > bl4->tmax5" << std::endl;
	std::cout <<"tref = " << tref << " \t \t bl4->tmax5 = " << bl4->tmax5 << std::endl;
	goto pnu45;
      }
      k=k+1;
      bl2->crois[k]=tref;
      bl2->ind[k]=i;
      bl2->jnd[k]=-1;
    pnu45:
      i1=i-1;
      //	  if (i == 1) 
      if (i == 0) {
	std::cout <<"i == 0 near pnu45" << std::endl;
	//goto pnu40;
	continue;
      }
      //   1326	c ici on ne calcule que les G4interactions nn impliquant le projectile !!! (2/02)
      if (i1 > bl3->ia1) {
	i1=bl3->ia1;
      }
      //   1328	c**********************************************************************               
      for(G4int j = 0; j < i1; j++) { //do 41 j=1,i1
	//   1330	c    no collisions before the first collision of the incident particle
	//   1331	c      do 41 j=1,1
	//call time (i,j)
	time (i, j);
	if (bl1->ta < 0.) {
	  continue;
	}
	if(bl1->ta > bl4->tmax5) {
	  continue;
	}
	eij=am(bl1->p1[i]+bl1->p1[j],bl1->p2[i]+bl1->p2[j],bl1->p3[i]+bl1->p3[j],bl1->eps[i]+bl1->eps[j]);
	if (eij < 1925.) {
	  continue;
	}
	isos=bl1->ind2[i]+bl1->ind2[j];
	if (31.*(bl3->rab2) > totalCrossSection(eij,0,isos)) {
	  continue;
	}
	k=k+1;
	//if (j == 1) 
	if (j == 0) {
	  kcol=kcol+1;
	}
	bl2->crois[k]=bl1->ta;
	bl2->ind[k]=i;
	bl2->jnd[k]=j;
      }
      continue;
    }
  }
  else {
    // deutons
    if (kindstruct->kindf7 < 6) { //then
      // here for incoming pions:
      std::cout <<"Incoming pions..." << std::endl;
      for(G4int i = bl3->ia1+1; i < ia; i++) { //do i=ia1+1,ia
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
	if (tref < bl4->tmax5) {
	  k=k+1;
	  bl2->crois[k]=tref;
	  bl2->ind[k]=i;
	  bl2->jnd[k]=-1;
	}
      }
      new2(y1[0], y2[0], y3[0], q1[0], q2[0], q3[0], q4[0], 1, 0);

      //   modif a.b. 21/06/2002: should check at least one valid collision
      //   1361	c    with incoming pion.
      //   1362	c      kcol=1
      if(k != 0) {
	kcol = 1;
      }
    }
    else {
      for(G4int i = 0; i < bl3->ia1; i++) { //do 38 i=1,ia1
	bl5->nesc[i] = 1;
	if (i != ilm) { 
	  goto pnu36;
	}
	//   1371	      tref=ref(x1(i),x2(i),x3(i),p1(i),p2(i),p3(i),eps(i),r22)
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
	bl5->nesc[i] = 0;
	npproj[i] = 0;
	goto pnu37;
      pnu36:
	t1 = bl3->x1[i]*(bl1->p1[i])+bl3->x2[i]*(bl1->p2[i])+bl3->x3[i]*(bl1->p3[i]);                      
	t2 = bl1->p1[i]*(bl1->p1[i])+bl1->p2[i]*(bl1->p2[i])+bl1->p3[i]*(bl1->p3[i]);               
	t3 = t1/t2;
	t4 = bl3->x1[i]*(bl3->x1[i])+bl3->x2[i]*(bl3->x2[i])+bl3->x3[i]*(bl3->x3[i]); 
	//   1379	c incoming nucleons enter potential at maximum radius (modif. 13/06/01)
	t5 = t3*t3 + ((ws->rmaxws)*(ws->rmaxws) - t4)/t2;
	if (t5 < 0.) {
	  continue;
	}
	tref = (-1*t3 - sqrt(t5))*(bl1->eps[i]);  
	if (tref > bl4->tmax5) {
	  continue;
	}
	npproj[i] = 1;
      pnu37:
	k = k + 1; 
	bl2->crois[k] = tref; 
	bl2->ind[k] = i; 
	bl2->jnd[k] = -1;
      }
      kcol = 1;
      
      for(G4int i = bl3->ia1+1; i < ia; i++) { //do  39 i=ia1+1,ia
	npproj[i] = 0;
	//   1394	      tref=ref(x1(i),x2(i),x3(i),p1(i),p2(i),p3(i),eps(i),r22)
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
	if (tref < bl4->tmax5) { //then
	  k=k+1;
	  bl2->crois[k]=tref;
	  bl2->ind[k]=i;
	  bl2->jnd[k]=-1;
	} //endif

	time (i, ilm);
	if (bl1->ta < 0.) {
	  continue;
	}
	if (bl1->ta > bl4->tmax5) {
	  continue;
	}
	eij=am(bl1->p1[i]+bl1->p1[ilm],bl1->p2[i]+bl1->p2[ilm],bl1->p3[i]+bl1->p3[ilm],bl1->eps[i]+bl1->eps[ilm]);
	if (eij < 1925.) {
	  continue;
	}
	isos=bl1->ind2[i]+bl1->ind2[ilm];                                             
	if (31.*(bl3->rab2) > totalCrossSection(eij,0,isos)) {
	  continue;
	}								
	k=k+1;                                                             
	kcol=kcol+1; 
	bl2->crois[k]=bl1->ta;                                                      
	bl2->ind[k]=i; 
	bl2->jnd[k]=ilm;
      }
    }
  }
  // deutons
  if (kcol != 0) {
    goto pnu48;
  }
  nopart = -1;
  //   1419	c pour eviter renvoi des resultats du run precedent cv 7/7/98
  iarem = bl3->ia2;
  izrem = iz2;
  esrem = 0.0;
  erecrem = 0.0;
  // fin ajout cv
  //   1425	      return
  std::cout <<"kcol == 0" << std::endl;
  goto pnureturn;

  // initialization  at the beginning of the run
 pnu48:
  timi = 0.0;
  tim = 0.0;
  ncol = 0;

  // compteur des collisions a deux corps (call collis acceptes par pauli)
  ncol_2c = 0;

  // npion=0
  mrnn = 0;
  mrnd = 0;
  mrdd = 0;
  mrdn = 0;
  mrdp = 0;
  mrpd = 0;
  mcdd = 0;
  mpaul2 = 0;
  mpaul1 = 0;

  // approx. (not considering escaping protons of incident clusters) 11/03 a.b.
  itch = iz1 + iz2 - 1;
  for(G4int i = 0; i < ia; i++)  { //do 47 i=1,ia
    bl5->tlg[i] = 0.0;
    nc[i] = 0;
  }
  itt = 1;

  // tableau des energies a l'initialisation
  if(varavat->kveux == 1) {
    for(G4int i = bl3->ia1+1; i < ia; i++) {
      varavat->epsd[i] = bl1->eps[i];
    }
    iflag20 = 0;
    iflag40 = 0;
    iflag60 = 0;
  }

  // search for the smallest positive t(a,b)
  // pour tests, =0 G4interdit les  reflexions avant un avatar du projectile,
  //             =1 comme avant (reflexions autorisees). (a.b. 3/2002)
  irst_avatar = 1;

 pnu449:
  std::cout <<"pnu 449..." << std::endl;
  next = 1;
  indic[next] = 1;
 pnu44:
  if(next == 0) {
    goto pnu449;
  }
  
  idep = indic[next] + 1;
  tau = bl2->crois[idep-1];
  if(idep > k) {
    std::cout <<"going to pnu448" << std::endl;
    goto pnu448;
  }
  for(G4int i = idep; i < k; i++) { //do 42 i=idep,k
    if (bl2->crois[i] > tau) {
      continue;
    }
    tau = bl2->crois[i]; 
    next = next + 1;
    indic[next] = i;
  }
 pnu448:
  imin = indic[next];
  l1 = bl2->ind[imin];
  l2 = bl2->jnd[imin];

  // test le 20/3/2003: tue sinon le dernier avatar?
  if (k == 0) {
    std::cout <<"going to pnu230" << std::endl;
    goto pnu230;
  }
  k = k - 1;
  next = next - 1;

  // correction s.vuillier 25/1/96 decalage temps correct
  if (imin > k) {
    std::cout <<"going to pnu46..." << std::endl;
    goto pnu46;
  }

  for(G4int i = imin; i < k; i++) { //do 43 i=imin,k
    bl2->crois[i] = bl2->crois[i+1];
    bl2->ind[i] = bl2->ind[i+1];
    bl2->jnd[i] = bl2->jnd[i+1];
  }
 pnu46:
  tim = timi + tau;

  // tableau des energies a t=20,40,60 fm/c
  if(varavat->kveux == 1) {
    if(iflag20 == 0 && tim >= 20.0) {
      iflag20 = 1;
      for(G4int i = 0; i < ia; i++) {
	if(bl5->nesc[i] == 0) {
	  if(jparticip[i] == 1) {
	    varavat->eps2[i] = bl1->eps[i];
	  }
	  else {
	    varavat->eps2[i] = 0.0;
	  }
	}
	else {
	  varavat->eps2[i] = 0.;
	}
      }
    }
    
    if(iflag40 == 0 && tim >= 40.0) {
      iflag40 = 1;
      for(G4int i = 0; i < ia; i++) {
	if(bl5->nesc[i] == 0) {
	  if(jparticip[i] == 1) {
	    varavat->eps4[i] = bl1->eps[i];
	  }
	  else {
	    varavat->eps4[i] = 0.0;
	  }
	}
	else {
	  varavat->eps4[i] = 0.0;
	}
      }
    }
    
    if(iflag60 == 0 && tim >= 60.) {
      iflag60=1;
      for(G4int i = 0; i < ia; i++) {
	if(bl5->nesc[i] == 0) {
	  if(jparticip[i] == 1) {
	    varavat->eps6[i] = bl1->eps[i];
	  }
	  else {
	    varavat->eps6[i] = 0.0;
	  }
	}
	else {
	  varavat->eps6[i] = 0.0;
	}
      }
    }
  }

  // modif: pas de reflexions avant au moins un avatar du (des) nucleon incident
  // celui-ci ne peut etre qu'une collision nn (ou pin)

  if((irst_avatar == 0) && (l2 == -1)) {
    goto pnu44;
  }

  irst_avatar = irst_avatar+1;

  if (tim < temfin) {
    std::cout <<"tim = " << tim << " temfin = " << temfin << std::endl;
    std::cout <<"tim < temfin, going to pnu49" << std::endl;
    goto pnu49;
  }

  goto pnu255;
 pnu49:
  if (k == 0) {
    std::cout <<"k == 0. going to pnu255" << std::endl;
    goto pnu255;
  }
  // l1 va a la surface du noyau:
  if (l2 == -1) {
    std::cout <<"l2 == -1. going to pnu220" << std::endl;
    goto pnu220;
  }

  if((k4-1) <= 0) {
    std::cout <<"k4-1 <= 0. going to pnu255" << std::endl;
    goto pnu803;
  }
  if((k4-1) > 0) {
    std::cout <<"k4-1 > 0. going to pnu255" << std::endl;
    goto pnu830;
  }

  // l1 est un delta: 
 pnu830: 
  if(l2 == 0) {
    goto pnu220;
  }
  // interaction pi(l1-ia)-nucleon(l2)
  if(l1 > ia) {
    goto pnu801;
  }
 pnu803:
  std::cout <<"pnu803:" << std::endl;
  // pas de collision entre 2 non participants:
  if(jparticip[l1] == 0 && jparticip[l2] == 0) {
    goto pnu44;
  }

  // parameters for the next colliding pair
  t[9] = bl1->eps[l1] + bl1->eps[l2]; //t(10)->t[9] 
  t0 = 1.0/t[9]; // t(10)->t[9]
  b1 = (bl1->p1[l1] + bl1->p1[l2])*t0;
  b2 = (bl1->p2[l1] + bl1->p2[l2])*t0;
  b3 = (bl1->p3[l1] + bl1->p3[l2])*t0;
  s = (1.0 - b1*b1 - b2*b2 - b3*b3)*t[9]*t[9]; //t(10)->t[9]
  sq = sqrt(s);

  if(sq < 1925.5) {
    goto pnu44;
  }

  bl1->ta = tau/bl1->eps[l1];
  x1l1 = bl3->x1[l1] + bl1->p1[l1]*(bl1->ta);
  x2l1 = bl3->x2[l1] + bl1->p2[l1]*(bl1->ta);
  x3l1 = bl3->x3[l1] + bl1->p3[l1]*(bl1->ta);
  bl1->ta = tau/bl1->eps[l2];
  x1l2 = bl3->x1[l2] + bl1->p1[l2]*(bl1->ta);
  x2l2 = bl3->x2[l2] + bl1->p2[l2]*(bl1->ta);
  x3l2 = bl3->x3[l2] + bl1->p3[l2]*(bl1->ta);

  // test on the minimum distance of approach
  t[10] = x1l1 - x1l2; //t(11)->t[10]
  t[11] = x2l1 - x2l2; //t(12)->t[11]
  t[12] = x3l1 - x3l2; //t(13)->t[12]
  t[13] = t[10]*t[10] + t[11]*t[11] + t[12]*t[12]; //t(N)->t[N-1]
  t[14] = b1*t[10] + b2*t[11] + b3*t[12]; //t(N)->t[N-1]
  t[15] = b1*b1 + b2*b2 + b3*b3; //t(16)->t[15]
  bb2 = t[13] + t[14]*t[14]/(1.0 - t[15]); //t(N)->t[N-1]

  // Replaced goto structure:
  // if (k3 == 1) go to 260
  // if (k4 == 0) go to 260
  mg=bl1->ind1[l1]+bl1->ind1[l2];
  if((k3 != 1) && (k4 != 0) && (mg == 1)) {
    isos=bl1->ind2[l1]+bl1->ind2[l2];
    // if (mg != 1) go to 260
    ldel=l2;
    if(mg-bl1->ind1[l1] == 0) {
      ldel=l1;
    }
    bl6->xx10 = sqrt(pow(bl1->eps[ldel],2) - pow(bl1->p1[ldel],2) - pow(bl1->p2[ldel],2) - pow(bl1->p3[ldel],2));
    bl6->isa = bl1->ind2[ldel];
    bmax2 = totalCrossSection(sq,mg,isos)/31.415926;
    if (k5 == 0 && mg != 0) {
      bmax2 = bmax2 - lowEnergy(sq,mg,isos)/31.415926;
    }
    // go to 261
  }
  else {
    bmax2 = totalCrossSection(sq,mg,isos)/31.41592;
  }

  if (bb2 < bmax2) {
    goto pnu220;
  }
  if (k == 0) {
    goto pnu230;
  }
  goto pnu44;
  // loop while((bb2 >= bmax2) && (k != 0)) (PK)
  // evaluation of the positions at time = tim
 pnu220:
  std::cout <<"Starting pnu220..." << std::endl;
  
  timi = tim;
  if(varavat->kveux == 1) { //then
    iavat = iavat + 1;
    varavat->timeavat[iavat] = tim;
    varavat->l1avat[iavat] = l1;
    varavat->l2avat[iavat] = l2;
    varavat->energyavat[iavat] = sq;

    if(l1 <= ia) {
      varavat->jpartl1[iavat] = jparticip[l1];
    }
    else {
      varavat->jpartl1[iavat] = 0;
    }

    if(l2 > 0) {
      varavat->jpartl2[iavat] = jparticip[l2];
    }
    else {
      varavat->jpartl2[iavat] = 0;
    }
  }

  // gel des nucleons non participants sur le premier avatar (nn)=(l1,1)      
  if (irst_avatar == 1) {
    for(G4int i = 0; i < l1; i = i + l1-0) {
      bl1->ta = tau/bl1->eps[i];                                                
      bl3->x1[i] = bl3->x1[i] + bl1->p1[i]*(bl1->ta);                                     
      bl3->x2[i] = bl3->x2[i] + bl1->p2[i]*(bl1->ta);                                    
      bl3->x3[i] = bl3->x3[i] + bl1->p3[i]*(bl1->ta);
    }
    for(G4int i = 0; i < k; i++) {
      bl2->crois[i] = bl2->crois[i] + tau;
    }
  }
  else {
    for(G4int i = 0; i < ia; i++) {
      bl1->ta = tau/bl1->eps[i];
      bl3->x1[i] = bl3->x1[i] + bl1->p1[i]*(bl1->ta);
      bl3->x2[i] = bl3->x2[i] + bl1->p2[i]*(bl1->ta);
      bl3->x3[i] = bl3->x3[i] + bl1->p3[i]*(bl1->ta);
    }
  }

  if(npion != 0) {
    for(G4int i = 0; i < npion; i++) {
      bl1->ta=tau/q4[i];
      y1[i] = y1[i] + q1[i]*(bl1->ta);
      y2[i] = y2[i] + q2[i]*(bl1->ta);
      y3[i] = y3[i] + q3[i]*(bl1->ta);
    }
  }

  if(l2 == 0) {
    std::cout <<"l2 == 0, going to pnu805..." << std::endl;
    goto pnu805;
  }
  // Candidate: if(l2!=0)...

  // reflexions sur le potentiel, sortie eventuelle de la particule:
  if (l2 == -1) {
    std::cout <<"l2 == -1, going to pnu600..." << std::endl;
    goto pnu600;
  }

  if(l1 > ia) {
    std::cout <<"l1 > ia, going to pnu831..." << std::endl;
    goto pnu831;
  }

  // collision of particles l1 and l2
  ich1 = bl1->ind1[l1];
  ich2 = bl1->ind1[l2];
  ich3 = bl1->ind2[l1];
  ich4 = bl1->ind2[l2];
  aml1 = sqrt(pow(bl1->eps[l1],2) - pow(bl1->p1[l1],2) - pow(bl1->p2[l1],2) - pow(bl1->p3[l1],2));
  aml2 = sqrt(pow(bl1->eps[l2],2) - pow(bl1->p1[l2],2) - pow(bl1->p2[l2],2) - pow(bl1->p3[l2],2));
  gl1 = bl1->eps[l1]/aml1;
  gl2 = bl1->eps[l2]/aml2;
  // l-conservation
  if (k6 == 1) {
    t[30] = (aml1*(bl3->x1[l1]) + aml2*(bl3->x1[l2]))/(aml1 + aml2); //t(31)->t[30]
    t[31] = (aml1*(bl3->x2[l1]) + aml2*(bl3->x2[l2]))/(aml1 + aml2); //t(32)->t[31]
    t[32] = (aml1*(bl3->x3[l1]) + aml2*(bl3->x3[l2]))/(aml1 + aml2); //t(33)->t[32]
    tt31 = bl3->x1[l1] - bl3->x1[l2];
    tt32 = bl3->x2[l1] - bl3->x2[l2];
    tt33 = bl3->x3[l1] - bl3->x3[l2];
    t[33] = (aml2*(bl1->p1[l1]) - aml1*(bl1->p1[l2]))/(aml1 + aml2); //t(34)->t[33]
    t[34] = (aml2*(bl1->p2[l1]) - aml1*(bl1->p2[l2]))/(aml1 + aml2); //t(35)->t[34]
    t[35] = (aml2*(bl1->p3[l1]) - aml1*(bl1->p3[l2]))/(aml1 + aml2); //t(36)->t[35]
    tt34 = bl1->p1[l1] + bl1->p1[l2];
    tt35 = bl1->p2[l1] + bl1->p2[l2];
    tt36 = bl1->p3[l1] + bl1->p3[l2];
  }

  // l-conservation
  t[20] = bl1->p1[l1]; //t(21)->t[20]
  t[21] = bl1->p2[l1]; //t(22)->t[21]
  t[22] = bl1->p3[l1]; //t(23)->t[22]
  t[23] = bl1->eps[l1]; //t(24)->t[23]
  t[24] = bl1->p1[l2]; //t(25)->t[24]
  t[25] = bl1->p2[l2]; //t(26)->t[25]
  t[26] = bl1->p3[l2]; //t(27)->t[26]
  t[27] = bl1->eps[l2]; //t(28)->t[27]

  // info delta ou nucleon:
  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[l1];
    varavat->del2avat[iavat] = bl1->ind1[l2];                                               
  }
  G4double minus_b1 = -1.0*b1;
  G4double minus_b2 = -1.0*b2;
  G4double minus_b3 = -1.0*b3;
  loren(&(bl1->p1[l1]), &(bl1->p2[l1]), &(bl1->p3[l1]), &minus_b1, &minus_b2, &minus_b3, &(bl1->eps[l1]));
  loren(&(bl1->p1[l2]), &(bl1->p2[l2]), &(bl1->p3[l2]), &minus_b1, &minus_b2, &minus_b3, &(bl1->eps[l2]));

  std::cout <<"Calling collis" << std::endl;
  collis(&(bl1->p1[l1]), &(bl1->p2[l1]), &(bl1->p3[l1]),
	 &(bl1->eps[l1]), &(bl1->p1[l2]), &(bl1->p2[l2]), &(bl1->p3[l2]), &(bl1->eps[l2]),
	 &(t[11]), &(t[12]), &(t[13]), &(t[14]), &np, &ip, &k2, &k3, &k4, &k5,
	 &(bl1->ind1[l1]), &(bl1->ind1[l2]), &(bl1->ind2[l1]), &(bl1->ind2[l2]));
  std::cout <<"Collis call complete. " << std::endl;
  
  loren(&(bl1->p1[l1]), &(bl1->p2[l1]), &(bl1->p3[l1]), &b1, &b2, &b3, &(bl1->eps[l1]));
  loren(&(bl1->p1[l2]), &(bl1->p2[l2]), &(bl1->p3[l2]), &b1, &b2, &b3, &(bl1->eps[l2]));

  if (bl1->ind1[l1] == 1) {
    goto pnu243;
  }
  xbl1 = pauliBlocking(l1, rbl, pbl);
  standardRandom(&rndm,&(hazard->igraine[10])); 
  if (rndm > (1.0 - xbl1)) {
    goto pnu248;
  }
 pnu243:
  if (bl1->ind1[l2] == 1) {
    goto pnu241;
  }
  xbl2 = pauliBlocking(l2, rbl, pbl);
  standardRandom(&rndm,&(hazard->igraine[10]));
  if (rndm > (1.0 - xbl2)) {
    goto pnu248;
  }
  goto pnu241;
 pnu248:
  mpaul1 = mpaul1 + 1;
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 1;
  }
  // restitution de l1 et l2 si rejet de la col. par pauli:
  bl1->p1[l1] = t[20]; //t(21)->t[20]
  bl1->p2[l1] = t[21]; //t(22)->t[21]
  bl1->p3[l1] = t[23]; //t(23)->t[22]
  bl1->eps[l1] = t[23]; //t(24)->t[23]
  bl1->p1[l2] = t[24]; //t(25)->t[24]
  bl1->p2[l2] = t[25]; //t(26)->t[25]
  bl1->p3[l2] = t[26]; //t(27)->t[26]
  bl1->eps[l2] = t[27]; //t(28)->t[27]
  bl1->ind1[l1] = ich1;
  bl1->ind1[l2] = ich2;
  bl1->ind2[l1] = ich3;
  bl1->ind2[l2] = ich4;
  if (k == 0) {
    goto pnu230;
  }
  for(G4int i = 0; i < k; i++) {
    bl2->crois[i] = bl2->crois[i] - tau;
  }

  // pour le temps de calcul (a.b. 02/2002)
  goto pnu44;

 pnu241:

  // la premiere collision a deux corps ne peut pas baisser l'energie
  // du nucleon de recul (bloque par pauli dans un noyau cible froid).
  // (ici, toujours l2 < l1)
  ncol_2c = ncol_2c + 1;
  if(ncol_2c == 1) {
    for (G4int icomp = 0; icomp < bl3->ia1; icomp++) {
      // test on the first collision modified 4/07/2001 for direct and exchange.
      if(icomp == l1 || icomp == l2) {
	xavant = min(t[23],t[27]); //t(24)->t[23], t(28)->t[27]
	xapres = min(bl1->eps[l1],bl1->eps[l2]);
	if(xapres <= xavant) {
	  if(varavat->kveux == 1) {
	    varavat->bloc_cdpp[iavat] = 1;
	  }
	  goto pnu248;
	}
      }
    }

    // pour le ntuple des avatars:
    if(varavat->kveux == 1) { //then
      varavat->r1_first_avat[0] = bl3->x1[0]; //(N)->[N-1]
      varavat->r1_first_avat[1] = bl3->x2[0]; //(N)->[N-1]
      varavat->r1_first_avat[2] = bl3->x3[0]; //(N)->[N-1]
    } //endif
  }
  else {
    // les collisions suivantes ne penvent conduire a un noyau de a nucleons
    // sous l'energie de fermi et dans une config. d'energie inferieure a
    // efer-(ia2-nbalttf)*tf).
    egs = 0.0;
    nbalttf = 0;
    for(G4int i = 0; i < ia; i++) {
      if(bl5->nesc[i] == 0) {
	if(sqrt(pow(bl1->p1[i],2)+pow(bl1->p2[i],2)+pow(bl1->p3[i],2)) < bl10->pf) { //then
	  nbalttf = nbalttf + 1;
	  egs = egs + bl1->eps[i] - fmp;
	}
      }
    }

    if(egs < (efer-(bl3->ia2-nbalttf)*tf)) {
      if(varavat->kveux == 1) {
	varavat->bloc_cdpp[iavat] = 1;
      }
      goto pnu248;
    }
  }

  if(varavat->kveux == 1) {
    varavat->bloc_cdpp[iavat] = 0;
    varavat->bloc_paul[iavat] = 0;
  }

  jparticip[l1] = 1;
  std::cout <<"Particle " << l1 << " is now participant." << std::endl;
  jparticip[l2] = 1;
  std::cout <<"Particle " << l2 << " is now participant." << std::endl;

  if (ws->nosurf <= 0) {
    // surface
    pppp = sqrt(pow(bl1->p1[l1],2) + pow(bl1->p2[l1],2) + pow(bl1->p3[l1],2));
    rrrr = sqrt(pow(bl3->x1[l1],2) + pow(bl3->x2[l1],2) + pow(bl3->x3[l1],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[l1] = bl3->x1[l1]*rcorr/rrrr;
	bl3->x2[l1] = bl3->x2[l1]*rcorr/rrrr;
	bl3->x3[l1] = bl3->x3[l1]*rcorr/rrrr;
      }
    }
    pppp = sqrt(pow(bl1->p1[l2],2) + pow(bl1->p2[l2],2) + pow(bl1->p3[l2],2));
    rrrr = sqrt(pow(bl3->x1[l2],2) + pow(bl3->x2[l2],2) + pow(bl3->x3[l2],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[l2] = bl3->x1[l2]*rcorr/rrrr;
	bl3->x2[l2] = bl3->x2[l2]*rcorr/rrrr;
	bl3->x3[l2] = bl3->x3[l2]*rcorr/rrrr;
      }
    }
  }
  
  if (np == 0) {
    goto pnu240;
  }

  npion = npion + 1;
  loren(&t[11], &t[12], &t[13], &b1, &b2, &b3, &t[14]); //t(N)->t[N-1]
  q1[npion] = t[11]; //t(12)->t[11]
  q2[npion] = t[12]; //t(13)->t[12]
  q3[npion] = t[13]; //t(14)->t[13]
  q4[npion] = t[14]; //t(15)->t[14]
 pnu240:
  ncol = ncol + 1;
  if (l2 != 1) {
    goto pnu870;
  }

  // critere pour la leading particle: avant impulsion longitudinale max l=1
  // change en fevrier 2002: leading part. = energie totale max (l=1)
  if (bl1->p3[l2] > bl1->p3[l1]) { 
    goto pnu870;
  }

  // attention, il faut mieux faire et selectionner la plus grande energie
  // des particules participantes (jparticip()=1) et dans le noyau (nesc()=0)!

  xr1 = bl1->p1[l1];
  xr2 = bl1->p2[l1];
  xr3 = bl1->p3[l1];
  xr4 = bl1->eps[l1];
  xr5 = bl3->x1[l1];
  xr6 = bl3->x2[l1];
  xr7 = bl3->x3[l1];
  xr8 = gl1;
  ixr1 = bl1->ind1[l1];
  ixr2 = bl1->ind2[l1];
  ixr3 = ich1;
  bl1->p1[l1] = bl1->p1[l2];
  bl1->p2[l1] = bl1->p2[l2];
  bl1->p3[l1] = bl1->p3[l2];
  bl1->eps[l1] = bl1->eps[l2];
  bl3->x1[l1] = bl3->x1[l2];
  bl3->x2[l1] = bl3->x2[l2];
  bl3->x3[l1] = bl3->x3[l2];
  gl1 = gl2;
  bl1->ind1[l1] = bl1->ind1[l2];
  bl1->ind2[l1] = bl1->ind2[l2];
  ich1 = ich2;
  bl1->p1[l2] = xr1;
  bl1->p2[l2] = xr2;
  bl1->p3[l2] = xr3;
  bl1->eps[l2] = xr4;
  bl3->x1[l2] = xr5;
  bl3->x2[l2] = xr6;
  bl3->x3[l2] = xr7;
  gl2 = xr8;
  bl1->ind1[l2] = ixr1;
  bl1->ind2[l2] = ixr2;
  ich2 = ixr3;

  if(ich1 + ich2 - bl1->ind1[l1] - bl1->ind1[l2] != 0 || (ich1 + ich2) != 1) {
    goto pnu870;
  }
  if (k == 0) {
    goto pnu870;
  }

  for(G4int i = 0; i < k; i++) {
    if((bl2->ind[i] != 1) || (bl2->jnd[i] != 0)) {
      if((bl2->ind[i] != l1) || (bl2->jnd[i] != 0)) {
	continue;
      }
      bl2->ind[i] = 1;
      break;
    }

    bl2->ind[i] = l1;
    break;
  }
  pnu870:
  bl5->tlg[l1] = th*(bl1->eps[l1])/sqrt(pow(bl1->eps[l1],2)-pow(bl1->p1[l1],2)-pow(bl1->p2[l1],2)-pow(bl1->p3[l1],2));
  bl5->tlg[l2] = th*(bl1->eps[l2])/sqrt(pow(bl1->eps[l2],2)-pow(bl1->p1[l2],2)-pow(bl1->p2[l2],2)-pow(bl1->p3[l2],2));
  nc[l1] = nc[l1] + 1;
  nc[l2] = nc[l2] + 1;
  led = 0;

  if((ich1+ich2-bl1->ind1[l1]-bl1->ind1[l2]) < 0) {
    mrnd = mrnd + 1;
  }
  if((ich1+ich2-bl1->ind1[l1]-bl1->ind1[l2]) == 0) {
    if((ich1+ich2-1) < 0) {
      mrnn = mrnn + 1;
    }
    if((ich1+ich2-1) == 0) {
      mrdd = mrdd + 1;
      led = 1;
    }
    if((ich1+ich2-1) == 0) {
      mcdd = mcdd + 1;
      led = 1;
    }
  }
  if((ich1+ich2-bl1->ind1[l1]-bl1->ind1[l2]) > 0) {
    mrdn = mrdn + 1;
  }


  // reevaluation of the times t(a,b) for (a or b)=(l1 or l2)

  // l-conservation 
  if (k6 == 1) {
    aml1 = am(bl1->p1[l1],bl1->p2[l1],bl1->p3[l1],bl1->eps[l1]);
    aml2 = am(bl1->p1[l2],bl1->p2[l2],bl1->p3[l2],bl1->eps[l2]);
    t[36] = (aml2*(bl1->p1[l1]) - aml1*(bl1->p1[l2]))/(aml1+aml2); //t(37)->t[36]
    t[37] = (aml2*(bl1->p2[l1]) - aml1*(bl1->p2[l2]))/(aml1+aml2); //t(38)->t[37]
    t[38] = (aml2*(bl1->p3[l1]) - aml1*(bl1->p3[l2]))/(aml1+aml2); //t(39)->t[38]
    t[39] = sqrt(t[33]*t[33] + t[34]*t[34] + t[35]*t[35]); //t(N)->t[N-1]
    t[40] = sqrt(t[36]*t[36] + t[37]*t[37] + t[38]*t[38]); //t(N)->t[N-1]
    rhopi = tt31*t[33] + tt32*t[34] + tt33*t[35]; //t(N)->t[N-1]
    t[42] = tt31 - rhopi*t[33]/pow(t[39],2); //t(N)->t[N-1]
    t[43] = tt32 - rhopi*t[34]/pow(t[39],2); //t(N)->t[N-1]
    t[44] = tt33 - rhopi*t[35]/pow(t[39],2); //t(N)->t[N-1]
    t[45] = sqrt(t[42]*t[42] + t[43]*t[43] + t[44]*t[44]); //t(N)->t[N-1]
    t[42] = t[42]/t[45]; //t(N)->t[N-1]
    t[43] = t[43]/t[45]; //t(N)->t[N-1]
    t[45] = t[45]/t[46];
    cif = (t[33]*t[36] + t[34]*t[37] + t[35]*t[38])/t[39]/t[40]; //t(N)->t[N-1]

    // trouble with forward scattering 22/3/95
    if(fabs(cif) > 1.) {
      cif = sign(1.0,cif);
    }
    sif = sqrt(1.0 - cif*cif);
    t[36] = (t[33]*cif/t[39] + t[42]*sif)*t[40]; //t(N)->t[N-1]
    t[37] = (t[34]*cif/t[39] + t[43]*sif)*t[40]; //t(N)->t[N-1]
    t[38] = (t[35]*cif/t[39] + t[44]*sif)*t[40]; //t(N)->t[N-1]
    tri = sqrt(tt31*tt31 + tt32*tt32 + tt33*tt33);
    cchi = rhopi/tri/t[39]; //t(40)->t[39]
    schi = sqrt(1.0 - cchi*cchi);
    c1 = cif*cchi - sif*schi;
    c2 = sif*cchi + cif*schi;
    tt31 = (c1*t[33]/t[39] + c2*t[42])*tri*t[39]/t[40]; //t(N)->t[N-1]
    tt32 = (c1*t[34]/t[39] + c2*t[43])*tri*t[39]/t[40]; //t(N)->t[N-1]
    tt33 = (c1*t[35]/t[39] + c2*t[44])*tri*t[39]/t[40]; //t(N)->t[N-1]
    bl3->x1[l1] = t[30] + aml2*tt31/(aml1 + aml2); //t(31)->t[30]
    bl3->x2[l1] = t[31] + aml2*tt32/(aml1 + aml2); //t(32)->t[30]
    bl3->x3[l1] = t[32] + aml2*tt33/(aml1 + aml2); //t(33)->t[32]
    bl3->x1[l2] = t[30] - aml1*tt31/(aml1 + aml2); //t(31)->t[30]
    bl3->x2[l2] = t[31] - aml1*tt32/(aml1 + aml2); //t(32)->t[31]
    bl3->x3[l2] = t[32] - aml1*tt33/(aml1 + aml2); //t(33)->t[32]
    bl1->p1[l1] = aml1*tt34/(aml1 + aml2) + t[36]; //t(37)->t[36]
    bl1->p2[l1] = aml1*tt35/(aml1 + aml2) + t[37]; //t(38)->t[37]
    bl1->p3[l1] = aml1*tt36/(aml1 + aml2) + t[38]; //t(39)->t[38]
    bl1->eps[l1] = w(bl1->p1[l1],bl1->p2[l1],bl1->p3[l1],aml1);
    bl1->p1[l2] = aml2*tt34/(aml1 + aml2) - t[36]; //t(37)->t[36]
    bl1->p2[l2] = aml2*tt35/(aml1 + aml2) - t[37]; //t(38)->t[37]
    bl1->p3[l2] = aml2*tt36/(aml1 + aml2) - t[38]; //t(39)->t[38]
    bl1->eps[l2] = w(bl1->p1[l2],bl1->p2[l2],bl1->p3[l2],aml2);
  }
  // l-conservation

  if(k != 0) {
    kd = 0;
    ccr = tau;
    for(G4int i = 0; i < k; i++) {
      i20 = i - kd;
      if (k4 != 2 || led != 1) {
	if((bl2->ind[i] == l1) || (bl2->ind[i] == l2) || (bl2->jnd[i] == l1) || (bl2->jnd[i] == l2)) {
	  kd = kd + 1;
	  continue;
	}
      }
      else {
	if(bl2->jnd[i] == 0) {
	  if (bl2->ind[i] == l1 && bl1->ind1[l1] == 1) {
	    bl2->crois[i]=(bl2->crois[i]-ccr)*(bl1->eps[l1])/sqrt(pow(bl1->eps[l1],2)-pow(bl1->p1[l1],2)-pow(bl1->p2[l1],2)-pow(bl1->p3[l1],2))/gl1+ccr;
	  }
	  if (bl2->ind[i] == l2 && bl1->ind1[l2] == 1) {
	    bl2->crois[i]=(bl2->crois[i]-ccr)*(bl1->eps[l2])/sqrt(pow(bl1->eps[l2],2)-pow(bl1->p1[l2],2)-pow(bl1->p2[l2],2)-pow(bl1->p3[l2],2))/gl2+ccr;
	  }
	}
      }
    
      bl2->crois[i20]=bl2->crois[i]-ccr;
      bl2->ind[i20]=bl2->ind[i];
      bl2->jnd[i20]=bl2->jnd[i];
      continue;
    }
  
    k = k - kd;
  }

  newt(l1,l2);
  std::cout <<"Line 3318..." << std::endl;
  tref=ref(bl3->x1[l1], bl3->x2[l1], bl3->x3[l1], bl1->p1[l1], bl1->p2[l1], bl1->p3[l1], bl1->eps[l1],r22);

  if(tref <= bl4->tmax5) {
    k=k+1;
    bl2->crois[k]=tref;
    bl2->ind[k]=l1;
    bl2->jnd[k]=-1;
  }

  tref=ref(bl3->x1[l2], bl3->x2[l2], bl3->x3[l2], bl1->p1[l2], bl1->p2[l2], bl1->p3[l2], bl1->eps[l2],r22);

  if(tref <= bl4->tmax5) {
    k=k+1;
    bl2->crois[k]=tref;
    bl2->ind[k]=l2;
    bl2->jnd[k]=-1;
  }

  if (k4 == 2) {
    goto pnu848;
  }

  if(k < 0) {
    goto pnu230;
  }
  if(k == 0) {
    goto pnu230;
  }
  if(k > 0) {
    goto pnu449;
  }

 pnu848:
  if (npion == 0) {
    goto pnu844;
  }
  if (bl1->ind1[l1] == 1) {
    goto pnu843;
  }
  for(G4int k20 = 0; k20 < npion; k20++) {
    new3((y1[k20]), (y2[k20]), (y3[k20]), (q1[k20]), (q2[k20]), (q3[k20]), (q4[k20]), k20, l1);
  }

 pnu843:
  if(bl1->ind1[l2] != 1) {
    for(G4int k20 = 0; k20 < npion; k20++) {
      new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], -k20, l2);
    }
  }
 pnu844:

  if(bl1->ind1[l1]+bl1->ind1[l2] <= ich1+ich2) {
    goto pnu849;
  }
  if(bl1->ind1[l1]-ich1 != 1) {
    goto pnu820;
  }
  lnew = l1;
  goto pnu821;
 pnu820:
  if(bl1->ind1[l2]-ich2 != 1) {
    goto pnu849;
  }
  lnew = l2;

 pnu821:
  standardRandom(&rndm,&(hazard->igraine[16]));
  // largeur variable du delta (phase space factor G4introduced 4/2001)
  amlnew = sqrt(pow(bl1->eps[lnew],2)-pow(bl1->p1[lnew],2)-pow(bl1->p2[lnew],2)-pow(bl1->p3[lnew],2)); 
  geff = bl1->eps[lnew]/amlnew;
  qqq = sqrt((pow(amlnew,2) - pow((fmp+fmpi),2))*(pow(amlnew,2) - pow((fmp-fmpi),2)))/(2.0*amlnew);
  psf = pow(qqq,3)/(pow(qqq,3) + 5832000.0);
  tdel = -hc/(g0*psf)*log(rndm)*geff;                                  

  if(tdel <= bl4->tmax5) {
    k = k + 1;
    bl2->crois[k] = tdel;
    bl2->ind[k] = lnew;
    bl2->jnd[k] = 0;
  }

 pnu849:
  if (k == 0) { 
    goto pnu230;
  }
  goto pnu449;

  // decay of the delta particle                                       p-n09780
 pnu805:
  npion = npion + 1;
  ichd = bl1->ind2[l1];
  t[30] = bl1->p1[l1];   //t(31)->t[30]
  t[31] = bl1->p2[l1]; //t(32)->t[31]
  t[32] = bl1->p3[l1]; //t(33)->t[32]
  t[33] = bl1->eps[l1]; //t(34)->t[33]
  var_ab = pow(bl1->eps[l1],2) - pow(bl1->p1[l1],2) - pow(bl1->p2[l1],2) - pow(bl1->p3[l1],2);
  ym[npion] = 0.0;

  if(var_ab > 0.0) {
    ym[npion] = sqrt(var_ab);
  }

  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[l1];
    varavat->energyavat[iavat] = ym[npion];
  }

  std::cout <<"calling decay2" << std::endl;
  decay2(&(bl1->p1[l1]), &(bl1->p2[l1]), &(bl1->p3[l1]), &(bl1->eps[l1]), &q1[npion], &q2[npion], &q3[npion],
	 &q4[npion], &ym[npion], &fmp, &fmpi, &(bl9->hel[l1]));
  std::cout <<"end of decay2 call" << std::endl;
  // decay
  if (bl1->ind2[l1]*(bl1->ind2[l1]) == 9) {
    std::cout <<"going to pnu806" << std::endl;
    goto pnu806;
  }
  
  standardRandom(&rndm, &(hazard->ial));
  if (rndm < 0.333333333) {
    std::cout <<"rndm < (1/3) going to pnu837." << std::endl;
    std::cout <<"npion (used as an index of array ipi (ipi[npion])) = " << npion << std::endl;
    std::cout <<"l1 = " << l1 << " l2 = " << l2 << std::endl;
    goto pnu837;
  }
  
  ipi[npion]=0;
  goto pnu809;

 pnu837:
  ipi[npion]=bl1->ind2[l1]*2;
  bl1->ind2[l1]=-1*(bl1->ind2[l1]);
  goto pnu809;
 pnu806:
  bl1->ind2[l1]=bl1->ind2[l1]/3;
  ipi[npion]=2*(bl1->ind2[l1]);
 pnu809: // continue
  bl1->ind1[l1]=0;
  bl5->tlg[l1]=0.;

  // escape ?
  if (bl5->nesc[l1] > 0) {
    goto pnu850;
  }

  iteste = 0;
  xpb = pauliBlocking(l1, rbl, pbl);
  standardRandom(&rndm,&(hazard->igraine[10])); 

  // pauli blocking?
  if (rndm <= xpb) {
    goto pnu1848;
  }
  // le decay ne peut conduire a un noyau de a nucleons
  // sous l'energie de fermi et dans une config. d'energie inferieure a
  // c  efer-(ia2-nbalttf)*tf).
  egs = 0.0;
  nbalttf = 0;
  iteste = 1;
  for(G4int i = 0; i < ia; i++) {
    if(bl5->nesc[i] == 0) {
      if(sqrt(pow(bl1->p1[i],2)+pow(bl1->p2[i],2)+pow(bl1->p3[i],2)) < bl10->pf) {
	nbalttf = nbalttf + 1;
	egs = egs + bl1->eps[i] - fmp;
      }
    }
  }
  if(egs >= (efer-(bl3->ia2-nbalttf)*tf)) {
    goto pnu850;
  }

  // attention, logique negative!!! liberer le goto si on veut supprimer la
  // sequence precedente (cdpp sur delta-> pi n)

  if(varavat->kveux == 1) {
    varavat->bloc_cdpp[iavat] = 1;
  }

  // it is blocked!      
 pnu1848:
  mpaul2 = mpaul2 + 1;
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 1;
  }

  // largeur variable du delta (phase space factor G4introduced 4/2001)
  // (180.**3 = 5832000.)
  qqq = sqrt((pow(ym[npion],2) - pow((fmp+fmpi),2))*(pow(ym[npion],2) - pow((fmp-fmpi),2)))/(2.*ym[npion]);
  psf = pow(qqq,3)/(pow(qqq,3)+5832000.0);
  tdel = hc*t[33]/(g0*psf*ym[npion]); //t(34)->t[33]                                 

  if (iteste == 0) {
    tdel = tdel*xpb/(1.000001-xpb);
  }

  if(tdel <= bl4->tmax5) {
    k = k + 1;
    bl2->crois[k] = tdel;
    bl2->ind[k] = l1;
    bl2->jnd[k] = 0;
  }

  bl1->p1[l1] = t[30]; //t(31)->t[30]
  bl1->p2[l1] = t[31]; //t(32)->t[31]
  bl1->p3[l1] = t[32]; //t(33)->t[32]
  bl1->eps[l1] = t[33]; //t(34)->t[33]
  bl1->ind1[l1] = 1;
  bl1->ind2[l1] = ichd;
  npion = npion - 1;

  if (k == 0) {
    goto pnu230;
  }

  for(G4int i = 0; i < k; i++) {
    bl2->crois[i] = bl2->crois[i] - tau;
  }
  goto pnu449;

  // valid decay of the delta
 pnu850:
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 0;
    varavat->bloc_cdpp[iavat] = 0;
  }

  if (ws->nosurf <= 0) {
    // surface
    pppp = sqrt(pow(bl1->p1[l1],2) + pow(bl1->p2[l1],2) + pow(bl1->p3[l1],2));
    rrrr = sqrt(pow(bl3->x1[l1],2) + pow(bl3->x2[l1],2) + pow(bl3->x3[l1],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[l1] = bl3->x1[l1]*rcorr/rrrr;
	bl3->x2[l1] = bl3->x2[l1]*rcorr/rrrr;
	bl3->x3[l1] = bl3->x3[l1]*rcorr/rrrr;
      }
    }
  }

  ncol = ncol + 1;
  mrdp = mrdp + 1;
  y1[npion] = bl3->x1[l1];
  y2[npion] = bl3->x2[l1];
  y3[npion] = bl3->x3[l1];

  if (k == 0) {
    goto pnu4047;
  }
  
  kd = 0;
  ccr = tau;
  for(G4int i = 0; i < k; i++) {
    i20 = i - kd;

    if((bl2->ind[i] == l1) || (bl2->jnd[i] == l1)) {
      kd = kd + 1;
    }
    else {
      bl2->crois[i20] = bl2->crois[i] - ccr;
      bl2->ind[i20] = bl2->ind[i];
      bl2->jnd[i20] = bl2->jnd[i];
    }
  }
  k = k - kd;
 pnu4047:
  if (bl5->nesc[l1] != 0) {
    goto pnu845;
  }

  new1(l1);
  k = k + 1;
  bl2->crois[k] = ref(bl3->x1[l1], bl3->x2[l1], bl3->x3[l1], bl1->p1[l1], bl1->p2[l1], bl1->p3[l1], bl1->eps[l1],r22);
  bl2->ind[k] = l1;
  bl2->jnd[k] = -1;

  if(npion > 1) {
    n20 = npion - 1;
    for(G4int k20 = 0; k20 < n20; k20++) {
      new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], k20, l1);
    }
  }

 pnu845: 
  new2(y1[npion], y2[npion], y3[npion], q1[npion], q2[npion], q3[npion], q4[npion], npion, l1);
  if(k == 0) {
    goto pnu230;
  }
  goto pnu449;

  // pion-nucleon collision
 pnu801:
  lp = l1 - ia;
  dis1 = bl3->x1[l2]-y1[lp] + (bl1->p1[l2]/bl1->eps[l2] - q1[lp]/q4[lp])*tau;
  dis2 = bl3->x2[l2]-y2[lp] + (bl1->p2[l2]/bl1->eps[l2] - q2[lp]/q4[lp])*tau;
  dis3 = bl3->x3[l2]-y3[lp] + (bl1->p3[l2]/bl1->eps[l2] - q3[lp]/q4[lp])*tau;
  dist = dis1*dis1 + dis2*dis2 + dis3*dis3;
  t[9] = bl1->eps[l2] + q4[lp]; //t(10)->t[9]
  t0 = 1.0/t[9]; //t(10)->t[9]
  b1 = (bl1->p1[l2] + q1[lp])*t0;
  b2 = (bl1->p2[l2] + q2[lp])*t0;
  b3 = (bl1->p3[l2] + q3[lp])*t0;
  s = (1.0 - b1*b1 - b2*b2 - b3*b3)*t[9]*t[9]; //t(10)->t[9]
  sq = sqrt(s);
  cg = 4+bl1->ind2[l2]*ipi[lp];

  if(sq > 3000.0) {
    goto pnu832;
  }
  if(31.41592*dist > pionNucleonCrossSection(sq)*cg/6.) {
    goto pnu832;
  }
  goto pnu220;

 pnu832:
  if (k == 0) {
    goto pnu230;
  }
  goto pnu44;
 pnu831:
  standardRandom(&rndm, &(hazard->igraine[18]));
  geff = t[9]/sq; //t(10)->t[9]
  gg = g0;
  if (sq > 1500.0) {
    gg=200.0;
  }

  // largeur variable du delta (phase space factor G4introduced 4/2001)
  // (180.**3 = 5832000.)
  qqq = sqrt((pow(sq,2) - pow((fmp+fmpi),2))*(pow(sq,2) - pow((fmp-fmpi),2)))/(2.0*sq);
  psf = pow(qqq,3)/(pow(qqq,3)+5832000.);
  tdel = -hc/(gg*psf)*log(rndm)*geff;         

  bl1->ind1[l2] = 1;
  bl1->ind2[l2] = bl1->ind2[l2] + ipi[lp];
  nc[l2] = nc[l2] + 1;
  bl1->eps[l2] = t[9]; //t(10)->t[9]
  bl1->p1[l2] = bl1->p1[l2] + q1[lp];
  bl1->p2[l2] = bl1->p2[l2] + q2[lp];
  bl1->p3[l2] = bl1->p3[l2] + q3[lp]; 

  // ce nucleon (ici delta) devient un participant:
  jparticip[l2] = 1;
  std::cout <<"Particle " << l2 << " is now participant." << std::endl;
  
  if (ws->nosurf <= 0) {
    // surface
    pppp = sqrt(pow(bl1->p1[l2],2) + pow(bl1->p2[l2],2) + pow(bl1->p3[l2],2));
    rrrr = sqrt(pow(bl3->x1[l2],2) + pow(bl3->x2[l2],2) + pow(bl3->x3[l2],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[l2] = bl3->x1[l2]*rcorr/rrrr;
	bl3->x2[l2] = bl3->x2[l2]*rcorr/rrrr;
	bl3->x3[l2] = bl3->x3[l2]*rcorr/rrrr;
      }
    }
    // fin surface
  }

  // difference with standard cascade :
  // the delta is located at the nucleon site to avoid problems
  // with the reflexion on the wall
  if(lp != npion) {
    lp1 = lp+1;
    for(G4int i10 = lp1; i10 < npion; lp1++) {
      ipi[i10-1] = ipi[i10];
      ym[i10-1] = ym[i10];
      q1[i10-1] = q1[i10];
      q2[i10-1] = q2[i10];
      q3[i10-1] = q3[i10];
      q4[i10-1] = q4[i10];
      y1[i10-1] = y1[i10];
      y2[i10-1] = y2[i10];
      y3[i10-1] = y3[i10];
    }
  }
  npion = npion-1;
  ncol = ncol+1;
  mrpd = mrpd+1;

  if(k != 0) {
    kd = 0;
    ccr = tau;
    for(G4int i = 0; i < k; i++) {
      i20 = i - kd;

      if((bl2->ind[i] == l1) || (bl2->ind[i] == l2) || (bl2->jnd[i] == l1) || (bl2->jnd[i] == l2)) {
	kd = kd + 1;
      }
      else {
	bl2->crois[i20] = bl2->crois[i] - ccr;
	bl2->ind[i20] = bl2->ind[i];
	bl2->jnd[i20] = bl2->jnd[i];
      }
    }

    k = k - kd;
    for(G4int i10 = 0; i10 < k; i10++) {
      if(bl2->ind[i10] <= l1) {
	continue;
      }
      bl2->ind[i10] = bl2->ind[i10] - 1;
    }
  }
  
  new1(l2);

  if(tdel <= bl4->tmax5) {
    k = k + 1;
    bl2->crois[k] = tdel;
    bl2->ind[k] = l2;
    bl2->jnd[k] = 0;
  }

  k = k + 1;
  bl2->crois[k] = ref(bl3->x1[l2], bl3->x2[l2], bl3->x3[l2], bl1->p1[l2], bl1->p2[l2], bl1->p3[l2], bl1->eps[l2],r22);
  bl2->ind[k] = l2;
  bl2->jnd[k] = -1;
  goto pnu449;

  // reflection on or transmission through the potential wall
 pnu600:
  // deutons pas bien compris ici cv ?
  if (npproj[l1] == 0) {
    std::cout <<"npproj[l1] == 0, going to pnu608" << std::endl;
    goto pnu608;
  }
  
  if (bl1->ind1[l1] != 0) {
    std::cout << "wrong reentering particle" << std::endl;
  }
  
  if (bl3->x1[l1]*(bl1->p1[l1])+bl3->x2[l1]*(bl1->p2[l1])+bl3->x3[l1]*(bl1->p3[l1]) > 0.) {
    std::cout <<"wrong reentering particle" << std::endl;
  }
  
  var_ab = pow(bl1->p1[l1],2) + pow(bl1->p2[l1],2) + pow(bl1->p3[l1],2);
  gpsg = 0.0;
  if (var_ab > 0.0) {
    gpsg = sqrt((pow(bl1->eps[l1]+v0,2)-pm2)/var_ab); 
  }
  
  bl1->p1[l1] = gpsg*(bl1->p1[l1]);                                                
  bl1->p2[l1] = gpsg*(bl1->p2[l1]);                                               
  bl1->p3[l1] = gpsg*(bl1->p3[l1]);                                               
  bl1->eps[l1] = bl1->eps[l1] + v0;
  npproj[l1] = 0;
  bl5->nesc[l1] = 0;

  // reevaluation of the times tab after entrance of 2nd,..nucleon 
  // of the projectile (goto 602 instead of 607 modif. 13/06/01)
  std::cout <<"Going to 602" << std::endl;
  goto pnu602;

  // deutons
  // pour un non participant la transmission est impossible:
 pnu608:
  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[l1];
    varavat->energyavat[iavat] = bl1->eps[l1] - fmp;
  }
  if(jparticip[l1] == 0) {
    std::cout <<"jparticip[l1] == 0. going to pnu601." << std::endl;
    goto pnu601;
  }
  if(varavat->kveux == 1) {
    varavat->go_out[iavat]=1;
  }
  if (bl1->ind1[l1] == 0) {
    std::cout <<"bl1->ind1[l1] == 0. going to pnu605." << std::endl;
    goto pnu605;
  }
  fm = am(bl1->p1[l1],bl1->p2[l1],bl1->p3[l1],bl1->eps[l1]);
  pot = v1;
  goto pnu606;

 pnu605:
  fm = fmp;
  pot = v0;

 pnu606:
  tp = transmissionProb(bl1->eps[l1]-fm,bl1->ind2[l1],itch,bl3->r2,v0);
  if(varavat->kveux == 1) {
    varavat->energyavat[iavat] = bl1->eps[l1] - fm;
  }
  standardRandom(&rndm,&(hazard->igraine[10]));

  if (rndm > tp) {
    std::cout <<"rndm > tp, going to pnu601" << std::endl;
    goto pnu601;
  }

  // ici la particule l1 s'échappe du noyau:
  bl5->nesc[l1] = 1;
  nbquit = nbquit + 1;
  itch = itch - (1 + bl1->ind2[l1])/2;
  var_ab = pow(bl1->p1[l1],2) + pow(bl1->p2[l1],2) + pow(bl1->p3[l1],2);
  gpsg = 0.0;
  if(var_ab > 0.0) {
    gpsg = sqrt((pow(bl1->eps[l1]-pot,2) - fm*fm)/(var_ab));
  }
  bl1->p1[l1] = gpsg*(bl1->p1[l1]);
  bl1->p2[l1] = gpsg*(bl1->p2[l1]);
  bl1->p3[l1] = gpsg*(bl1->p3[l1]);
  bl1->eps[l1] = bl1->eps[l1] - pot;

  // comptage des particules hors du noyau (7/6/2002):
  // (remnant minimum=1 nucleon)
  if(nbquit >= (ia-1)) {
    goto pnu255;
  }
  goto pnu602;

  // here no transmission possible
 pnu601:
  std::cout <<"pnu601:" << std::endl;

  pspr=bl3->x1[l1]*(bl1->p1[l1])+bl3->x2[l1]*(bl1->p2[l1])+bl3->x3[l1]*(bl1->p3[l1]);
  if(varavat->kveux == 1) {
    varavat->go_out[iavat]=0;
  }

  // surface: modif a.b. pour tenir compte du rayon variable du noyau.
  // (x2cour remplace r22 le rayon**2 fixe du noyau)
  x2cour = pow(bl3->x1[l1],2) + pow(bl3->x2[l1],2) + pow(bl3->x3[l1],2);
  bl1->p1[l1] = bl1->p1[l1] - 2.0*(bl3->x1[l1])*pspr/x2cour;
  bl1->p2[l1] = bl1->p2[l1] - 2.0*(bl3->x2[l1])*pspr/x2cour;
  bl1->p3[l1] = bl1->p3[l1] - 2.0*(bl3->x3[l1])*pspr/x2cour;
  // fin modif surface a.b.             

 pnu602:
  std::cout <<"pnu602:" << std::endl;
  
  if(k != 0) {
    std::cout <<"k != 0" << std::endl;
    kd = 0;
    ccr = tau;
    for(G4int i = 0; i < k; i++) {
      i20 = i - kd;

      if((bl2->jnd[i] == l1) || ((bl2->ind[i] == l1) && (bl2->jnd[i] != 0))) {
	kd = kd + 1;
	continue;
      }
      
      bl2->crois[i20] = bl2->crois[i] - ccr;
      bl2->ind[i20] = bl2->ind[i];
      bl2->jnd[i20] = bl2->jnd[i];
    }
    k = k - kd;
    std::cout <<"k != 0 done..." << std::endl;
  }

  if (bl5->nesc[l1] == 1) {
    std::cout <<"bl5->nesc[l1] == 1. going to pnu613." << std::endl;
    goto pnu613;
  }

  std::cout <<"Calling new1." << std::endl;
  new1(l1);
  std::cout <<"End of new1 call..." << std::endl;
  
  std::cout <<"Calling ref." << std::endl;
  tref = ref(bl3->x1[l1],bl3->x2[l1],bl3->x3[l1],bl1->p1[l1],bl1->p2[l1],bl1->p3[l1],bl1->eps[l1],r22);
  std::cout <<"End of ref call..." << std::endl;
  
  if (tref > bl4->tmax5) {
    std::cout <<"tref > bl4->tmax5. going to pnu615." << std::endl;
    goto pnu615;
  }
  k = k + 1;
  bl2->crois[k] = tref;
  bl2->ind[k] = l1;
  bl2->jnd[k] = -1;
 pnu615:
  if (npion == 0) {
    std::cout <<"npion == 0. goig to pnu613. " << std::endl;
    goto pnu613;
  }
  if (bl1->ind1[l1] == 1) {
    std::cout <<"bl1->ind1[l1] == 1. goig to pnu613. " << std::endl;
    goto pnu613;
  }
  for(G4int k20 = 0; k20 < npion; k20++) { //do 614 k20=1,npion
    std::cout <<"Calling new3 on line 3876" << std::endl;
    new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], k20, l1);
    std::cout <<"Call to new3 complete! " << std::endl;
  }
 pnu613:
  if (k == 0) {
    std::cout <<"k == 0. going to pnu230" << std::endl;
    goto pnu230;
  }

  std::cout <<"going to pnu 449..." << std::endl;
  goto pnu449;

  // decay of the surviving deltas
 pnu230:
 pnu255:
  std::cout <<"pnu255: " << std::endl;
  if (k3 == 1) {
    std::cout <<"k3 == 1. going to pnu256" << std::endl;
    goto pnu256;
  }
  if (k4 == 0) {
    std::cout <<"k4 == 0. going to pnu256" << std::endl;
    goto pnu256;
  }

  std::cout <<"line 3911: neither of the above..." << std::endl;
  npidir = npion;
  for(G4int i = 0; i < ia; i++) {
    if (bl1->ind1[i] == 0) {
      std::cout <<"bl1->ind1[i] == 0. continuing..." << std::endl;
      continue;
    }
    npion = npion + 1;
    var_ab = pow(bl1->eps[i],2) - pow(bl1->p1[i],2) - pow(bl1->p2[i],2) - pow(bl1->p3[i],2);
    ym[npion] = 0.0;

    if(var_ab > 0.0) {
      ym[npion] = sqrt(var_ab);
    }
    xy1 = bl1->p1[i];
    xy2 = bl1->p2[i];
    xy3 = bl1->p3[i];
    xye = bl1->eps[i];
    if(varavat->kveux == 1) {
      iavat = iavat + 1;
      varavat->timeavat[iavat] = tim;
      varavat->l1avat[iavat] = i;
      varavat->l2avat[iavat] = -2;
      varavat->energyavat[iavat] = ym[npion];
      varavat->bloc_paul[iavat] = 0;
      varavat->bloc_cdpp[iavat] = 0;
      varavat->del1avat[iavat] = bl1->ind1[l1];
      varavat->jpartl1[iavat] = 1;
      varavat->jpartl2[iavat] = 0;
    }

    std::cout <<"Calling decay2..." << std::endl;
    decay2(&(bl1->p1[i]), &(bl1->p2[i]), &(bl1->p3[i]), &(bl1->eps[i]), &q1[npion], &q2[npion], &q3[npion],
	   &q4[npion], &ym[npion], &fmp, &fmpi, &(bl9->hel[i]));
    std::cout <<"End of decay2 call..." << std::endl;
    
    if(bl5->nesc[i] == 0) {
      idecf = 1;
    }

    if (ws->nosurf <= 0) {
      // surface
      if (bl5->nesc[i] == 0.) {
	pppp = sqrt(pow(bl1->p1[i],2) + pow(bl1->p2[i],2) + pow(bl1->p3[i],2));
	rrrr = sqrt(pow(bl3->x1[i],2) + pow(bl3->x2[i],2) + pow(bl3->x3[i],2));
	if (pppp <= bl10->pf) {
	  xv = pppp/bl10->pf;
	  rcorr = interpolateFunction(xv);
	  if (rrrr > rcorr) { //then
	    bl3->x1[i] = bl3->x1[i]*rcorr/rrrr;
	    bl3->x2[i] = bl3->x2[i]*rcorr/rrrr;
	    bl3->x3[i] = bl3->x3[i]*rcorr/rrrr;
	  }
	}
      }
      // fin surface
    }

    if (bl1->ind2[i]*(bl1->ind2[i]) == 9) {
      goto pnu280;
    }

    standardRandom(&rndm, &(hazard->ial));

    if (rndm*3. < 1.0) {
      goto pnu283;
    }
    ipi[npion] = 0;
    goto pnu285;

  pnu283:
    ipi[npion] = bl1->ind2[i]*2;
    bl1->ind2[i] = -1*(bl1->ind2[i]);
    goto pnu285;

  pnu280:
    bl1->ind2[i] = bl1->ind2[i]/3;
    ipi[npion] = 2*(bl1->ind2[i]);

  pnu285:
    y1[npion] = bl3->x1[i];
    y2[npion] = bl3->x2[i];
    y3[npion] = bl3->x3[i];
  }
  std::cout <<"out of the loop..." << std::endl;

  if(varavat->kveux == 1) {
    varavat->bavat = b;
    varavat->nopartavat = nopart;
    varavat->ncolavat = ncol;
    varavat->nb_avat = iavat;
  }

  // final properties of the incoming nucleon and of the remnant
  // before evaporation

  // tableau des energies a la fin (avatar.hbk)
  if(varavat->kveux == 1) {
    for(G4int i = 0; i < ia; i++) {
      if(bl5->nesc[i] == 0) {
	if(jparticip[i] == 1) {
	  varavat->epsf[i] = bl1->eps[i];
	}
	else {
	  varavat->epsf[i] = 0.0;
	}
      }
      else {
	varavat->epsf[i] = 0.0;
      }
    }
  }

 pnu256:
  std::cout <<"Starting pnu256..." << std::endl;
  elead = 0.0;
  lead = 0;
  npx = 0;
  erem = 0.;
  izrem = 0;
  inrem = 0;
  iarem = 0;
  rcm1 = 0.0;
  rcm2 = 0.0;
  rcm3 = 0.0;
  prem1 = 0.0;
  prem2 = 0.0;
  prem3 = 0.0;
  pout1 = 0.0;
  pout2 = 0.0;
  pout3 = 0.0;
  eout = 0.0;
  cmultn = 0.0;

  if (kindstruct->kindf7 <= 2) {
    if (ncol == 0 || nc[0] == 0) { // then nc(1)->nc[0]
      goto pnu9100;
    }
  }
  else {
    if (kindstruct->kindf7 <= 5) {
      if (ncol == 0) {
	goto pnu9100;
      }
    }
    else {
      // ici faisceau composite: modif a.b. 2/2002 pour tous les composites:
      nsum_col = 0;
      for(G4int i = 0; i < bl3->ia1; i++) {
	nsum_col = nsum_col + nc[i];
      }
      if (ncol == 0 || nsum_col == 0) { //then
	goto pnu9100;
      }
    }
  }

  goto pnu9101;
  // pour eviter renvoi des resultats du run precedent cv 20/11/98
 pnu9100:
  iarem = bl3->ia2;
  izrem = iz2;
  esrem = 0.0;
  erecrem = 0.0;
  nopart = -1;
  // fin ajout cv
  //return
  std::cout <<"Blah " << std::endl;
  goto pnureturn;

 pnu9101:
  nopart = 0;
  ekout = 0.0;

  for(G4int i = 0; i < ia; i++) {
    std::cout <<"Now we can get some particles out of this mess..." << std::endl;
    
    if(bl5->nesc[i] != 0) {
      xl1 = xl1-bl3->x2[i]*(bl1->p3[i]) + (bl3->x3[i])*(bl1->p2[i]);
      xl2 = xl2-bl3->x3[i]*(bl1->p1[i]) + bl3->x1[i]*(bl1->p3[i]);
      xl3 = xl3-bl3->x1[i]*(bl1->p2[i]) + bl3->x2[i]*(bl1->p1[i]);

      // ici ajout de pout cv le 5/7/95
      pout1 = pout1 + bl1->p1[i];
      pout2 = pout2 + bl1->p2[i];
      pout3 = pout3 + bl1->p3[i];
      eout = eout+bl1->eps[i] - fmp;
      ic33 = (bl1->ind2[i]+3)/2;

      nopart = nopart + 1;
      kind[nopart] = 3 - ic33;
      ep[nopart] = bl1->eps[i] - fmp;
      bmass[nopart] = fmp;
      ekout = ekout + ep[nopart];
      ptotl = sqrt(pow(bl1->p1[i],2) + pow(bl1->p2[i],2) + pow(bl1->p3[i],2));
      alpha[nopart] = bl1->p1[i]/ptotl;
      beta[nopart] = bl1->p2[i]/ptotl;
      gam[nopart] = bl1->p3[i]/ptotl;

      continue;
    }
    //  pnu254: 
    t[3] = bl3->x1[i]*(bl3->x1[i]) + bl3->x2[i]*(bl3->x2[i]) + bl3->x3[i]*(bl3->x3[i]); //t(4)->t[3]
    erem = erem + bl1->eps[i] - fmp;
    rcm1 = rcm1 + bl3->x1[i];
    rcm2 = rcm2 + bl3->x2[i];
    rcm3 = rcm3 + bl3->x3[i];
    prem1 = prem1 + bl1->p1[i];
    prem2 = prem2 + bl1->p2[i];
    prem3 = prem3 + bl1->p3[i];
    izrem = izrem + (1 + bl1->ind2[i])/2;
    iarem = iarem + 1;
  }

  //  correction pions 21/3/95 jc
  ichpion = 0;
  if(npion != 0) {
    for(G4int ipion = 0; ipion < npion; ipion++) {
      pout1 = pout1 + q1[ipion];
      pout2 = pout2 + q2[ipion];
      pout3 = pout3 + q3[ipion];
      eout = eout + q4[ipion];
      xl1 = xl1 - y2[ipion]*q3[ipion] + y3[ipion]*q2[ipion];
      xl2 = xl2 - y3[ipion]*q1[ipion] + y1[ipion]*q3[ipion];
      xl3 = xl3 - y1[ipion]*q2[ipion] + y2[ipion]*q1[ipion];
      ichpion = ichpion + ipi[ipion]/2;
      nopart = nopart + 1;
      kind[nopart] = 4 - ipi[ipion]/2;
      ptotl = sqrt(pow(q1[ipion],2) + pow(q2[ipion],2) + pow(q3[ipion],2));
      ep[nopart] = q4[ipion] - fmpi;
      bmass[nopart] = fmpi;
      ekout = ekout + ep[nopart];
      alpha[nopart] = q1[ipion]/ptotl;
      beta[nopart] = q2[ipion]/ptotl;
      gam[nopart] = q3[ipion]/ptotl;
    }
  }
  
  // fin correction pions sur impulsion et moment angulaire et charge
  // ici ajout de pfrem cv le 5/7/95
  pfrem1 = -pout1;
  pfrem2 = -pout2;
  pfrem3 = pinc - pout3;

  inrem = iarem - izrem;
  iejp = iz2 - izrem;
  iejn = bl3->ia2 - inrem - iz2;
  irem = inrem + izrem;

  // intrinsic momentum of the remnant (a.b. 05/2001): 
  // momentum of projectile minus momentum of all outgoing particles
  // minus angular momentum of the remnant computed
  // from the sum of all inside nucleons.
  //   2675	c      xl1=xl1-rcm2/irem*prem3+rcm3/irem*prem2
  //   2676	c      xl2=xl2-rcm3/irem*prem1+rcm1/irem*prem3
  //   2677	c      xl3=xl3-rcm1/irem*prem2+rcm2/irem*prem1
  //   2678	c                           here the remnant momentum is pin - sig(pout),
  //   2679	c   and the distance with respect to the barycenter of the actual target
  xl1 = xl1 - (rcm2/irem - x2_target)*pfrem3+(rcm3/irem - x3_target)*pfrem2;
  xl2 = xl2 - (rcm3/irem - x3_target)*pfrem1+(rcm1/irem - x1_target)*pfrem3;
  xl3 = xl3 - (rcm1/irem - x1_target)*pfrem2+(rcm2/irem - x2_target)*pfrem1;
  l = int(sqrt(xl1*xl1 + xl2*xl2 + xl3*xl3)/hc + 0.5);

  iej = bl3->ia2 - irem;

  eh5 = erem - pow(double(irem)/double(a2),1.666667)*efer;
  sepa = (bl3->ia2 - irem)*(v0 - tf);
  eh6 = eh5;

  // deutons ajout beproj ?????? on retire beproj (18/06/2002 ab cv)
  // eh5=erem-efer-beproj+(ia2-irem)*tf
  eh5 = erem - efer + (bl3->ia2 - irem)*tf;
  if (eh5 < 0.0) {
    eh5 = 0.00000001;
  }

  xlab = tlab - eout - eh5 - sepa;
  ecoreh5 = 0.0;

  if (iqe == 1) {
    eh5=0.00000001;
  }
  else {
    if (eh5 < 0.0) {
      if (npion == 0) {
	nopart = -1;
	std::cout <<"npion == 0" << std::endl;
	goto pnureturn;
      }
      else {
	ecoreh5 = -eh5;
	eh5 = 0.000000001;
      }
    }
  }
  if (idecf != 0 && eh5 < 0.0) {
    ecoreh5 = -eh5;
    eh5 = 0.000001;
  }

  iarem = irem;
  pfreml2 = pow(pfrem1,2) + pow(pfrem2,2) + pow(pfrem3,2);
  if (pfreml2 > 1.0e-12) {
    pfreml = sqrt(pfreml2);
    alrem = pfrem1/pfreml;
    berem = pfrem2/pfreml;
    garem = pfrem3/pfreml;
  }
  else {
    alrem = 0.0;
    berem = 0.0;
    garem = 1.0;
  }

  erecrem = pfreml2/(sqrt(pfreml2 + pow((fmp*iarem),2)) + fmp*iarem);

  if(iarem == 1) {
    erecrem = erecrem + eh5;
  }

  // correction recul
  erecg = erecrem + ecoreh5;
  // correction energie d'excitation pour une absorption (a.b., c.v. 2/2002)
  esrem = eh5;

  if (ekout < 0.001) {
    std::cout <<"ekout < 0.001" << std::endl;
    goto pnureturn;
  }
 
  // on ote l'instruction precedente car esrem toujours nulle 14/9/99

  if (erecg > 0.25) {
    fffc = (ekout - erecg)/ekout;
    if (fffc < 0.0) {
      fffc = 0.0;
    }

    for(G4int ipart = 0; ipart < nopart; ipart++) {
      ep[ipart] = ep[ipart]*fffc;
    }
  }

  // modif boudard juillet 99 (il faut tenir compte de la renormalisation
  // des energies pour les impulsions.)
  pfrem1 = 0.0;
  pfrem2 = 0.0;
  pfrem3 = pinc;
  for(G4int ipart = 0; ipart < nopart; ipart++) {
    xmodp = sqrt(ep[ipart]*(2.0*bmass[ipart] + ep[ipart]));
    pfrem1 = pfrem1 - alpha[ipart]*xmodp;
    pfrem2 = pfrem2 - beta[ipart]*xmodp;
    pfrem3 = pfrem3 - gam[ipart]*xmodp;
  }
  // fin modif a.b.

  pfreml2 = pow(pfrem1,2) + pow(pfrem2,2) + pow(pfrem3,2);
  erecrem = pfreml2/(sqrt(pfreml2 + pow((fmp*iarem),2)) + fmp*iarem);

  if (pfreml2 > 1.0e-12) {
    pfreml = sqrt(pfreml2);
    alrem = pfrem1/pfreml;
    berem = pfrem2/pfreml;
    garem = pfrem3/pfreml;
  }
  else {
    alrem = 0.0;
    berem = 0.0;
    garem = 1.0;
  }
  // fin  modif a.b. pour incl3 

  esrem = eh5;

  // if the remnant is a nucleon, no excitation energy
  if(iarem == 1) {
    esrem = 0.0;
  }

  std::cout <<"Reached end of routine..." << std::endl;

  goto pnureturn;

  std::cout <<"ia1 > 1 ! " << std::endl;
 pnureturn:
  (*ibert_p) = ibert;
  //(*f_p) = f;
  //f_p = &f;
  (*nopart_p) = nopart;
  //  (*kind_p) = kind;
//   (*ep_p) = ep;
//   (*alpha_p) = alpha; 
//   (*beta_p) = beta;
//   (*gam_p) = gam;
  (*izrem_p) = izrem;
  (*iarem_p) = iarem;
  (*esrem_p) = esrem; 
  (*erecrem_p) = erecrem;
  (*alrem_p) = alrem;
  (*berem_p) = berem;
  (*garem_p) = garem;
  (*bimpact_p) = bimpact;
  (*l_p) = l;
}


void G4Incl::collis(G4double *p1_p, G4double *p2_p, G4double *p3_p, G4double *e1_p, G4double *pout11_p, G4double *pout12_p, 
		    G4double *pout13_p, G4double *eout1_p, G4double *q1_p, G4double *q2_p, G4double *q3_p,
		    G4double *q4_p, G4int *np_p, G4int *ip_p, G4int *k2_p, G4int *k3_p, G4int *k4_p, 
		    G4int *k5_p, G4int *m1_p, G4int *m2_p, G4int *is1_p, G4int *is2_p)
{
  G4double p1 = (*p1_p);
  G4double p2 = (*p2_p);
  G4double p3 = (*p3_p);

  G4double e1 = (*e1_p);

  G4double pout11 = (*pout11_p);
  G4double pout12 = (*pout12_p);
  G4double pout13 = (*pout13_p);
  G4double eout1 = (*eout1_p);

  G4double q1 = (*q1_p);
  G4double q2 = (*q2_p);
  G4double q3 = (*q3_p);
  G4double q4 = -1*(*q4_p);

  G4int np = (*np_p);
  G4int ip = (*ip_p);

  G4int k2 = (*k2_p);
  G4int k3 = (*k3_p);
  G4int k4 = (*k4_p);
  G4int k5 = (*k5_p);

  G4int m1 = (*m1_p);
  G4int m2 = (*m2_p);

  G4int is1 = (*is1_p);
  G4int is2 = (*is2_p);

  // Variables:
  G4double a;
  G4double aaa;
  G4double aac;
  //  G4double alog;
  G4double alphac;
  //  G4double amax1;
  G4double apt;
  G4double argu;
  G4double b;
  G4double btmax;
  G4double cfi;
  G4double cpt;
  G4double ctet;
  G4double e3;
  G4double ecm;
  G4double ex[3];
  G4double ey[3];
  G4double ez[3];
  G4double f3;
  G4double f3max;
  G4double fi;
  G4double fracpn;
  G4double heli;
  G4int iexpi;
  G4int ii;
  G4int index;
  G4int index2;
  G4int isi;
  G4double pin;
  G4double pl;
  G4double pnorm;
  G4double pq;
  G4double psq;
  G4double qq[3];
  G4double qq4;
  G4double ranres;
  G4double rndm;
  G4double s;
  G4double s1;
  //  G4double sel;
  G4double sfi;
  G4double stet;
  G4double t;
  G4double x;
  G4double xkh;
  G4double xp1;
  G4double xp2;
  G4double xp3;
  G4double xx;
  G4double y;
  G4double yn;
  G4double z;
  G4double zn;
  G4double zz;

  // !!!  q4 = -1*q4;
  //   2837	c
  //   2838	c collis                                                                p-n17770
  //   2839	      subroutine collis(p1,p2,p3,e1,pout11,pout12,pout13,eout1,q1,q2,q3,p-n17780
  //   2840	     -q4,np,ip,k2,k3,k4,k5,m1,m2,is1,is2)                               p-n17790

  //   2841	      common/bl6/xx10,isa                                               p-n17800
  //   2842	      common/bl8/rathr,ramass                                           p-n17810
  //   2843	      common/hazard/ial,iy1,iy2,iy3,iy4,iy5,iy6,iy7,iy8,iy9,iy10,
  //   2844	     s               iy11,iy12,iy13,iy14,iy15,iy16,iy17,iy18,iy19
  //   2845	      common/bl9/hel(300),l1,l2
  //  2846	      dimension ex(3),ey(3),ez(3),qq(3)                                 p-n17820
  //  2847	      data xm,xm2,xmdel,ezero,xpi/938.2796,8.8037e5,1232.,1876.6,138./  p-n17830
  G4double xm = 938.2796;
  G4double xm2 = 8.8037e5;
  G4double xmdel = 1232.;
  G4double xpi = 138.;

  //   2848	c      data iy1,iy2,iy3,iy4,iy5,iy6,iy7,iy8,iy10,iy11,iy12,iy13/         p-n17840
  //   2849	c     1 12345,22345,32345,42345,52345,62345,72345,82345,34567,47059,21033p-n17850
  //   2850	c     1 12345,22345,32345,42345,52345,62345,72345,82341,34567,47059,21033p-n17850
  //   2851	c     2,32835/                                                           p-n17860
  //   2852	c      data iy9/15637/
  //  2853	      pcm(e,a,c)=0.5*sqrt((e**2-(a+c)**2)*(e**2-(a-c)**2))/e            p-n17870
  G4int iso = is1 + is2;
  np = 0;
  psq = p1*p1 + p2*p2 + p3*p3;
  pnorm = sqrt(psq);
  ecm = e1 + eout1;
  if(ecm < 1925.0) {
    goto collis160;
  }

  if (k3 == 1) {
    goto collis17;
  }
  
  if(ecm < (2065.0 + bl8->rathr)) {
    goto collis17;
  }

  if((k4-1) < 0) {
    goto collis18;
  }
  if((k4-1) == 0) {
    goto collis10;
  }
  if((k4-1) > 0) {
    goto collis20;
  }

 collis10:
  if((m1+m2-1) < 0) {
    goto collis19;
  }
  if((m1+m2-1) == 0) {
    goto collis14;
  }
  if((m1+m2-1) > 0) {
    goto collis13;
  }

 collis19:
  if((ecm-2170.4-bl8->ramass) <= 0) {
    goto collis17;
  }
  if((ecm-2170.4-bl8->ramass) > 0) {
    goto collis18;
  }

 collis20:
  if((m1+m2-1) < 0) {
    goto collis18;
  }
  if((m1+m2-1) == 0) {
    goto collis14;
  }
  if((m1+m2-1) > 0) {
    goto collis13;
  }

  // test on the recombination possibility for the n-delta system
 collis14:
  if (k5 == 0) {
    goto collis170;
  }

  standardRandom(&rndm, &(hazard->igraine[10]));
  s1 = lowEnergy(ecm, 1, iso);

  if(m1 != 0) {
    bl6->xx10=sqrt(e1*e1-psq);
    bl6->isa=is1;
  }
  else {
    bl6->xx10 = sqrt(pow(eout1,2)-psq);
    bl6->isa = is2;
  }

  s = s1 + srec(ecm,bl6->xx10,iso,bl6->isa);
  a = (s - s1)/s;

  if((rndm-a) <= 0) {
    goto collis170;
  }
  if((rndm-a) > 0) {
    goto collis17;
  }

  // test for the behaviour of the delta-delta system
 collis13:
  if (k5 == 0) {
    goto collis160;
  }
  goto collis17;

  // test on the inelasticity
 collis18:
  standardRandom(&rndm, &(hazard->igraine[0]));
  s = lowEnergy(ecm,0,iso);
  a = deltaProductionCrossSection(ecm,iso);
  a = s/(s+a);
  if(rndm > a) {
    goto collis100;
  }

  // elastic scattering
  // fit of the b parameter in the differential x-section:
  // taken from :j.c.,d.l'hote, j.vandermeulen,nimb111(1996)215
  // for pn :improvement of the backward scattering according
  // j.c et al, prc56(1997)2431
 collis17:
  pl = 0.5*ecm*sqrt(pow(ecm,2) - 4.0*xm2)/xm;
  x = 0.001*pl;
  if (iso == 0) {
    goto collis80;
  }
  if (pl > 2000.0) {
    goto collis81;
  }
  x = x*x;
  x = pow(x,4);
  b = 5.5e-6*x/(7.7+x);
  goto collis82;

 collis81:
  b = (5.34 + 0.67*(x-2.0))*1.e-6;
  goto collis82;

 collis80: 
  if (pl < 800.) {
    b = (7.16 - 1.63*x)*1.e-6;
    b = b/(1.0 + exp(-(x - 0.45)/0.05));
  }
  else {
    if (pl < 1100.) {
      b = (9.87 - 4.88*x)*1.e-6;
    }
    else {
      b = (3.68 + 0.76*x)*1.e-6;
    }
  }

 collis82:
  btmax = 4.0*psq*b;
  z = exp(-btmax);
  standardRandom(&rndm, &(hazard->igraine[1]));
  ranres = rndm; 
  y = 1.0 - rndm*(1.0 - z);
  t = log(y)/b;
  iexpi = 0;

  if (((m1+m2) == 0) && (iso == 0)) {
    apt = 1.0;                                                       
    if (pl > 800.) {
      apt = pow((800.0/pl),2);                                          
      cpt = amax1(6.23*exp(-1.79*x),0.3);
      alphac = 100.0*1.e-6;
      aaa = (1 + apt)*(1 - exp(-btmax))/b;
      argu = psq*alphac;

      if (argu >= 8) {
	argu = 0.0;
      }
      else {
	argu = exp(-4.0*argu);
      }

      aac = cpt*(1.0 - argu)/alphac;
      fracpn = aaa/(aac + aaa);
      standardRandom(&rndm, &(hazard->igraine[7]));
      if (rndm > fracpn) {
	z = exp(-4.0*psq*alphac);
	iexpi = 1;
	y = 1.0 - ranres*(10.0 - z);                                                 
	t = log(y)/alphac;
      }
    }
  }

  ctet = 1.0 + 0.5*t/psq;
  if(fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }      

  stet = sqrt(1.0 - pow(ctet,2));
  standardRandom(&rndm, &(hazard->igraine[2]));
  fi = 6.2832*rndm;
  cfi = cos(fi);
  sfi = sin(fi);
  xx = p1*p1 + p2*p2;
  zz = p3*p3;    

  if(xx >= (zz*1.0e-8)) {
    yn=sqrt(xx);
    zn=yn*pnorm;
    ez[0] = p1/pnorm; // ez(1) -> ez[0] and so on...
    ez[1] = p2/pnorm;
    ez[2] = p3/pnorm;
    ex[0] = p2/yn;
    ex[1] = -p1/yn;
    ex[3] = 0.0;     
    ey[0] = p1*p3/zn;
    ey[1] = p2*p3/zn;
    ey[2] = -xx/zn;
    p1 = (ex[0]*cfi*stet + ey[0]*sfi*stet + ez[0]*ctet)*pnorm;
    p2 = (ex[1]*cfi*stet + ey[1]*sfi*stet + ez[1]*ctet)*pnorm;
    p3 = (ex[2]*cfi*stet + ey[2]*sfi*stet + ez[2]*ctet)*pnorm;
  }
  else {
    p1 = p3*cfi*stet;
    p2 = p3*sfi*stet;
    p3 = p3*ctet;
  }
  pout11 = -p1;
  pout12 = -p2;
  pout13 = -p3;

  // backward scattering according the parametrization of ref
  // prc56(1997)1

  if(((m1+m2) != 1) || (iso == 0)) {
    standardRandom(&rndm, &(hazard->igraine[7]));
    apt = 1.0;
    if (pl > 800.0) {
      apt = pow((800.0/pl),2);
    }
    if ((iexpi == 1) || (rndm > (1./(1.+apt)))) {
      ii = is1;
      is1 = is2;
      is2 = ii;
    }
  }

  goto exitRoutine;

  // delta production
  // the production is not isotropic in this version
  // it has the same exp(b*t) structure as the nn elastic scatteringp-n19170
  // (formula 2.3 of j.cugnon et al, nucl phys a352(1981)505)
  // parametrization of b taken from ref. prc56(1997)2431
 collis100:
  if (k4 != 1) {
    goto collis101;
  }
  xmdel = 1232.0 + bl8->ramass;
  goto collis103;
 collis101: 
  //call ribm(rndm,iy10)
  standardRandom(&rndm, &(hazard->igraine[9]));

  y = tan(3.1415926*(rndm - 0.5));
  x = 1232.0 + 0.5*130.0*y + bl8->ramass;
  if (x < (xm+xpi+2.0)) {
    goto collis101;
  }
  if (ecm < (x+xm+1.)) {
    goto collis101;
  }

  // generation of the delta mass with the penetration factor
  // (see prc56(1997)2431)
  y = pow(ecm,2);
  q2 = (y - pow(1076.0,2))*(y - pow(800.0,2))/y/4.0;                                 
  q3 = pow((sqrt(q2)),3);                                                  
  f3max = q3/(q3 + pow(180.0,3));                                      
  y = pow(x,2);
  q2 = (y - pow(1076.0,2))*(y - pow(800.0,2))/y/4.0;                                 
  q3 = pow((sqrt(q2)),3);                                                 
  f3 = q3/(q3 + pow(180.0,3));

  standardRandom(&rndm, &(hazard->igraine[10]));
  if (rndm > (f3/f3max)) {
    goto collis101;
  }                          
  xmdel = x;

 collis103:
  pin = pnorm;
  pnorm = pcm(ecm,xm,xmdel);
  if (pnorm <= 0) {
    pnorm = 0.000001;
  }
  index = 0;
  index2 = 0;

  standardRandom(&rndm, &(hazard->igraine[3]));
  if (rndm < 0.5) {
    index = 1;
  }
  
  if (iso == 0) {
    standardRandom(&rndm, &(hazard->igraine[4]));
    if (rndm < 0.5) {
      index2 = 1;
    } 
  }

  standardRandom(&rndm, &(hazard->igraine[5]));
  x = 0.001*0.5*ecm*sqrt(pow(ecm,2) - 4.0*xm2)/xm;
  if(x < 1.4) {
    b = (5.287/(1.0 + exp((1.3 - x)/0.05)))*1.e-6;
  }
  else {
    b = (4.65 + 0.706*(x - 1.4))*1.e-6;
  }

  xkh = 2.0*b*pin*pnorm;
  ctet=1.0 + log(1.0 - rndm*(1.0 - exp(-2.0*xkh)))/xkh;
  if(fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }
  
  stet = sqrt(1.0 - pow(ctet,2));
  standardRandom(&rndm, &(hazard->igraine[6]));
  fi = 6.2832*rndm;
  cfi = cos(fi);
  sfi = sin(fi);

  // delta production: correction of the angular distribution 02/09/02
  xx = p1*p1 + p2*p2;
  zz = p3*p3;
  if(xx >= (zz*1.0e-8)) {
    yn = sqrt(xx);
    zn = yn*pin;
    ez[0] = p1/pin; // ez(1) -> ez[0] and so on...
    ez[1] = p2/pin;
    ez[2] = p3/pin;
    ex[0] = p2/yn;
    ex[1] = -p1/yn;
    ex[2] = 0.0;
    ey[0] = p1*p3/zn;
    ey[1] = p2*p3/zn;
    ey[2] = -xx/zn;
    xp1 = (ex[0]*cfi*stet + ey[0]*sfi*stet + ez[0]*ctet)*pnorm;
    xp2 = (ex[1]*cfi*stet + ey[1]*sfi*stet + ez[1]*ctet)*pnorm;
    xp3 = (ex[2]*cfi*stet + ey[2]*sfi*stet + ez[2]*ctet)*pnorm;
  }
  else {
    xp1 = pnorm*stet*cfi;
    xp2 = pnorm*stet*sfi;
    xp3 = pnorm*ctet;
  }
  // end of correction angular distribution of delta production

  e3 = sqrt(xp1*xp1 + xp2*xp2 + xp3*xp3 + xm*xm);
  if(k4 != 0) {
    goto collis161;
  }
  
  // decay of the delta particle (k4=0)
  np = 1;
  ip = 0;
  qq[0] = xp1; //qq(1) -> qq[0]
  qq[1] = xp2;
  qq[2] = xp3;
  qq4 = sqrt(xp1*xp1 + xp2*xp2 + xp3*xp3 + xmdel*xmdel);
  heli = pow(ctet,2);
  decay2(&qq[0],&qq[1],&qq[2],&qq4,&q1,&q2,&q3,&q4,&xmdel,&xm,&xpi,&heli);

  if(index != 0) {
    p1 = qq[0]; //qq(1) -> qq[0] and so on...
    p2 = qq[1];
    p3 = qq[2];
    pout11 = -xp1;
    pout12 = -xp2;
    pout13 = -xp3;
    eout1 = e3;
  }
  else {
    pout11 = qq[0]; //qq(1) -> qq[0] and so on...
    pout12 = qq[1];
    pout13 = qq[2];
    eout1 = e1;
    p1 = -xp1;
    p2 = -xp2;
    p3 = -xp3;
    e1 = e3;
  }

  if (iso == 0) {
    goto collis150;
  }
  if (rndm > 0.333333) {
    goto exitRoutine;
  }

  is1 = -is1;
  ip = -2*is1;

 collis150:
  if (index == 1) {
    goto exitRoutine;
  }
  if (rndm < 0.5) {
    goto collis152;
  }
  is1 = 1;
  is2 = 1;
  ip = -2;
  goto exitRoutine;

 collis152: 
  is1 = -1;
  is2 = -1;
  ip = 2;
  goto exitRoutine;

 collis160:
  pout11 = -p1;
  pout12 = -p2;
  pout13 = -p3;
  goto exitRoutine;

  // long-lived delta
 collis161:
  if(index != 1) {
    p1=xp1;
    p2=xp2;
    p3=xp3;
    eout1=e3;
    e1=ecm-eout1;
    m1=1;
  }
  else {
    p1=-xp1;
    p2=-xp2;
    p3=-xp3;
    eout1=e3;
    e1=ecm-eout1;
    m1=1;
  }

  // symmetrization of charges in pn -> n delta
  // the test on "index" above symetrizes the excitation of one 
  // of the nucleons with respect to the delta excitation
  // (see note 16/10/97)
  if (iso == 0) {
   if (index2 == 1) {
     isi = is1;
     is1 = is2;
     is2 = isi;
   }
   goto collis160;
  }

  bl9->hel[bl9->l1] = pow(ctet,2);
  standardRandom(&rndm, &(hazard->igraine[7]));
  if (rndm < 0.25) {
    goto collis160;
  }

  is1=3*is1*m1-(1-m1)*is1;
  is2=3*is2*m2-(1-m2)*is2;
  goto collis160;

  // recombination process
 collis170: 
  pnorm = pcm(ecm,xm,xm);
  standardRandom(&rndm, &(hazard->igraine[11]));
  ctet = -1.0 + 2.0*rndm;
  if(fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }
  stet = sqrt(1.0 - ctet*ctet);
  standardRandom(&rndm, &(hazard->igraine[12]));
  fi = 6.2832*rndm;
  cfi = cos(fi);
  sfi = sin(fi);
  p1 = pnorm*stet*cfi;
  p2 = pnorm*stet*sfi;
  p3 = pnorm*ctet;
  m1 = 0;
  m2 = 0;
  e1 = sqrt(p1*p1 + p2*p2 + p3*p3 + xm*xm);
  eout1 = ecm - e1;

  if (iso == 0) {
    goto collis160;
  }
  is1=iso/2;
  is2=iso/2;
  goto collis160;
  
 exitRoutine:
  (*p1_p) = pq;
  (*p2_p) = p2;
  (*p3_p) = p3;

  (*e1_p) = e1;

  (*pout11_p) = pout11;
  (*pout12_p) = pout12;
  (*pout13_p) = pout13;
  (*eout1_p) = eout1;

  (*q1_p) = q1;
  (*q2_p) = q2;
  (*q3_p) = q3;
  (*q4_p) = q4; // :::FIXME::: (*q4_p) = -1*q4;

  (*np_p) = np;
  (*ip_p) = np;

  (*k2_p) = k2;
  (*k3_p) = k2;
  (*k4_p) = k4;
  (*k5_p) = k5;

  (*m1_p) = m1;
  (*m2_p) = m2;

  (*is1_p) = is1;
  (*is2_p) = is2;
}

void G4Incl::decay2(G4double *p1_p, G4double *p2_p, G4double *p3_p, G4double *wp_p, G4double *q1_p, 
		    G4double *q2_p, G4double *q3_p, G4double *wq_p, G4double *xi_p, G4double *x1_p, G4double *x2_p, 
		    G4double *hel_p)
{            
  // This routine describes the anisotropic decay of a particle of mass
  // xi into 2 particles of masses x1,x2                             
  // the anisotropy is supposed to follow a 1+3*hel*(cos(theta))**2
  // law with respect to the direction of the incoming particle
  // in the input, p1,p2,p3 is the momentum of particle xi
  // in the output, p1,p2,p3 is the momentum of particle x1 , while
  // q1,q2,q3 is the momentum of particle x2

  // Temporary variables for input/output data:
  G4double p1 = (*p1_p);
  G4double p2 = (*p2_p);
  G4double p3 = (*p3_p);

  G4double q1 = (*q1_p);
  G4double q2 = (*q2_p);
  G4double q3 = (*q3_p);

  G4double wp = (*wp_p);
  G4double wq = (*wq_p);

  G4double xi = (*xi_p);
  G4double x1 = (*x1_p);
  G4double x2 = (*x2_p);
  G4double hel = (*hel_p);

  G4double rndm;

  G4double xe = wp;
  G4double b1=p1/xe;
  G4double b2=p2/xe;
  G4double b3=p3/xe;                                                          
  G4double xq=pcm(xi,x1,x2);                                                  
  G4double ctet, stet;

  G4double fi, cfi, sfi;
  G4double sal, cal;
  G4double t1, t2;
  G4double w1;
  G4double beta;

  do {
    standardRandom(&rndm, &(hazard->igraine[7]));
    ctet = -1.0 + 2.0*rndm;
    if(fabs(ctet) > 1.0) {
      ctet = sign(1.0,ctet);
    }
    stet = sqrt(1.0 - pow(ctet,2));
    standardRandom(&rndm, &(hazard->igraine[9]));
  } while(rndm > ((1.0 + 3.0*hel*pow(ctet,2))/(1.0 + 3.0*hel)));
  
  standardRandom(&rndm, &(hazard->igraine[8]));
  fi = 6.2832*rndm;
  cfi = cos(fi);
  sfi = sin(fi);
  beta = sqrt(b1*b1+b2*b2+b3*b3);

  sal = sqrt(pow(b1,2) + pow(b2,2))/beta;
  cal = b3/beta;

  if((beta >= 1.0e-10) || (sal >= 1.0e-6)) {
    t1 = ctet + cal*stet*sfi/sal;
    t2 = stet/sal;                                                       
    q1 = xq*(b1*t1 + b2*t2*cfi)/beta;
    q2 = xq*(b2*t1 - b1*t2*cfi)/beta;
    q3 = xq*(b3*t1/beta - t2*sfi);
  }
  else {
    q1 = xq*stet*cfi;                                                    
    q2 = xq*stet*sfi;                                                    
    q3 = xq*ctet;
  }
  
  hel = 0.0;                                                       
  w1 = q1*q1 + q2*q2 + q3*q3;
  wq = sqrt(w1 + x2*x2);
  p1 = -q1;
  p2 = -q2;
  p3 = -q3;
  wp = sqrt(w1+x1*x1);
  loren(&q1, &q2, &q3, &b1, &b2, &b3, &wq);
  loren(&p1, &p2, &p3, &b1, &b2, &b3, &wq);

  // Return calculated values:
  (*p1_p) = p1;
  (*p2_p) = p2;
  (*p3_p) = p3;

  (*q1_p) = q1;
  (*q2_p) = q2;
  (*q3_p) = q3;

  (*wp_p) = wp;
  (*wq_p) = wq;

  (*xi_p) = xi;
  (*x1_p) = x1;
  (*x2_p) = x2;
  (*hel_p) = hel;
}

void G4Incl::time(G4int i, G4int j)
{
  // time 
  G4double t[10];

  t[0] = bl1->p1[i]/bl1->eps[i] - bl1->p1[j]/bl1->eps[j]; // t(1)->t[0] 
  t[1] = bl1->p2[i]/bl1->eps[i] - bl1->p2[j]/bl1->eps[j]; // t(2)->t[1] and so on ...
  t[2] = bl1->p3[i]/bl1->eps[i] - bl1->p3[j]/bl1->eps[j]; 
  t[3] = bl3->x1[i] - bl3->x1[j];
  t[4] = bl3->x2[i] - bl3->x2[j];
  t[5] = bl3->x3[i] - bl3->x3[j];
  t[6] = t[0]*t[3] + t[1]*t[4] + t[2]*t[5];
  t[9] = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];

  if(t[9] <= 1.0e-10) {
    bl1->ta = 100000;
  }
  else {	
    bl1->ta = -t[6]/t[9];
  }

  bl3->rab2 = t[3]*t[3] + t[4]*t[4] + t[5]*t[5] + bl1->ta*t[6];
}

void G4Incl::newt(G4int l1, G4int l2)
{
  G4int ig, id, kg, kd;
  G4int iy, ix;
  G4double E;

  G4int ia = bl3->ia1 + bl3->ia2;
  for(G4int i = 0; i < ia; i++) { // do 52 i=1,ia
    if (bl5->nesc[i] != 0) {
      continue;
    }
    //    if (i-l2) 53,52,54
    if((i-l2) < 0) {
      goto newt53;
    }
    if((i-l2) == 0) {
      continue;
    }
    if((i-l2) > 0) {
      goto newt54;
    }
  newt53: 
    ig=l2;
    id=i;
    kg=l1;
    kd=i;
    goto newt55;
  newt54:
    //      if (i-l1) 56,52,57
    if((i-l1) == 0) {
      continue;
    }
    if((i-l1) < 0) {
      goto newt56;
    }
    if((i-l1) > 0) {
      goto newt57;
    }
  newt56:
    kg=l1;
    kd=i; 
    goto newt58;
  newt57:
    kg=i;
    kd=l1;
  newt58:
    ig=i;
    id=l2;
  newt55:
    //  call time(ig,id)
    time(ig, id);
    if (bl1->ta < 0.) {
      goto newt50;
    }
    if(bl1->ta > bl4->tmax5) {
      goto newt50;
    }
    if (bl1->ta < bl5->tlg[l1]) { // tlg(12)->tlg[11]
      goto newt50;
    }
    if ((bl1->ind1[ig]+bl1->ind1[id]) > 0) {
      goto newt60;
    }
    E=am(bl1->p1[ig]+bl1->p1[id],bl1->p2[ig]+bl1->p2[id],bl1->p3[ig]+bl1->p3[id],bl1->eps[ig]+bl1->eps[id]);
    if (E < 1925.) {
      goto newt50;
    }
    iy=bl1->ind1[ig]+bl1->ind1[id];
    if (iy != 1) {
      goto newt61;
    }
    ix=ig*(bl1->ind1[ig])+id*(bl1->ind1[id]);
    bl6->xx10=am(bl1->p1[ix],bl1->p2[ix],bl1->p3[ix],bl1->eps[ix]);
    bl6->isa=bl1->ind2[ix];
  newt61:
    if ((31.*(bl3->rab2)) > totalCrossSection(E,iy,bl1->ind2[ig]+bl1->ind2[id])) {
      goto newt50;
    }
  newt60: // continue
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=ig;
    bl2->jnd[bl2->k]=id;
  newt50:
    time(kg,kd);
    if (bl1->ta < 0.) {
      continue;
    }
    if(bl1->ta > bl4->tmax5) {
      continue;
    }
    if (bl1->ta < bl5->tlg[10]) { //tlg(11)->tlg[10]
      continue;
    }
    if ((bl1->ind1[kg]+bl1->ind1[kd]) > 0) {
      goto newt62;
    }
    E=am(bl1->p1[kg]+bl1->p1[kd],bl1->p2[kg]+bl1->p2[kd],bl1->p3[kg]+bl1->p3[kd],bl1->eps[kg]+bl1->eps[kd]);
    if (E < 1925.) {
      continue;
    }
    iy=bl1->ind1[kg]+bl1->ind1[kd];
    if (iy != 1) {
      goto newt63;
    }
    ix=kg*(bl1->ind1[kg])+kd*(bl1->ind1[kd]);
    bl6->xx10=am(bl1->p1[ix],bl1->p2[ix],bl1->p3[ix],bl1->eps[ix]);
    bl6->isa=bl1->ind2[ix];
  newt63: 
    if ((31.*(bl3->rab2)) > totalCrossSection(E,iy,bl1->ind2[kg]+bl1->ind2[kd])) {
      continue;
    }
  newt62:
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=kg;
    bl2->jnd[bl2->k]=kd;
  }
}

void G4Incl::new1(G4int l1)
{
  G4int ia, iy, ix;
  G4double E;

  ia=bl3->ia1+bl3->ia2;
  for(G4int i = 0; i < ia; i++) {
    if (bl5->nesc[i] != 0) {
      continue;
    }
    //if(i-l1) 53,52,54
    if((i-l1) < 0) {
      goto new153;
    }
    if((i-l1) == 0) {
      continue;
    }
    if((i-l1) > 0) {
      goto new154;
    }
  new153: 
    time(i, l1);
    if(bl1->ta < 0.) {
      continue;
    }
    if(bl1->ta > bl4->tmax5) {
      continue;
    }
    if (bl1->ind1[i]+bl1->ind1[l1] > 0) {
      goto new160;
    }
    E=am(bl1->p1[i]+bl1->p1[l1],bl1->p2[i]+bl1->p2[l1],bl1->p3[i]+bl1->p3[l1],bl1->eps[i]+bl1->eps[l1]);
    if (E < 1925.) {
      continue;
    }
    iy=bl1->ind1[i]+bl1->ind1[l1];
    if (iy != 1) {
      goto new161;
    }
    ix=i*(bl1->ind1[i])+l1*(bl1->ind1[l1]);
    bl6->xx10=am(bl1->p1[ix],bl1->p2[ix],bl1->p3[ix],bl1->eps[ix]);
    bl6->isa=bl1->ind2[ix];
  new161: 
    if ((31.*(bl3->rab2)) > totalCrossSection(E ,iy,bl1->ind2[i]+bl1->ind2[l1])) {
      continue;
    }
  new160: //continue
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=l1;
    bl2->jnd[bl2->k]=i;
    continue;
  new154:
    time(i, l1);
    if(bl1->ta < 0.) {
      continue;
    }
    if(bl1->ta > bl4->tmax5) {
      continue;
    }
    if ((bl1->ind1[i]+bl1->ind1[l1]) > 0) {
      goto new170;
    }
    E=am(bl1->p1[i]+bl1->p1[l1],bl1->p2[i]+bl1->p2[l1],bl1->p3[i]+bl1->p3[l1],bl1->eps[i]+bl1->eps[l1]);
    if (E < 1925.) {
      continue;
    }
    iy=bl1->ind1[i]+bl1->ind1[l1];
    if (iy != 1) {
      goto new171;
    }
    ix=i*(bl1->ind1[i])+l1*(bl1->ind1[l1]);
    bl6->xx10=am(bl1->p1[ix],bl1->p2[ix],bl1->p3[ix],bl1->eps[ix]);
    bl6->isa=bl1->ind2[ix];
  new171: 
    if ((31.*(bl3->rab2)) > totalCrossSection(E,iy,bl1->ind2[i]+bl1->ind2[l1])) {
      continue;
    }
  new170: 
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=i;
    bl2->jnd[bl2->k]=l1;
  }
}

void G4Incl::new2(G4double y1, G4double y2, G4double y3, G4double q1, G4double q2, G4double q3, 
		  G4double q4, G4int npion, G4int l1)
{
  G4double t[10];

  G4int ia = bl3->ia1 + bl3->ia2;
  for(G4int i = 0; i < ia; i++) {
    if (bl5->nesc[i] != 0) {
      continue;
    }
    if(i == l1) {
      continue;
    }
    if(bl1->ind1[i] == 1) {
      continue;
    }

    t[0] = bl1->p1[i]/bl1->eps[i] - q1/q4;
    t[1] = bl1->p2[i]/bl1->eps[i] - q2/q4;
    t[2] = bl1->p3[i]/bl1->eps[i] - q3/q4;
    t[3] = bl3->x1[i] - y1;
    t[4] = bl3->x2[i] - y2;
    t[5] = bl3->x3[i] - y3;
    t[6] = t[0]*t[3] + t[1]*t[4] + t[2]*t[5];

    if(t[6] > 0.0) {
      continue;
    }

    t[9] = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
    bl1->ta = -1 * t[6]/t[9];
    if(bl1->ta > bl4->tmax5) {
      continue;
    }
    G4double xx2 = t[3]*t[3] + t[4]*t[4] + t[5]*t[5] + (bl1->ta)*t[6];
    G4double E = sqrt(pow((bl1->eps[i]+q4),2) - pow((bl1->p1[i]+q1),2) - pow((bl1->p2[i]+q2),2) - pow((bl1->p3[i]+q3),2));
    if ((31.0*xx2) > pionNucleonCrossSection(E)) {
      continue;
    }
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=ia+npion;
    bl2->jnd[bl2->k]=i;
  }
}

void G4Incl::new3(G4double y1, G4double y2, G4double y3, G4double q1, G4double q2, G4double q3, 
		  G4double q4, G4int npion, G4int l1)
{
  G4double t[10];
  G4double E, xx2;
  G4int ia;

  if(bl5->nesc[l1] > 0) {
    return;
  }

  t[0] = bl1->p1[l1]/bl1->eps[l1] - q1/q4;
  t[1] = bl1->p2[l1]/bl1->eps[l1] - q2/q4;
  t[2] = bl1->p3[l1]/bl1->eps[l1] - q3/q4;
  t[3] = bl3->x1[l1] - y1;
  t[4] = bl3->x2[l1] - y2;
  t[5] = bl3->x3[l1] - y3;
  t[6] = t[0]*t[3] + t[1]*t[4] + t[2]*t[5];

  if(t[6] > 0.0) {
    return;
  }

  t[9] = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  bl1->ta = -1 * t[6]/t[9];
  if(bl1->ta > bl4->tmax5) {
    return;
  }
  if (bl1->ta < bl5->tlg[l1]) {
    return;
  }
  xx2 = t[3]*t[3] + t[4]*t[4] + t[5]*t[5] + (bl1->ta)*t[6];
  E = sqrt(pow((bl1->eps[l1]+q4),2) - pow((bl1->p1[l1]+q1),2) - pow((bl1->p2[l1]+q2),2) - pow((bl1->p3[l1]+q3),2));
  if ((31.0*xx2) > pionNucleonCrossSection(E)) {
    return;
  }
  
  bl2->k = bl2->k + 1;
  bl2->crois[bl2->k] = bl1->ta;
  ia = bl3->ia1+bl3->ia2;
  bl2->ind[bl2->k] = ia + npion;
  bl2->jnd[bl2->k] = l1;
}

void G4Incl::loren(G4double *q1, G4double *q2, G4double *q3, G4double *b1, G4double *b2, G4double *b3, G4double *E)
{
  // Transforms momentum q and energy E from a frame moving with
  // velocity beta

  G4double bb2 = (*b1) * (*b1) + (*b2) * (*b2) + (*b3) * (*b3);
  G4double bq = (*b1) * (*q1) + (*b2) * (*q2) + (*b3) * (*q3);
  G4double gam2 = 1.0/(1.0 - bb2);
  G4double gam = sqrt(gam2);
  G4double c = gam2/(gam + 1.0);
  G4double g = c * bq + gam*(*E);
  (*E) = gam * ((*E) + bq);
  (*q1) = (*q1) + (*b1) * g;
  (*q2) = (*q2) + (*b2) * g; 
  (*q3) = (*q3) + (*b3) * g;
}

G4double G4Incl::pauliBlocking(G4int l, G4double xr, G4double pr)
{
  //   G4int l = (*l_p);
  //   G4double xr = (*xr_p);
  //   G4double pr = (*pr_p);
  //   G4double f = (*f_p);

  // This subroutine calculates the occupation in phase space around
  // nucleon l , by counting the particles in a volume around l the
  // volume is the product of a sphere of radius xr in r-space by a
  // sphere of radius pr in momentum space average is taken on the spin
  // only

  //   3756	      common/bl1/p1(300),p2(300),p3(300),eps(300),ind1(300),ind2(300),tap-n24950
  //   3757	      common/bl3/r1,r2,x1(300),x2(300),x3(300),ia1,ia2,rab2             p-n24960
  //   3758	      common/bl5/tlg(300),nesc(300)                                     p-n24970
  //   3759	      common/saxw/ xx(30,500),yy(30,500),ss(30,500),nbpG4inter,imat
  //   3760	      common/ws/r0,adif,rmaxws,drws,nosurf,xfoisa,npaulstr,bmax

  G4double pmod, pr2;
  G4double xr2, rdeq, dx2, dp2;
  G4double rs, vol;
  G4int nl;
  G4int ia;

  if (ws->npaulstr == 2) {
    return 0.0;
  }

  if (ws->npaulstr == 1) {
    // pauli strict
    pmod = sqrt(pow(bl1->p1[l],2) + pow(bl1->p2[l],2) + pow(bl1->p3[l],2));
    if (pmod < 270.0) {
      return 1.0;
    }
    else {
      return 0.0;
    }
  }
  else {
    // Statistic Pauli blocking
    xr2 = xr*xr;
    pr2 = pr*pr;
    vol = pow((40.0*3.1415926/3.0),2) * (pow((xr*pr)/(2.0*3.1415926*197.13),3));
    rs = sqrt(bl3->x1[l]*bl3->x1[l] + bl3->x2[l]*bl3->x2[l] + bl3->x3[l]*bl3->x3[l]);
    if (ws->nosurf <= 0) {
      // modifs a.b.: r2 -> rmaxws pour la densite en w.s.
      rdeq = ws->rmaxws;
    }
    else {
      rdeq = ws->r0;
    }

    if ((rs - xr) <= rdeq) {
      if ((rs + xr) > rdeq) {
	vol = vol*0.5*(rdeq-rs+xr)/xr;
      }

      ia = bl3->ia1 + bl3->ia2;
      nl = 0;

      for(G4int i = 0; i < ia; i++) {
	dx2 = pow((bl3->x1[l]-bl3->x1[i]),2) + pow((bl3->x2[l]-bl3->x2[i]),2) + pow((bl3->x3[l]-bl3->x3[i]),2);
	dp2 = pow((bl1->p1[l]-bl1->p1[i]),2) + pow((bl1->p2[l]-bl1->p2[i]),2) + pow((bl1->p3[l]-bl1->p3[i]),2);
      
	if((bl5->nesc[i] > 0) || (bl1->ind1[i] > 0) || (bl1->ind2[i] != bl1->ind2[l]) || (dx2 > xr2) || (dp2 > pr2)) {
	  if(((nl - 1)/vol/2.0) > 1.0) {
	    return 1.0;
	  }
	  else {
	    return ((nl - 1)/vol/2.0);
	  }
	}
	nl = nl + 1;
      }
    }
    else {
      return 0.0;
    }
  }

  return 0.0; // The algorithm is not supposed to reach this point.
}

G4double G4Incl::lowEnergy(G4double E, G4double m, G4double i)
{
  // fit by j.vandermeulen
  // low enrgy fit of j.c., d. l'hote, j. vdm, nim b111(1996)215
  // i = 2,0,-2  for pp,pn,nn
  // m = 0,1,2 for nucleon-nucleon,nucleon-delta,delta,delta

  G4double scale = 1.0;
  G4double plab = E*sqrt(E*E-3.52e6)/1876.6;
  G4double p1 = 0.001*plab;
  G4double alp;

  if(plab > 2000.0) {
    // goto sel13;
    // sel13:
    return ((77.0/(p1 + 1.5))*scale);
  }

  if((m-1) < 0) {
    if (i == 0) {
      if (plab < 800.0) {
	if (plab < 450.0) {
	  alp = log(p1);
	  return(6.3555*exp(-3.2481*alp - 0.377*alp*alp));
	}
	else {
	  return((33.0 + 196.0*sqrt(fabs(pow((p1 - 0.95),5))))*scale);
	}
      }
      else {
	return(31.0/sqrt(p1)*scale);
      }
    }
  }

  if (plab < 800.0) {
    if (plab < 440.0) {
      return(34.0*pow((p1/0.4),(-2.104)));
    }
    else {
      return((23.5 + 1000.*pow((p1 - 0.7),4))*scale);
    }
  }
  else if(plab > 2000.0) {
    return ((77.0/(p1 + 1.5))*scale);
  }
  else {
    return((1250.0/(50.0 + p1) - 4.0*pow((p1 - 1.3),2))*scale);
  }
}

G4double G4Incl::totalCrossSection(G4double E, G4int m, G4int i)
{
  // total cross-sections
  // i=2,0,-2  for pp,pn,nn
  // m=0,1,2 for nucleon-nucleon,nucleon-delta,delta,delta

  G4double stotResult;
  G4double sine;

  if((m-1) < 0) {
    sine = deltaProductionCrossSection(E,i);
  }

  if((m-1) == 0) {
    sine = srec(E,(bl6->xx10),i,(bl6->isa));
  }

  if((m-1) > 0) {
    sine = 0.0;
  }

  stotResult = sine + lowEnergy(E,m,i);
  return stotResult;
}

G4double G4Incl::srec(G4double Ein, G4double d, G4int i, G4int isa)
{
  G4double E = Ein;
  G4double s;
  G4double x, y;
  G4double srecResult;

  if (i*i == 16) {
    return 0.0;
  }

  if(E <= (938.3 + d)) {
    return 0.0;
  }
  else {
    if(E < (938.3 + d + 2.0)) {
      E = 938.3 + d + 2.0;
    }
    s = E*E;
    x = (s - 3.523e6)/(s - pow((938.3 + d),2));
    y = s/(s - pow((d - 938.3),2));
    srecResult = 0.5*x*y*deltaProductionCrossSection(E, i);
    srecResult = srecResult*(32.0 + i*i*(isa*isa - 5))/64.0;
    srecResult = srecResult/(1.0 + 0.25*i*i);
    srecResult = 3.0*srecResult;  //pi absorption increased also for internal pions (7/3/01)

    return srecResult;
  }
}

G4double G4Incl::deltaProductionCrossSection(G4double E, G4int i)
{
  // delta production cross-sections
  // fit by j.vandermeulen
  // i = 2,0,-2  for pp,pn,nn

  G4double scali = 1.0;
  G4double plab;
  G4double p1;
  G4double sproResult;

  G4double EE = E - bl8->rathr;

  if(EE*EE-3.53e6 < 0) {
    return 0.0;
  }
  else {
    plab = EE*sqrt(EE*EE-3.52e6)/1876.6;
    p1 = 0.001*plab;
    if (plab > 800.0) {
      //goto spro1;
      //  spro1:
      if (i*i == 4) {
	//	goto spro10;
	//spro10:
	if (plab < 2000.0) {
	  //    goto spro11;
	  //  spro11:
	  if (plab < 1500.0) {
	    //      goto spro12;
	    // spro12:
	    sproResult = 23.5 + 24.6/(1.0 + exp(-10.0*p1 + 12.0)) - 1250.0/(p1+50.0)+4.0*pow((p1-1.3),2);
	    return (sproResult*scali);
	  }
	  else {
	    sproResult = 41.0 + 60.0*(p1 - 0.9)*exp(-1.2*p1) - 1250.0/(p1+50.0) + 4.*pow((p1 - 1.3),2);
	    return (sproResult*scali);
	  }
	}
	else {
	  return ((41.0 + (60.0*p1 - 54.0)*exp(-1.2*p1) - 77.0/(p1 + 1.5))*scali);
	}
      }
      else {
	if (plab < 2000.0) {
	  //goto spro2;
	  //  spro2:
	  if (plab < 1000.0) {
	    //goto spro3;
	    //  spro3:
	    return ((33.0 + 196.0*sqrt(pow(fabs(p1 - 0.95),5)) - 31.1/sqrt(p1))*scali);
	  }
	  else {
	    return ((24.2 + 8.9*p1 - 31.1/sqrt(p1))*scali);
	  }
	}
	else {
	  return ((42.-77./(p1+1.5))*scali);
	}
      }
    }
    // plab <= 800.0
    else {
      return 0.0;
    }
  }
}

G4double G4Incl::pionNucleonCrossSection(G4double x)
{
  //   sigma(pi+ + p) in the (3,3) region
  //   new fit by j.vandermeulen + constant value above the (3,3)
  //   resonance

  G4double y = x*x;
  G4double q2 = (y-pow(1076.,2))*(y-pow(800.,2))/y/4.;
  G4double q3, f3;
  G4double spn;

  if(q2 <= 0) {
    return 0.0;
  }
  else {
    q3 = pow((sqrt(q2)),3);
    f3 = q3/(q3+pow(180.0,3));
    spn = 326.5/(pow(((x - 1215.0 - bl8->ramass)*2.0/110.0),2)+1.0);
    spn = spn*(1.0 - 5.0 * (bl8->ramass/1215.0));
    return (spn*f3);
  }
}

G4double G4Incl::transmissionProb(G4double E, G4double iz, G4double izn, G4double r, G4double v0)
{
  // transmission probability for a nucleon of kinetic energy
  // E on the edge of the well of depth v0 (nr approximation)
  // iz is the isospin of the nucleon,izn the instanteneous charge
  // of the nucleus and r is the target radius

  G4double x;
  G4double barr = 0.0;

  // We need enough energy to escape from the potential well.
  if (E > v0) {
    x = sqrt(E*(E-v0));
    barr = 4.*x/(E+E-v0+x+x);
    if (iz > 0) {
      G4double b = izn*1.44/r;
      G4double px = sqrt((E-v0)/b);

      if (px < 1.0) {
	G4double g = izn/137.03*sqrt(2.*938.3/(E-v0))*(acos(px)-px*sqrt(1.-px*px));
	if (g > 35.){
	  barr=0.0;
	}
	else {
	  barr = barr*exp(-2.*g);
	}
	return barr;
      }
      else {
	return barr;
      }
    }
    else {
      return barr;
    }
  }
  else {
    return barr;
  }
}

// Saxw, Ws, Bl10
G4double G4Incl::ref(G4double x1, G4double x2, G4double x3, G4double p1, G4double p2, G4double p3, G4double E, G4double r2)
{
  const G4double  pf = 270.339 , pf2 = 73083.4;

  // surface : modif de ref
  // ref=time necessary for a nucleon to reach the surface

  G4double t2 = p1*p1+p2*p2+p3*p3;
  G4double p = sqrt(t2);
  G4double r = r2;
  G4double xv;
  G4double t1, t3, t4, t5;
  G4double s;

  if (ws->nosurf <= 0) {
    //  modif pour w.s.:
    xv = p/pf;
    r = interpolateFunction(xv);
    r = r*r;
    if (t2 > pf2) {
      r = pow(ws->rmaxws,2);
    }
  }

  do {
    t4 = x1*x1+x2*x2+x3*x3;
    if (t4 > r) {
      std::cout <<"r = " << r << "t4 = " << t4 << std::endl;
      s = sqrt(r*0.99/t4);
      x1 = x1*s;
      x2 = x2*s;
      x3 = x3*s;
    }
  }
  while(t4 > r);
  
  t1 = x1*p1+x2*p2+x3*p3;
  t3 = t1/t2;
  t5 = t3*t3+(r-t4)/t2;   

  std::cout <<"ref: t1 = " << t1 << std::endl;
  std::cout <<"ref: t2 = " << t2 << std::endl;
  std::cout <<"ref: t3 = " << t3 << std::endl;
  std::cout <<"ref: t4 = " << t4 << std::endl;
  std::cout <<"ref: t5 = " << t5 << std::endl;
  if (t5 > 0) {
    return ((-t3+sqrt(t5))*E);
  }
  else {
    return 10000.0;
  }
}

void G4Incl::forceAbsor(G4int nopart, G4double iarem, G4double izrem, G4double esrem, G4double erecrem,
			G4double alrem, G4double berem, G4double garem, G4double jrem)
{ 
  //  4341        C------------------------------------------------------------------------------
  //   4342             SUBROUTINE FORCE_ABSOR(nopart,F,IAREM,IZREM,ESREM,ERECREM,
  //   4343            s  ALREM,BEREM,GAREM,JREM) 

  //   4346             DIMENSION F(15)
  //   4347             REAL*4 ia1,iz1
  //   4348             
  //   4349             COMMON/hazard/ial,IY(19)
  //   4350       C Dialogue with INCL for nucleus density and parameters.
  //   4351             COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
  //   4352       C RMS espace R, espace P, Fermi momentum and energy for light gauss nuc.      
  //   4353             COMMON/light_gaus_nuc/rms1t(9),pf1t(9),pfln(9),tfln(9),vnuc(9)

  G4int itg = 0;
  G4double sep = 0.0;
  G4double ia1, iz1;
  G4double del = 0.0;
  G4double bimpact;
  G4double bmax, bmaxt;
  G4double proba, proba_trans;
  G4double alea;

  if(nopart != -1) {
    return;
  }

  const G4double fmp = 938.2796;
  bl3->ia2 = calincl->f[0]; // f(1) -> f[0]
  sep = 6.8309;

  if(bl3->ia2 <= 4) {
    if(bl3->ia2 == 2) {
      itg = 6 - 1;
    }
    if(bl3->ia2 == 3 && calincl->f[1] == 1) {
      itg = 7 - 1;
    }
    if(bl3->ia2 == 3 && calincl->f[1] == 2) {
      itg = 8 - 1;
    }
    if(bl3->ia2 == 4) {
      itg = 9 - 1;
    }
    sep = light_gaus_nuc->vnuc[itg] - light_gaus_nuc->tfln[itg]; // :::BUG::: Off-by-one!!!
  }

  if((calincl->f[2] >= 10.0) && (calincl->f[2] <= 100.0)) {
    if(calincl->f[6] == 1.0) {
      bl3->ia1 = 1.0;
      iz1 = 1.0;
      G4double fmpinc = 938.2796;
      G4double pbeam2 = calincl->f[2]*(calincl->f[2] + 2.0*fmpinc);
      bmaxt = bmax;
      proba_trans = coulombTransm(calincl->f[2],ia1,iz1,calincl->f[0],calincl->f[1]);
      proba = forceAbs(1,calincl->f[0],calincl->f[1],calincl->f[2],bmaxt,proba_trans);

      standardRandom(&alea,&(hazard->igraine[4]));
      if(alea > proba) {
        return;
      }

      iarem = calincl->f[0] + bl3->ia1;
      izrem = calincl->f[1] + iz1;
      
      del = sqrt(pow(((calincl->f[0] + 1.0)*fmpinc + calincl->f[2]),2) - pbeam2);
      erecrem = pbeam2/((calincl->f[0] + 1.0)*fmpinc+calincl->f[2] + del);

      esrem = calincl->f[2] + sep - erecrem;

      alrem = 0.00001;
      berem = 0.0;
      garem = 0.99999;
      bimpact = 0.0;
      jrem = 0.0;
      nopart = 0;
      return;
    }
    else if((calincl->f[6] == 2) && (calincl->f[2] >= 20.0)) {
      bl3->ia1 = 1.0;
      iz1 = 0.0;
      G4double fmpinc = 938.2796;
      G4double pbeam2 = calincl->f[2]*(calincl->f[2] + 2.0*fmpinc);
      bmaxt = bmax;
      proba_trans = coulombTransm(calincl->f[2],ia1,iz1,calincl->f[0],calincl->f[1]);
      proba = forceAbs(1,calincl->f[0],calincl->f[1],calincl->f[2],bmaxt,proba_trans);

      standardRandom(&alea,&(hazard->igraine[4]));
      if(alea > proba) {
        return;
      }

      iarem = calincl->f[0] + bl3->ia1;
      izrem = calincl->f[1] + iz1;
      
      del = sqrt(pow(((calincl->f[0]+1.)*fmpinc+calincl->f[2]),2)-pbeam2);
      erecrem = pbeam2/((calincl->f[0] + 1.0)*fmpinc + calincl->f[2] + del);

      esrem = calincl->f[2] + sep - erecrem;

      alrem = 0.00001;
      berem = 0.0;
      garem = 0.99999;
      bimpact = 0.0;
      jrem = 0.0;
      nopart = 0;
      return;
    }
  } //      end if
}

G4double G4Incl::forceAbs(G4double iprojo, G4double at, G4double zt, G4double ep, G4double bmax, G4double pt)
{
  // Results of xabs2 and sig_reac
  G4double sig_exp, sig_incl;
  G4double proba;
  
  G4double ap,zp,A,Z,E; 
  A=at;
  Z=zt;
  E=ep;
  double sig_g=31.41592654*bmax*bmax;
  if(iprojo == 1) {
    ap = 1.0;
    zp = 1.0;
  }
  else {
    ap=1.0;
    zp=0.0;
  }

  sig_exp = xabs2(zp, ap, zt, at, ep);
  sig_incl = crossSection(iprojo, ep, at);

  proba = (sig_exp-pt*sig_incl)/(pt*(sig_g - sig_incl));
  std::cout << "sig_exp : " << sig_exp << " \t \t sig_incl: " << sig_incl << std::endl;
  if(proba <= 0.0) {
    proba = 0.0;
  }
  if(proba > 1.0) {
    proba = 1.0;
  }

  return proba;
}

G4double G4Incl::xabs2(G4double zp, G4double ap, G4double zt, G4double at, G4double ep)
{	                                   
  G4double sig = 0.0;
  
  G4double Const, xzt, xat, Const1, t1, gcm, bcm, plab, ecmp, ecmt, rela, ecm, rm, bigr, bigb;
  G4double xm, x1, sl, phst, ce, term1, delta, beta, twxsec;
  G4double xsec;

  const G4double dp0 = 0.e0, dp1 = 1.e0, dp2 = 2.e0, dp3 = 3.e0, dph = 0.5e0;
  const G4double dp10 = 1.e1, dpth = dp1/dp3, dppi = 3.1415926535898;

  // absoprption xsec revised version rkt-97/5 neutron data from
  // barashenkov this gives absorption xsec for given zp,ap,zt,at,e
  // (mev/nucleon) arguement changed to mev; then e=ep/ap mev/nucleon
  // can be used for neutrons also.  this has coulomb as ours

  G4double E = ep/ap;

  //   nucleon-nucleon inelastc xsec not included here                       
  if ((nint(ap*at) == 1) || (nint(zp+zt) == 1)) {                    
    return dp0;
  }
  G4double rp = radius(ap);                                                     
  G4double rt = radius(at);                                                     
  G4double vp = (dp1 + dpth)*dppi*pow(rp,3);                                          
  G4double vt = (dp1 + dpth)*dppi*pow(rt,3);                                          
  G4double density = dph*((ap/vp) + (at/vt));                                     
  Const=1.75e0*density/8.824728e-02;                                 

  if ((zt < zp) || ((zt == zp) && (at < ap))) {
    xzt = zp;                                                          
    xat = ap;                                                          
    zp = zt;                                                           
    ap = at;                                                           
    zt = xzt;                                                          
    at = xat;                                                          
  }

  if (nint(ap) == 1) {
    Const=2.05;
  }                                   
  if ((nint(zp) == 2) && (nint(ap) == 4)) {
    Const1 = 2.77 - at*8.0e-03 + (at*at)*1.8e-05;                          
  }
  if (nint(zp) == 3) {
    Const=Const/3.0;
  }                                 
  t1=40.0;                                                          
  if (nint(zp) == 0) {
    if ((nint(at) >= 11) && (nint(at) < 40)) {
      t1=30.0;                 
    }
    if (nint(zt) == 14) {
      t1=35.0;
    }                                      
    if (nint(zt) == 26) {
      t1=30.0;
    }                                      
  }

  gcm = (ap*(dp1 + E/938.0) + at)/(pow((pow(ap,2) + pow(at,2) + dp2*ap*(E + 938.0)*at/938.e0),dph));
  bcm = sqrt(dp1-dp1/(pow(gcm,2)));
  plab = ap*sqrt(dp2*938.0*E+E*E);                                    
  ecmp = gcm*(E+938.0)*ap - bcm*gcm*plab - ap*938.0;                     
  ecmt = gcm*938.0*at - at*938.0;
  rela = ecmp + ecmt;                                                    
  ecm = rela;                                                          
  if (ecm < (0.1*rela)) {
    ecm = 0.1*rela;
  }                             
  rm = (197.32/137.01)*zp*zt/ecm;                                      
  bigr = rp + rt + 1.2*(pow(ap,dpth) + pow(at,dpth))/(pow(ecm,dpth));                   
  bigb = 1.44*zp*zt/bigr;                                              
  if ((nint(zp) == 1) && (nint(at) > 56)) {
    bigb = 0.90*bigb;
  }              
  if ((nint(ap) > 56) && (nint(zt) == 1)) {
    bigb = 0.90*bigb;
  }              
  if ((nint(ap) == 1) && (nint(at) == 12)) {
    bigb = 3.5*bigb;
  }               
  if (nint(ap) == 1) {
    if ((nint(at) <= 16) && (nint(at) >= 13)) {
      bigb = (at/7.)*bigb;                                             
    }
    if (nint(zt) == 12) {
      bigb = 1.8*bigb;
    }                               
    if (nint(zt) == 14) {
      bigb = 1.4*bigb;
    }                               
    if (nint(zt) == 20) {
      bigb = 1.3*bigb;
    }                               
  }
  if ((nint(ap) == 1) && (nint(at) < 4)) {
    bigb = 21.0*bigb;
  }               
  if ((nint(ap) < 4) && (nint(at) == 1)) {
    bigb = 21.0*bigb;
  }               
  if ((nint(ap) == 1) && (nint(at) == 4)) {
    bigb = 27.0*bigb;
  }               
  if ((nint(ap) == 4) && (nint(at) == 1)) {
    bigb = 27.0*bigb;
  }               
  if ((nint(zp) == 0) || (nint(zt) == 0)) {
    bigb = dp0;
  }	                     
  xsec = dp10*dppi*bigr*bigr*(dp1-bigb/ecm);			                        
  xm=1.0;                                                             
  if (nint(zp) == 0) {
    if (nint(at) < 200) {
      x1 = 2.83 - 3.1e-02*at + 1.7e-04*at*at;                              
      if (x1 <= 1) {
	x1=1.0;
      }	                                           
      sl=dp1;                                                        
      if (nint(at) == 12) {
	sl=1.6;
      }                                    
      if (nint(at) < 12) {
	sl=0.6;
      }                                    
      xm = (1 - x1*exp(-E/(sl*x1)));                                     
      if (E < 20) {
	//	cout <<"e,xm= " << e << " " << xm << endl;
      }
    }                            
    else {                                                            
      xm = (1-0.3*exp(-(E-1)/15))*(1 - exp(-(E-0.9)));
    }
  }
  if ((nint(zp) == 2) && (nint(ap) == 4)) {
    Const = Const1 - 0.8/(1 + exp((250.0-E)/75.0));                      
  }
  if ((nint(zp) == 1) && (nint(ap) == 1)) {
    if (nint(at) > 45) {
      t1 = 40.0 + at/dp3;                                               
    }
    if (nint(at) < 4) {
      t1 = 55;                                                         
    }
    Const = 2.05 - 0.05/(1 + exp((250.0-E)/75.0));                     
    if (nint(at) < 4) {
      Const = 1.7;
    }                                    
    if (nint(zt) == 12) {
      t1=40.0;                                                     
      Const=2.05 - dp3/(1.0 + safeExp((E - 20.0)/dp10));                    
    }
    if (nint(zt) == 14) {
      t1 = 40.0;                                                     
      Const = 2.05 - 1.75/(1.0 + safeExp((E - 20.0)/dp10));                 
    }                                    
    if (nint(zt) == 18) {
      t1 = 40.0;                                                      
      Const = 2.05 - dp2/(1.0 + safeExp((E - 20.0)/dp10));                    
    }                 
    if (nint(zt) == 20) {
      t1 = 40.0;                                                     
      Const = 2.05 - dp1/(1.0 + safeExp((E - 40.0)/dp10));                    
      Const = Const - 0.25/(1 + exp((250.0 - E)/75.0));                  
    }
    if (nint(zt) >= 35) {
      phst = (nint(zt)-35.e0)/260.e0;	                                 
      Const = Const - phst/(1 + exp((250.0 - E)/75.0));                    
    }
  }

  if ((nint(zp) == 0) && (nint(ap) == 1)) {
    Const = 2*(0.134457/density);                                      
    if ((nint(at) > 140) && (nint(at) <200)) {
      Const = Const - 1.5*(at - dp2*zt)/at;                                                       
    }
    if (nint(at) < 60) {
      Const = Const - 1.5*(at - dp2*zt)/at;
    }            
    if (nint(at) <= 40) {
      Const = Const + 0.25/(dp1 + safeExp(-(170.0 - E)/100.0));                                                         
    }
    if (nint(zt) > 82) {
      Const = Const - zt/(at - zt);                      
    }
    if (nint(zt) >= 82) {
      Const = Const - dp2/(1.0 + safeExp((E - 20.0)/20.0));
    }  
    if ((nint(zt) <= 20) && (nint(zt) >= 10)) {
      Const = Const - dp1/(dp1 + safeExp((E - 20.0)/dp10));
    }                                                
  }
  
  ce = Const * (1.0 - exp(-E/t1)) - 0.292*exp(-E / 792) * cos(0.229*pow(E,0.453));    
  term1 = pow((at*ap),dpth)/(pow(at,dpth) + pow(ap,dpth));                           
  delta = 1.615*term1 - 0.873*ce;                                      
  delta = delta + 0.140*term1/(pow(ecm,dpth));                                 
  delta = delta + 0.794*(at - dp2*zt)*zp/(at*ap);                        
  delta = -delta;                   
  beta = 1.0;                                                           
  twxsec = dp10*dppi*1.26e0*1.26e0*beta*pow((0.873e0*pow(ap,dpth) + 0.873e0*pow(at,dpth)-delta),2);   
  sig = twxsec*(dp1-bigb/ecm)*xm;
  if (sig < dp0) {
    sig = dp0;
  }

  if(isnan(sig)) {
    sig = 0.0;
  }

  return sig;
}

void G4Incl::standardRandom(G4double *rndm, G4int *seed)
{
  // IBM standard random number generator

  (*seed) = (*seed) * 65539;

  if((*seed) < 0) {
    (*seed) = (*seed) + 2147483647+1;
  }

  (*rndm) = (*seed) * 0.4656613e-9;
}

void G4Incl::gaussianRandom(G4double *rndm)
{
  // Gaussian random number generator
  
  G4double tempRandom = 0.0, random = 0.0, randomShuffle = 0.0;
 
  do {
    random = 0.0;
    
    for(G4int i = 0; i < 12; i++) {
      standardRandom(&tempRandom, &(hazard->ial));
      random = random + tempRandom;
    }

    random = random - 6.0;
  } while(random*random > 9);

  // Shuffle the random seeds 
  standardRandom(&randomShuffle, &(hazard->igraine[10]));

  if(randomShuffle > 0.5) {
    standardRandom(&tempRandom, &(hazard->ial));
  }

  (*rndm) = random;
}

G4double G4Incl::safeExp(G4double x)
{                                                
  if (x < -80.0) {
    x = -80.0;                                         
  }
  if (x > 80.0) {
    x = 80.0;
  }                                           

  return exp(x);
}

G4double G4Incl::radius(G4double A)
{                                               
  const G4double dp1 = 1.0, dp3 = 3.0; 
  const G4double dp5 = 5.0, dpth = dp1/dp3;
  const G4int naSize = 23;
  const G4int rmsSize = naSize;

  const G4int na[naSize] = {1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26};
  
  const G4double rms[rmsSize] = {0.85,2.095,1.976,1.671,2.57,2.41,2.519,2.45,
				 2.42,2.471,2.440,2.58,2.611,2.730,2.662,2.727,
				 2.900,3.040,2.969,2.94,3.075,3.11,3.06};
             
  G4double fact = sqrt(dp5/dp3);                                                
  G4int ia = int(floor(A+0.4));                                                        
  G4double result = fact * (0.84 * pow(A,dpth) + 0.55);                               
  for(G4int i = 0; i < naSize; i++) {
    if (ia == na[i]) {
      result = fact*rms[i];
    }
  }
  
  return result;                                                            
}

G4double G4Incl::crossSection(G4int projectile, G4double E, G4double A) 
{
  const G4double coefp[4][3] = {{ -5.9260e-9, 6.89450e-6, -6.0980e-6},
				{  2.1544e-6,-1.84800e-3, -5.9820e-4},
				{ -2.5900e-4, 0.17595e+0,  1.1741e+0},
				{  1.1504e-3, 2.82810e+0,-28.7300e+0}};
 
  const G4double coefn[4][3] = {{1.6105e-9,3.3985e-7,1.4678e-5},
				{-5.35e-7,-3.465e-4,-0.01633e+0},
				{-4.755e-6,0.07608e+0,2.8135e+0},
				{-3.622e-3,3.5924e+0,-38.294e+0}};
 
  const G4double coef2p[5][3] = {{6.8108e-9,-2.163e-7,2.1898e-6},
				 {-2.187e-6,7.8331e-5,-7.164e-4},
				 {2.3651e-4,-9.690e-3,0.076424e+0},
				 {-9.195e-3,0.5030e+0,-2.4979e+0},
				 {-0.01087e+0,2.6494e+0,-2.5173e+0}}; 

  G4double apow[3], epow[5];
  G4int ii, jj;

  if(A >= 27.0) {
    ii = 3;        
    jj = 4;
  }
  else {
    ii = 3;
    jj = 5;
  }

  for(int i = 1; i <= ii; i++) {
    apow[i-1] = pow(A,(ii-i));
  }

  for(int j = 1; j <= jj; j++) {
    epow[j-1] = pow(E,(jj-j));
  }

  double result = 0.0;

  if(A >= 27.0) { 
    if(projectile == 1) {
      for(G4int i = 0; i < ii; i++) {
	for(G4int j = 0; j < jj; j++) {
	  result = result + coefp[j][i]*apow[i]*epow[j];
	}
      }
    }
    else {
      for(G4int i = 0; i < ii; i++ ) {
	for(G4int j = 0; j < jj; j++) {
	  result = result + coefn[j][i]*apow[i]*epow[j];
	}
      }
    }
  }
  else {
    for(G4int i = 0; i < ii; i++) {
      for(G4int j = 0; j < jj; j++) {
	result = result + coef2p[j][i]*apow[i]*epow[j];
      }
    }
  }

  return result;
}

G4double G4Incl::coulombTransm(G4double E, G4double fm1, G4double z1, G4double fm2, G4double z2)
{
  G4double eta,rho;
  const G4double c2 = 0.00516;
  const G4double c3 = 0.007165;
  const G4double uma = 938.0;
  G4double ml;
  
  G4double ecm = E*fm2/(fm1+fm2);
  G4double fm = fm1*fm2*uma/(fm1+fm2);
  G4double r1 = 0.0;
  G4double r2 = 0.0;
  if(fm1 >= 2.0) {
    r1 = 1.2*pow(fm1,0.33333333);
  }
  r2 = 1.2*pow(fm2,0.33333333);
  eta = c2*z1*z2*sqrt(fm/ecm);
  rho = c3*(r1+r2)*sqrt(fm*ecm);

  return clmb1(rho,eta,&ml);
}

G4double G4Incl::clmb1(G4double rho, G4double eta, G4double *ml) 
{                                        
  const G4double dp1 = 1.0, dp2 = 2.e0, dp4 = 4.e0, dph = 0.5, dp5 = 5.0; 

  const G4double prm1 = 69.06;                                          
  const G4int ln0 = 81, lt0 = 21;                                        
  const G4int ln1 = 61, lt1 = 61;                                        
  static G4double psi0[ln0], trans0[lt0][ln0], x0[lt0], f0[ln0], psi1[ln1];
  static G4double trans1[lt1][ln1], x1[lt1], f1[ln1];                                

  const G4double pi = 3.14159;                                               
  const G4double c0 = 0.11225, c1 = dph, gamma = 0.5772157e0, s3 = 0.2020569;
  const G4double s4 = 0.08232323;
  static G4double y = dp2*eta;                                                         
  static G4double psi = rho*y;                                                         

  static G4int i0, j0;
  
  static G4double prob;
  static G4double dumm, x, cx;
  static G4double t, t1, t2, t3;
  static G4double f, g;
  static G4double temp, temp1, temp2;
  static G4double xk, delp0, delp1, delx0, delx1;

  if (rho > y) {                                                 
    if (psi > dp4 && psi < 50.0) {                           
      prob = clmb2(rho,eta,&dumm);                                       
    }
    else {
      x = exp(log(eta)/6.0);
      prob = sqrt(dp1 - y*x/(c0 + c1 * pow(x,3) + rho * x));
    } 
    (*ml) = 0;
  }                                                            
  else {
    x = rho/y;                                                          
    if (psi <= psi0[0]) {
      t = min(pi*y,prm1);                                              
      cx = t/(exp(t) - dp1);                                             
      t1 = cos(rho) * (dp1-0.75*pow(psi,2) + dp5 * x * pow(psi,2)) - dph * psi * rho * sin(rho);
      t2 = dp1 + dph * psi * (dp1 - x/6.0);                                   
      if (eta > dp1) {                                          
	t3 = log(psi)+dp2*gamma - dp1 + dp1/(12.e0*pow(eta,2))+dp1/(12.e1*pow(eta,4));
      }                                                        
      else {                                                          
	t3 = log(dp2*rho) + gamma - dp1/(dp1 + pow(eta,2)) + s3*pow(eta,2) + s4*pow(eta,4);   
      }
      g = t1 + psi*t2*t3;                                                
      f = cx*rho*t2;                                                    
      prob = cx/(pow(g,2)+pow(f,2));                                           
      (*ml) = 3;
    }                                                          
    else if (psi <= psi0[ln0-1]) {
      if (x <= x0[0]) {
	temp = log(psi/psi0[0]);
	j0 = 1 + int(temp/delp0);                                        
	j0 = min(max(j0,1),(ln0-1));                                     
	temp = temp - (j0-1)*delp0;                                      
	t = f0[j0] + (f0[j0+1] - f0[j0])*temp/delp0;                       
	xk = x*sqrt(psi);                                              
	prob = (dp1+3.33e-1*x+3.e-1*xk+1.e-1*pow(xk,2))*exp(t);            
	t = min(pi*y,prm1);                                            
	cx = t/(exp(t)-dp1);                                           
	prob = cx/pow(prob,2);                                             
	(*ml) = 1;
      }                                                        
      else {                                                          
	temp1 = log(x/x0[0]);
	i0 = min(max(1 + int(temp1/delx0),1),lt0-1);                                     
	temp1 = temp1 - (i0 - 1)*delx0;                                    
	temp2 = log(psi/psi0[0]);
	j0 = min(max(1+int(temp2/delp0),1),ln0-1);                                     
	temp2 = temp2-(j0-1)*delp0;                                    
	t1 = trans0[i0][j0] + (trans0[i0+1][j0] - trans0[i0][j0]) * temp1/delx0;
	t2 = trans0[i0][j0+1] + (trans0[i0+1][j0+1] - trans0[i0][j0+1]) * temp1/delx0;
	prob = exp(t1 + (t2 - t1)*temp2/delp0);                                   
	(*ml)=2;                                                        
      }
    }
    else if (psi <= psi1[ln1-1]) {
      if (x <= x1[0]) {
	temp = log(psi/psi1[0]);
	j0 = min(max(1+int(temp/delp1), 1), ln1-1);
	t = f1[j0]+(f1[j0+1]-f1[j0])*(temp - (j0 - 1)*delp1)/delp1;                       
	xk = x*sqrt(psi);                                              
	prob = (dp1+3.33e-1*x+3.0e-1*xk+1.e-1*pow(xk,2))*exp(t);            
	t = min(pi*y,prm1);                                            
	cx = t/(exp(t)-dp1);                                           
	prob = cx/pow(prob,2);                                             
	(*ml) = 1;
      }                                                        
      else {                                                         
	temp1 = log(x/x1[0]);
	i0 = min(max(1+int(temp1/delx1),1),lt1-1);                                    
	temp1 = temp1-(i0-1)*delx1;                                    
	temp2 = log(psi/psi1[0]);
	j0 = min(max(1+int(temp2/delp1),1),ln1-1);                                     
	temp2 = temp2 - (j0-1)*delp1;                                    
	t1 = trans1[i0][j0] + (trans1[i0+1][j0] - trans1[i0][j0])*temp1/delx1;
	t2 = trans1[i0][j0+1] + (trans1[i0+1][j0+1] - trans1[i0][j0+1])*temp1/delx1;                                                      
	prob = exp(t1 + (t2-t1)*temp2/delp1);                                    
	(*ml)=2;                                                        
      }
    }
    else {
      prob = clmb2(rho,eta,&dumm);                                       
      (*ml) = 4;                                              
    }
  }

  return prob;
}

G4double G4Incl::clmb2(G4double rho, G4double eta, G4double *t1)
{
  const G4double dp0 = 0.0, dp1 = 1.0, dp2 = 2.0, dp3 = 3.0;
  const G4double dph = 0.5, dpth = dp1/dp3;                                                   
  const G4int ln = 102;

  const G4double t0[ln] = {0.0, dp0,.1083,.1369,.1572,.1736,.1876,.2,.2113,
			   0.2216,.2312,.2403,.2489,.2571,.265,.2725,.2798,
			   .2869,.2938,.3006,.3071,.3136,.3199,.3261,.3322,
			   .3382,.3442,.3499,.3557,.3615,.3672,.3729,.3785,
			   .3841,.3897,.3952,.4008,.4063,.4118,.4173,.4228,
			   .4283,.4338,.4393,.4448,.4504,.4559,.4615,.4671,
			   .4728,.4784,.4841,.4899,.4957,.5015,.5074,.5133,
			   .5193,.5253,.5315,.5376,.5439,.5503,.5567,.5632,
			   .5698,.5765,.5833,.5903,.5973,.6045,.6118,.6193,
			   .6269,.6346,.6426,.6507,.659,.6675,.6763,.6853,
			   .6945,.704,.7139,.724,.7345,.7453,.7566,.7683,
			   .7805,.7932,.8065,.8205,.8352,.8508,.8673,.8849,
			   .9038,.9243,.9467,.9715,dp1};                                            
  const G4double x1 = 0.01;
  const G4double xi = 100;                                        
  static G4double x,temp,prob;
  static G4int i;

  x = dp1/(dp1 + sqrt(dph*rho*eta));                                     
  if (x < x1) {                                                 
    temp = t0[2] * pow((x/x1),dpth);
  }                                       
  else {
    i = int(floor(xi*x));                                                          
    i = i + 1;
    if(i == 101) {
      i = 100;
    }
    i = max(min(i,ln-1),2); // 2->1
    temp = t0[i] + (t0[i+1] - t0[i]) * (x - i * x1)/x1;                          
  }
  (*t1) = dp1 - temp;                                                       
  prob = dp1 - dp2 * (*t1) * eta/rho;                                            

  return (max(prob,dp0));
}

// Utilities

G4double G4Incl::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Incl::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4Incl::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Incl::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Incl::nint(G4double number)
{
  G4double intpart;
  G4double fractpart;
  fractpart = modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return int(floor(number));
    }
    else {
      return int(ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return int(floor(number));
    }
    else {
      return int(ceil(number));
    }
  }
  
  return 0;
}

G4double G4Incl::callFunction(G4int functionChoice, G4double r)
{
  if(functionChoice == wsaxFunction) {
    return wsax(r);
  }
  else if(functionChoice == derivWsaxFunction) {
    return derivWsax(r);
  }
  else if(functionChoice == dmhoFunction) {
    return dmho(r);
  }
  else if(functionChoice == derivMhoFunction) {
    return derivMho(r);
  }
  else if(functionChoice  == derivGausFunction) {
    return derivGaus(r);
  }
  else if(functionChoice == densFunction) {
    return dens(r);
  }

  return 0.0;
}

G4double G4Incl::am(G4double a, G4double b, G4double c, G4double d) 
{
  return sqrt(d*d-a*a-b*b-c*c);
}

G4double G4Incl::pcm(G4double e, G4double a, G4double c)
{
  return (0.5*sqrt((pow(e,2)-pow((a+c),2))*(pow(e,2)-pow((a-c),2)))/e);
}

G4double G4Incl::sign(G4double a, G4double b)
{
  if(b >= 0) {
    return utilabs(a);
  }
  if(b < 0) {
    return (-utilabs(a));
  }

  std::cout <<"Error: sign function failed. " << std::endl;
  return a; // The algorithm is never supposed to reach this point.
}

G4double G4Incl::utilabs(G4double a)
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

  std::cout <<"Error: utilabs function failed. " << std::endl;

  return a;
}

G4double G4Incl::amax1(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else if(b > a) {
    return b;
  }
  else if(a == b) {
    return a;
  }

  return a; // The algorithm is never supposed to reach this point.
}

G4double G4Incl::w(G4double a, G4double b, G4double c, G4double d)
{
  return (sqrt(a*a+b*b+c*c+d*d));
}

G4int G4Incl::idnint(G4double a)
{
  G4int value = 0;

  G4int valueCeil = int(ceil(a));
  G4int valueFloor = int(floor(a));

  if(abs(value - valueCeil) < abs(value - valueFloor)) {
    return valueCeil;
  }
  else {
    return valueFloor;
  }
}
