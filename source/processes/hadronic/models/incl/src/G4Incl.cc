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
// $Id: G4Incl.cc,v 1.37 2010-12-15 07:41:31 gunter Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4Incl.hh"
#include <iostream>
#include "Randomize.hh"
#include "G4InclRandomNumbers.hh"
#include "G4Ranecu.hh"

G4Incl::G4Incl()
{
  verboseLevel = 0;

  // Set functions to be used for integration routine.  
  wsaxFunction = 0;
  derivWsaxFunction = 1;
  dmhoFunction = 2;
  derivMhoFunction = 3;
  derivGausFunction = 4;
  densFunction = 5;

  // Default: no Fermi break-up
  useFermiBreakup = false;

  // Default: no projectile spectators
  useProjSpect = false;

  randomGenerator = new G4InclGeant4Random();
  //  randomGenerator = new G4Ranecu();
  clearState();
}

G4Incl::G4Incl(G4Hazard *aHazard, G4Dton *aDton, G4Saxw *aSaxw, G4Ws *aWs)
{
  verboseLevel = 0;

  // Default: no Fermi break-up
  useFermiBreakup = false;

  // Default: no projectile spectators
  useProjSpect = false;

  // Set functions to be used for integration routine.  
  wsaxFunction = 0;
  derivWsaxFunction = 1;
  dmhoFunction = 2;
  derivMhoFunction = 3;
  derivGausFunction = 4;
  densFunction = 5;

  // Set input data for testing.
  hazard = aHazard;
  dton = aDton;
  saxw = aSaxw;
  ws = aWs;

  //  randomGenerator = new G4Ranecu();
  randomGenerator = new G4InclGeant4Random();
  clearState();
}

G4Incl::G4Incl(G4Hazard *aHazard, G4InclInput *aCalincl, G4Ws *aWs, G4Mat *aMat, G4VarNtp *aVarntp)
{
  verboseLevel = 0;

  // Default: no Fermi break-up
  useFermiBreakup = false;

  // Default: no projectile spectators
  useProjSpect = false;

  // Set functions to be used for integration routine.    
  wsaxFunction = 0;
  derivWsaxFunction = 1;
  dmhoFunction = 2;
  derivMhoFunction = 3;
  derivGausFunction = 4;
  densFunction = 5;

  // Set input data for INCL run.  
  hazard = aHazard;
  calincl = aCalincl;
  ws = aWs;
  mat = aMat;
  varntp = aVarntp;

  randomGenerator = new G4InclGeant4Random();
  //  randomGenerator = new G4Ranecu();
  light_gaus_nuc = new G4LightGausNuc();
  light_nuc = new G4LightNuc();
  spl2 = new G4Spl2();
  saxw = new G4Saxw();
  dton = new G4Dton();
  bl1 = new G4Bl1();
  bl2 = new G4Bl2();
  bl3 = new G4Bl3();
  bl4 = new G4Bl4();
  bl5 = new G4Bl5();
  bl6 = new G4Bl6();
  bl8 = new G4Bl8();
  bl9 = new G4Bl9();
  bl10 = new G4Bl10();
  kindstruct = new G4Kind();
  bev = new G4Bev();
  paul = new G4Paul();
  varavat = new G4VarAvat();
  varavat->kveux = 0;

  be = new G4VBe();
  ps = new G4InclProjSpect();
  fermi = new G4InclFermi();
  qvp = new G4QuadvectProjo();

  volant = new G4Volant();
  volant->iv = 0;
  evaporationResult = new G4VarNtp();
  evaporationResult->ntrack = -1;

  // Initialize evaporation.
  abla = new G4Abla(hazard, volant, evaporationResult);
  abla->initEvapora();

  theLogger = 0;
  clearState();
}

G4Incl::~G4Incl()
{
  delete be;
  delete ps;
  delete fermi;
  delete qvp;

  delete randomGenerator;
  delete light_gaus_nuc;
  delete light_nuc;
  delete spl2;
  delete saxw;
  delete dton;
  delete bl1;
  delete bl2;
  delete bl3;
  delete bl4;
  delete bl5;
  delete bl6;
  delete bl8;
  delete bl9;
  delete bl10;
  delete kindstruct;
  delete bev;
  delete paul;
  delete varavat;

  delete abla;
  delete evaporationResult;
  delete volant;
}

void G4Incl::setUseFermiBreakUp(G4bool useIt)
{
  useFermiBreakup = useIt;
}

void G4Incl::setUseProjectileSpectators(G4bool useIt)
{
  useProjSpect = useIt;
}

/**
 *Methods for debugging.
 */

void G4Incl::dumpParticles()
{
  G4int ia = bl3->ia1 + bl3->ia2;
  G4cout <<"Nucleons: (number of nucleons = " << ia << ")" << G4endl;
  for(G4int i = 0; i <= ia; i++) {
    G4cout <<"x1(" << i << ") = " << bl3->x1[i] << G4endl;
    G4cout <<"x2(" << i << ") = " << bl3->x2[i] << G4endl;
    G4cout <<"x3(" << i << ") = " << bl3->x3[i] << G4endl;
    G4cout <<"p1(" << i << ") = " << bl1->p1[i] << G4endl;
    G4cout <<"p2(" << i << ") = " << bl1->p2[i] << G4endl;
    G4cout <<"p3(" << i << ") = " << bl1->p3[i] << G4endl;
    G4cout <<"eps(" << i << ") = " << bl1->eps[i] << G4endl;
  }
}

G4double G4Incl::energyTest(G4int i)
{
  return am(bl1->p1[i]+bl1->p1[i],bl1->p2[i]+bl1->p2[i],bl1->p3[i]+bl1->p3[i],bl1->eps[i]+bl1->eps[i]);
}

void G4Incl::dumpBl5(std::ofstream& dumpOut)
{
  dumpOut <<"Dumping G4Bl5:" << G4endl;
  for(G4int i = 0; i < 300; i++) {
    dumpOut <<"i = " << i << " nesc[i] = " << bl5->nesc[i] << G4endl;
  }
}

void G4Incl::dumpSaxw(std::ofstream& dumpOut)
{
  dumpOut << "Dumping G4Saxw" << G4endl;
  dumpOut << "saxw->k    = " << saxw->k    << G4endl;
  dumpOut << "saxw->n    = " << saxw->n    << G4endl;
  dumpOut << "saxw->imat = " << saxw->imat << G4endl;
  for(G4int i = 0; i < 30; i++) {
    dumpOut <<"i = " << i << " x = " << saxw->x[i][0] << " y = " << saxw->y[i][0] << " s = " << saxw->s[i][0] << G4endl; 
  }
}

void G4Incl::dumpBl1(std::ofstream& dumpOut)
{
  dumpOut <<"Dumping Bl1: " << G4endl;
  dumpOut <<"bl1->ta = " << bl1->ta << G4endl;
  for(G4int i = 0; i <= bl2->k; i++) {
    dumpOut <<"i = " << i;
    dumpOut <<" bl1->p1 = "   << bl1->p1[i]    << " bl1->p2 = "  << bl1->p2[i] <<" bl1->p3 = " << bl1->p3[i];
    dumpOut <<" bl1->eps = "  << bl1->eps[i];
    dumpOut <<" bl1->ind1 = " << bl1->ind1[i]  <<" bl1->ind2 = " << bl1->ind2[i] << G4endl;
  }
  dumpOut <<"End of Bl1 dump." << G4endl << G4endl;
}

void G4Incl::dumpBl2(std::ofstream& dumpOut)
{
  dumpOut <<"Dumping Bl2: bl2->k = " << bl2->k << G4endl;
  for(G4int i = 0; i <= bl2->k; i++) {
    dumpOut <<"i = " << i;
    dumpOut <<" bl2->ind " << bl2->ind[i] << " bl2->jnd = " << bl2->jnd[i] <<" bl2->crois = " << bl2->crois[i] << G4endl;
  }
  dumpOut <<"End of Bl2 dump." << G4endl << G4endl;}


void G4Incl::dumpBl3(std::ofstream& dumpOut)
{
  dumpOut <<"Dumping Bl3:" << G4endl;
  dumpOut <<"r1   = " << bl3->r1   << " r2 = "  << bl3->r2  << G4endl;
  dumpOut <<"ia1  = " << bl3->ia1  << " ia2 = " << bl3->ia2 << G4endl;
  dumpOut <<"rab2 = " << bl3->rab2                          << G4endl;

  for(G4int i = 0; i <= bl2->k; i++) {
    dumpOut <<"i = "        << i;
    dumpOut <<" bl3->x1 = " << bl3->x1[i] << " bl3->x2 = " << bl3->x2[i] <<" bl3->x3 = " << bl3->x3[i] << G4endl;
  }
  dumpOut <<"End of Bl2 dump." << G4endl << G4endl;
}

// End debug functions

void G4Incl::setVerboseLevel(G4int level)
{
  verboseLevel = level;
  if(verboseLevel > G4InclUtils::silent) {
    G4cout <<";; G4Incl: Setting verbose level to " << verboseLevel << G4endl;
  }
  abla->setVerboseLevel(level);
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

void G4Incl::setInput(G4InclInput *newCalincl)
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

/**
 * INCL main routines for event processing.
 */

void G4Incl::processEventIncl(G4InclInput *input)
{
  //G4cout <<"Starting event " << eventnumber << G4endl;
  if(input == 0) {
    G4cerr <<"G4Incl fatal error: NULL pointer passed as input!" << G4endl;
    return;
  }
  calincl = input;

  const G4double uma = 931.4942;
  const G4double melec = 0.511;
  const G4double fmp = 938.2796;

  G4double pcorem = 0.0;
  G4double pxrem = 0.0;
  G4double pyrem = 0.0;
  G4double pzrem = 0.0;

  G4double ap = 0.0, zp = 0.0, mprojo = 0.0, pbeam = 0.0;

  varntp->clear();

  if(calincl->bulletType() == -12) {
    be->ia_be = std::abs(calincl->bulletType());
    be->iz_be = 6;
  } else if(calincl->bulletType() == -666) {
    be->iz_be = calincl->extendedProjectileZ();
    be->ia_be = calincl->extendedProjectileA();
  }

  if(calincl->isExtendedProjectile() == false && calincl->bulletType() < -max_a_proj) {
  //  if(calincl->bulletType() < -max_a_proj) {
    G4cout <<"max a of composite projectile is: " << max_a_proj << G4endl;
    exit(0);
  }
  if(calincl->bulletType() < 0) {
    //    calincl->bulletType() = std::floor(calincl->bulletType() + 0.1); WTF???
    be->pms_be=100.;
    G4int i_tabled=0;
    if(be->iz_be == 3 && be->ia_be == 6) {
      be->rms_be=2.56;
      be->bind_be=32.0;
      i_tabled=1;
    } else if(be->iz_be == 3 && be->ia_be == 7) { // TODO: Check the values!
      be->rms_be=2.56;
      be->bind_be=32.0;
      i_tabled=1;
    } else if(be->iz_be == 3 && be->ia_be == 8) {
      be->rms_be=2.40;
      be->bind_be=39.25;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 7) {
      be->rms_be=2.51;
      be->bind_be=58.17;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 9) {
      be->rms_be=2.51;
      be->bind_be=58.17;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 10) {
      be->rms_be=2.45;
      be->bind_be=64.75;
      i_tabled=1;
    } else if(be->iz_be == 5 && be->ia_be == 10) {
      be->rms_be=2.45;
      be->bind_be=64.75;
      i_tabled=1;
    } else if(be->iz_be == 5 && be->ia_be == 11) {
      be->rms_be=2.40;
      be->bind_be=76.21;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 9) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 10) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 11) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && calincl->bulletType() == -12) { // Special Carbon case
      G4cout <<"Carbon 12 (special) selected." << G4endl;
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 12) {
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 7 && be->ia_be == 16) {
      be->rms_be=2.73;
      be->bind_be=127.62;
      i_tabled=1;
    } else {
      G4cout <<"Warning: No rms and binding for projectile ion A = " << be->ia_be << " Z = " << be->iz_be << G4endl;
      be->rms_be=2.44;
      be->bind_be=92.17;
      G4cout <<"Warning: Using probably bad values rms = " << be->rms_be << " binding = " << be->bind_be << G4endl;
      i_tabled=1;      
    }
      
    if(i_tabled == 0) {
      G4cout <<"This heavy ion (a,z)= " << be->ia_be << " " << be->iz_be << " is not defined as beam in INCL" << G4endl;
      exit(0);
    }
      
    //    G4cout <<"z projectile, rms_r, rms_p (gaussian model)" << be->iz_be << " " << be->rms_be << " " << be->pms_be << G4endl;
    //    G4cout <<"binding energy (mev):" << be->bind_be << G4endl;
    //    G4cout <<"fermi-breakup dresner below a=" << calincl->f[11] << G4endl;
  }      
  //  G4cout <<"Target Mass and Charge: " << calincl->targetA() << " " << calincl->targetZ() << G4endl;
  //  calincl->f[10] = 0; // No clusters

  if(calincl->bulletType() == -12) {  // C12 special case
    mprojo=fmp*std::abs(calincl->bulletType()) - be->bind_be;
    pbeam=std::sqrt(calincl->bulletE()*(calincl->bulletE()+2.*mprojo));
    ap=std::abs(calincl->bulletType());
    zp=be->iz_be;
  } else if(calincl->bulletType() == -666) { // Generic extended projectile
    mprojo=fmp*be->ia_be - be->bind_be;
    pbeam=std::sqrt(calincl->bulletE()*(calincl->bulletE()+2.*mprojo));
    ap=be->ia_be;
    zp=be->iz_be;
  }
  // pi+
  if(calincl->bulletType() == 3) { 
    mprojo = 139.56995;
    ap = 0.0;
    zp = 1.0;
  }

  // pi0
  if(calincl->bulletType() == 4) {
    mprojo = 134.9764;
    ap = 0.0;
    zp = 0.0;
  }

  // pi-
  if(calincl->bulletType() == 5) {
    mprojo = 139.56995;
    ap = 0.0;
    zp = -1.0;
  }

  // coulomb en entree seulement pour les particules ci-dessous

  // proton
  if(calincl->bulletType() == 1) {
    mprojo = 938.27231;
    ap = 1.0;
    zp = 1.0;
  }

  // neutron  
  if(calincl->bulletType() == 2) {
    mprojo = 939.56563;
    ap = 1.0;
    zp = 0.0;
  }

  // deuteron
  if(calincl->bulletType() == 6) {
    mprojo = 1875.61276;
    ap = 2.0;
    zp = 1.0;
  }

  // triton
  if(calincl->bulletType() == 7) {
    mprojo = 2808.95;
    ap = 3.0;
    zp = 1.0;
  }

  // He3
  if(calincl->bulletType() == 8) {
    mprojo = 2808.42;
    ap = 3.0;
    zp = 2.0;
  }

  // Alpha
  if(calincl->bulletType() == 9) {
    mprojo = 3727.42;
    ap = 4.0;
    zp = 2.0;
  }

  // Carbon
  if(calincl->bulletType() == -12) {
    mprojo = 6.0*938.27231 + 6.0*939.56563;
    ap = 12.0;
    zp = 6.0;
  }

  pbeam = std::sqrt(calincl->bulletE()*(calincl->bulletE() + 2.0*mprojo));         

  G4double at = calincl->targetA();
       
  G4int ibert = 1;

  G4int nopart = 0;
  G4int izrem = 0;
  G4int iarem = 0;
  G4double esrem = 0.0;
  G4double erecrem = 0.0;
  G4double berem = 0.0;
  G4double garem = 0.0;
  G4double bimpac = 0.0;
  G4int jrem = 0;
  G4double xjrem = 0.0, yjrem = 0.0, zjrem = 0.0;
  G4double alrem = 0.0;

  // Coulomb barrier
  
  G4double probaTrans = 0.0;
  G4double rndm = 0.0;

  if((calincl->bulletType() == 1) || (calincl->bulletType() >= 6)) {
    //    probaTrans = coulombTransm(calincl->bulletE(),apro,zpro,calincl->targetA(),calincl->targetZ());
    probaTrans = coulombTransm(calincl->bulletE(),ap,zp,calincl->targetA(),calincl->targetZ());
    standardRandom(&rndm, &(hazard->ial));
    if(rndm <= (1.0 - probaTrans)) {
      varntp->ntrack = -1;
      return;
    }
  }

  //  G4cout <<"Before PNU:" << G4endl;
  //  randomGenerator->printSeeds();
  // Call the actual INCL routine:
  pnu(&ibert, &nopart,&izrem,&iarem,&esrem,&erecrem,&alrem,&berem,&garem,&bimpac,
      &jrem, &xjrem, &yjrem, &zjrem);
  //  G4cout <<"After PNU:" << G4endl;
  //  randomGenerator->printSeeds();
  G4double mrem = int(zjrem/197.328); // CHECK
  if (mrem > jrem) mrem=jrem;
  if (mrem < -jrem) mrem=-jrem;

//   nopart=1;
//   kind[0]=1;
//   ep[0]=799.835;
//   alpha[0]=0.08716;
//   beta[0]=0.;
//   gam[0]=0.99619;
//   izrem=82;
//   iarem=208;
//   esrem=200.;
//   erecrem=0.18870;
//   alrem=-0.47101;
//   berem=0.;
//   garem=0.88213;
//   bimpac=2.;
  forceAbsor(&nopart, &iarem, &izrem, &esrem, &erecrem, &alrem, &berem, &garem, &jrem);
  G4double aprf = double(iarem);    // mass number of the prefragment
  G4double jprf = 0.0;                // angular momentum of the prefragment

  // Mean angular momentum of prefragment                                  
  jprf = 0.165 * std::pow(at,(2.0/3.0)) * aprf*(at - aprf)/(at - 1.0);                               
  if (jprf < 0) {
    jprf = 0.0;
  }

  // check m.de jong, ignatyuk, schmidt nuc.phys a 613, pg442, 7th line
  jprf = std::sqrt(2*jprf);

  jprf = jrem;
  varntp->jremn = jrem;      // jrem copie dans le ntuple

  G4double numpi = 0;  // compteurs de pions, neutrons protons
  G4double multn = 0; 
  G4double multp = 0;

  // ecriture dans le ntuple des particules de cascade (sauf remnant)      
  varntp->ntrack = nopart;          // nombre de particules pour ce tir
  if(varntp->ntrack >= VARNTPSIZE) {
    if(verboseLevel > 2) {
      G4cout <<"G4Incl error: Output data structure not big enough." << G4endl;
    }
  }
  varntp->massini = iarem;
  varntp->mzini = izrem;
  varntp->exini = esrem;
  varntp->bimpact = bimpac;
  
  //  three ways to compute the mass of the remnant: 
  //                -from the output of the cascade and the canonic mass
  //                -from energy balance (input - all emitted energies)
  //                -following the approximations of the cugnon code (esrem...)
  G4double mcorem = mprojo + calincl->bulletE() + abla->pace2(double(calincl->targetA()),double(calincl->targetZ()))
    + calincl->targetA()*uma - calincl->targetZ()*melec;

  G4double pxbil = 0.0;
  G4double pybil = 0.0;
  G4double pzbil = 0.0;         

  if(nopart > -1) {
    // Fill the projectile spectator variables
    varntp->masp = ps->a_projspec;
    varntp->mzsp = ps->z_projspec;
    varntp->exsp = ps->ex_projspec;
    varntp->spectatorP1 = ps->p1_projspec;
    varntp->spectatorP2 = ps->p2_projspec;
    varntp->spectatorP3 = ps->p3_projspec;
    varntp->spectatorT = ps->t_projspec;
    for(G4int j = 0; j <= nopart; j++) {
      if(ep[j] < 0.0) { // Workaround to avoid negative energies (and taking std::sqrt of a negative number).
	G4cout <<"G4Incl: Not registering particle with energy: " << ep[j] << G4endl;
	continue;
      }
      if(kind[j] == 0) continue; // Empty particle rows are sometimes produced by lurking indexing problems. We can simply skip those "bad" entries...
      if(gam[j] > CLHEP::pi) {
	if(verboseLevel > 2) {
	  G4cout <<"G4Incl: Just avoided floating point exception by using an ugly hack..." << G4endl;
	}
	continue; // Avoid floating point exception
      }

      varntp->itypcasc[j] = 1; // Particle was produced by the cascade
      // Spectators of composite projectiles (7/2006, AB)
      // (kind is negative in that case)
      if(kind[j] <= 0) { // Particle is a projectile spectator that comes directly from the cascade
	kind[j] *= -1;
	varntp->itypcasc[j]=-1;
	//	G4cout <<"Spectator registered!" << G4endl;
	//	continue;
      }
	
      // kind(): 1=proton, 2=neutron, 3=pi+, 4=pi0, 5=pi -      
      if(kind[j] == 1) { 
	varntp->avv[j] = 1;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+1876.5592)); // cugnon
	multp = multp + 1;
	mcorem = mcorem - ep[j] - 938.27231;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Proton produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 2) { 
	varntp->avv[j] = 1;
	varntp->zvv[j] = 0;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+1876.5592)); // cugnon
	multn = multn + 1;
	mcorem = mcorem - ep[j] - 939.56563;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Neutron produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }
      
      if(kind[j] == 3) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 139.56995;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi+ produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 4) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = 0;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 134.9764;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi0 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 5) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = -1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 139.56995;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi+ produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 6) { 
	varntp->avv[j] = 2;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*1874.34)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2806.359;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Deuteron produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 7) { 
	varntp->avv[j] = 3;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*2806.359)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2806.359;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Triton produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 8) { 
	varntp->avv[j] = 3;
	varntp->zvv[j] = 2;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*2807.119)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2807.119;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: He3 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 9) { 
	varntp->avv[j] = 4;
	varntp->zvv[j] = 2;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*3724.818)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 3724.818;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: He3 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      varntp->enerj[j] = ep[j];
      varntp->tetlab[j] = 180.0*std::acos(gam[j])/3.141592654;
      varntp->philab[j] = 180.0*std::atan2(beta[j],alpha[j])/3.141592654;
      pxbil = pxbil + varntp->plab[j]*alpha[j];
      pybil = pybil + varntp->plab[j]*beta[j];
      pzbil = pzbil + varntp->plab[j]*gam[j];

      if(verboseLevel > 3) {
	G4cout <<"Momentum: " << varntp->plab[j]   << G4endl;
	G4cout <<"Theta: "    << varntp->tetlab[j] << G4endl;
	G4cout <<"Phi: "      << varntp->philab[j] << G4endl;
      }
    }

    // calcul de la masse (impulsion) du remnant coherente avec la conservation d'energie:
    pcorem = std::sqrt(erecrem*(erecrem + 2.0 * 938.2796 * iarem));   // cugnon
    mcorem = 938.2796 * iarem;                               // cugnon
    varntp->pcorem = pcorem;
    varntp->mcorem = mcorem;
    // Note: Il faut negliger l'energie d'excitation (ESREM) pour que le bilan 
    // d'impulsion soit correct a la sortie de la cascade.....et prendre la
    // masse MCOREM comme ci-dessus (fausse de ~1GeV par rapport aux tables...)        
    pxrem = pcorem * alrem;
    pyrem = pcorem * berem;
    pzrem = pcorem * garem;
    varntp->pxrem = pxrem;
    varntp->pyrem = pyrem;
    varntp->pzrem = pzrem;
    pxbil = pxbil + pxrem;
    pybil = pybil + pyrem;
    pzbil = pzbil + pzrem;

    // If on purpose, add here the spectator nuclei:	
    if(calincl->bulletType() < 0 && ps->a_projspec != 0) {
      pxbil=pxbil+ps->p1_projspec;
      pybil=pybil+ps->p2_projspec;
      pzbil=pzbil+ps->p3_projspec;
    }

    if((std::fabs(pzbil - pbeam) > 5.0) || (std::sqrt(std::pow(pxbil,2) + std::pow(pybil,2)) >= 3.0)) {
      if(verboseLevel > 3) {
	G4cout <<"Bad momentum conservation after INCL:" << G4endl;
	G4cout <<"delta Pz = " << std::fabs(pzbil - pbeam) << G4endl;
	G4cout <<"      Pt = " << std::sqrt(std::pow(pxbil, 2) + std::pow(pybil, 2)) << G4endl;
      }
    }
       
    volant->iv = 0;   // init du compteur des part evaporees
    varntp->kfis = 0;  //drapeau de fission copie dans le ntuple
    varntp->estfis = 0.0;
    varntp->izfis = 0;
    varntp->iafis = 0;

    // on recopie le remnant dans le ntuple
    varntp->ntrack = varntp->ntrack + 1;
    varntp->massini = iarem;
    varntp->mzini = izrem;
    varntp->exini = esrem;
    varntp->pxrem = pxrem;
    varntp->pyrem = pyrem;
    varntp->pzrem = pzrem;
    varntp->mcorem = mcorem;
    varntp->erecrem = pcorem;
    varntp->erecrem = erecrem;

#ifdef G4INCLDEBUG
    theLogger->fillHistogram1D("bimpact", varntp->bimpact);
#endif

#ifdef G4INCLDEBUG
    theLogger->fillHistogram1D("mzini", varntp->mzini);
#endif
  }
  if(nopart == -2) {
    varntp->ntrack = -2; //FIX: Error flag to remove events containing unphysical events (Ekin > Ebullet).
    evaporationResult->ntrack = -2; //FIX: Error flag to remove events containing unphysical events (Ekin > Ebullet).
  }
  else if(nopart == -1) {
    varntp->ntrack = -1;
    evaporationResult->ntrack = -1;
  }
  if(verboseLevel > 2) {
    G4cout << __FILE__ << ":" << __LINE__ << "Dump varntp after combining: " << G4endl;
    varntp->dump();
  }
}


void G4Incl::processEventInclAbla(G4InclInput *input, G4int eventnumber)
{
  //G4cout <<"Starting event " << eventnumber << G4endl;
  if(input == 0) {
    G4cerr <<"G4Incl fatal error: NULL pointer passed as input!" << G4endl;
    return;
  }
  calincl = input;

  const G4double uma = 931.4942;
  const G4double melec = 0.511;
  const G4double fmp = 938.2796;

  G4double pcorem = 0.0;
  G4double pxrem = 0.0;
  G4double pyrem = 0.0;
  G4double pzrem = 0.0;

  G4double ap = 0.0, zp = 0.0, mprojo = 0.0, pbeam = 0.0;

  varntp->clear();

  if(calincl->bulletType() == -12) {
    be->ia_be = std::abs(calincl->bulletType());
    be->iz_be = 6;
  } else if(calincl->bulletType() == -666) {
    be->iz_be = calincl->extendedProjectileZ();
    be->ia_be = calincl->extendedProjectileA();
  }

  if(calincl->isExtendedProjectile() == false && calincl->bulletType() < -max_a_proj) {
  //  if(calincl->bulletType() < -max_a_proj) {
    G4cout <<"max a of composite projectile is: " << max_a_proj << G4endl;
    exit(0);
  }
  if(calincl->bulletType() < 0) {
    //    calincl->bulletType() = std::floor(calincl->bulletType() + 0.1); WTF???
    be->pms_be=100.;
    G4int i_tabled=0;
    if(be->iz_be == 3 && be->ia_be == 6) {
      be->rms_be=2.56;
      be->bind_be=32.0;
      i_tabled=1;
    } else if(be->iz_be == 3 && be->ia_be == 7) { // TODO: Check the values!
      be->rms_be=2.56;
      be->bind_be=32.0;
      i_tabled=1;
    } else if(be->iz_be == 3 && be->ia_be == 8) {
      be->rms_be=2.40;
      be->bind_be=39.25;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 7) {
      be->rms_be=2.51;
      be->bind_be=58.17;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 9) {
      be->rms_be=2.51;
      be->bind_be=58.17;
      i_tabled=1;
    } else if(be->iz_be == 4 && be->ia_be == 10) {
      be->rms_be=2.45;
      be->bind_be=64.75;
      i_tabled=1;
    } else if(be->iz_be == 5 && be->ia_be == 10) {
      be->rms_be=2.45;
      be->bind_be=64.75;
      i_tabled=1;
    } else if(be->iz_be == 5 && be->ia_be == 11) {
      be->rms_be=2.40;
      be->bind_be=76.21;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 9) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 10) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 11) { // TODO: Check the values!
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && calincl->bulletType() == -12) { // Special Carbon case
      G4cout <<"Carbon 12 (special) selected." << G4endl;
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 6 && be->ia_be == 12) {
      be->rms_be=2.44;
      be->bind_be=92.17;
      i_tabled=1;
    } else if(be->iz_be == 7 && be->ia_be == 16) {
      be->rms_be=2.73;
      be->bind_be=127.62;
      i_tabled=1;
    } else {
      G4cout <<"Warning: No rms and binding for projectile ion A = " << be->ia_be << " Z = " << be->iz_be << G4endl;
      be->rms_be=2.44;
      be->bind_be=92.17;
      G4cout <<"Warning: Using probably bad values rms = " << be->rms_be << " binding = " << be->bind_be << G4endl;
      i_tabled=1;      
    }
      
    if(i_tabled == 0) {
      G4cout <<"This heavy ion (a,z)= " << be->ia_be << " " << be->iz_be << " is not defined as beam in INCL" << G4endl;
      exit(0);
    }
      
    //    G4cout <<"z projectile, rms_r, rms_p (gaussian model)" << be->iz_be << " " << be->rms_be << " " << be->pms_be << G4endl;
    //    G4cout <<"binding energy (mev):" << be->bind_be << G4endl;
    //    G4cout <<"fermi-breakup dresner below a=" << calincl->f[11] << G4endl;
  }      
  //  G4cout <<"Target Mass and Charge: " << calincl->targetA() << " " << calincl->targetZ() << G4endl;
  //  calincl->f[10] = 0; // No clusters

  if(calincl->bulletType() == -12) {  // C12 special case
    mprojo=fmp*std::abs(calincl->bulletType()) - be->bind_be;
    pbeam=std::sqrt(calincl->bulletE()*(calincl->bulletE()+2.*mprojo));
    ap=std::abs(calincl->bulletType());
    zp=be->iz_be;
  } else if(calincl->bulletType() == -666) { // Generic extended projectile
    mprojo=fmp*be->ia_be - be->bind_be;
    pbeam=std::sqrt(calincl->bulletE()*(calincl->bulletE()+2.*mprojo));
    ap=be->ia_be;
    zp=be->iz_be;
  }
  // pi+
  if(calincl->bulletType() == 3) { 
    mprojo = 139.56995;
    ap = 0.0;
    zp = 1.0;
  }

  // pi0
  if(calincl->bulletType() == 4) {
    mprojo = 134.9764;
    ap = 0.0;
    zp = 0.0;
  }

  // pi-
  if(calincl->bulletType() == 5) {
    mprojo = 139.56995;
    ap = 0.0;
    zp = -1.0;
  }

  // coulomb en entree seulement pour les particules ci-dessous

  // proton
  if(calincl->bulletType() == 1) {
    mprojo = 938.27231;
    ap = 1.0;
    zp = 1.0;
  }

  // neutron  
  if(calincl->bulletType() == 2) {
    mprojo = 939.56563;
    ap = 1.0;
    zp = 0.0;
  }

  // deuteron
  if(calincl->bulletType() == 6) {
    mprojo = 1875.61276;
    ap = 2.0;
    zp = 1.0;
  }

  // triton
  if(calincl->bulletType() == 7) {
    mprojo = 2808.95;
    ap = 3.0;
    zp = 1.0;
  }

  // He3
  if(calincl->bulletType() == 8) {
    mprojo = 2808.42;
    ap = 3.0;
    zp = 2.0;
  }

  // Alpha
  if(calincl->bulletType() == 9) {
    mprojo = 3727.42;
    ap = 4.0;
    zp = 2.0;
  }

  // Carbon
  if(calincl->bulletType() == -12) {
    mprojo = 6.0*938.27231 + 6.0*939.56563;
    ap = 12.0;
    zp = 6.0;
  }

  pbeam = std::sqrt(calincl->bulletE()*(calincl->bulletE() + 2.0*mprojo));         

  G4double at = calincl->targetA();
       
  G4int ibert = 1;

  G4int nopart = 0;
  G4int izrem = 0;
  G4int iarem = 0;
  G4double esrem = 0.0;
  G4double erecrem = 0.0;
  G4double berem = 0.0;
  G4double garem = 0.0;
  G4double bimpac = 0.0;
  G4int jrem = 0;
  G4double xjrem = 0.0, yjrem = 0.0, zjrem = 0.0;
  G4double alrem = 0.0;

  // Coulomb barrier
  
  G4double probaTrans = 0.0;
  G4double rndm = 0.0;

  if((calincl->bulletType() == 1) || (calincl->bulletType() >= 6)) {
    //    probaTrans = coulombTransm(calincl->bulletE(),apro,zpro,calincl->targetA(),calincl->targetZ());
    probaTrans = coulombTransm(calincl->bulletE(),ap,zp,calincl->targetA(),calincl->targetZ());
    standardRandom(&rndm, &(hazard->ial));
    if(rndm <= (1.0 - probaTrans)) {
      varntp->ntrack = -1;
      return;
    }
  }

  //  G4cout <<"Before PNU:" << G4endl;
  //  randomGenerator->printSeeds();
  // Call the actual INCL routine:
  pnu(&ibert, &nopart,&izrem,&iarem,&esrem,&erecrem,&alrem,&berem,&garem,&bimpac,
      &jrem, &xjrem, &yjrem, &zjrem);
  //  G4cout <<"After PNU:" << G4endl;
  //  randomGenerator->printSeeds();
  G4double mrem = int(zjrem/197.328); // CHECK
  if (mrem > jrem) mrem=jrem;
  if (mrem < -jrem) mrem=-jrem;

//   nopart=1;
//   kind[0]=1;
//   ep[0]=799.835;
//   alpha[0]=0.08716;
//   beta[0]=0.;
//   gam[0]=0.99619;
//   izrem=82;
//   iarem=208;
//   esrem=200.;
//   erecrem=0.18870;
//   alrem=-0.47101;
//   berem=0.;
//   garem=0.88213;
//   bimpac=2.;
  forceAbsor(&nopart, &iarem, &izrem, &esrem, &erecrem, &alrem, &berem, &garem, &jrem);
  G4double aprf = double(iarem);    // mass number of the prefragment
  G4double jprf = 0.0;                // angular momentum of the prefragment

  // Mean angular momentum of prefragment                                  
  jprf = 0.165 * std::pow(at,(2.0/3.0)) * aprf*(at - aprf)/(at - 1.0);                               
  if (jprf < 0) {
    jprf = 0.0;
  }

  // check m.de jong, ignatyuk, schmidt nuc.phys a 613, pg442, 7th line
  jprf = std::sqrt(2*jprf);

  jprf = jrem;
  varntp->jremn = jrem;      // jrem copie dans le ntuple

  G4double numpi = 0;  // compteurs de pions, neutrons protons
  G4double multn = 0; 
  G4double multp = 0;

  // ecriture dans le ntuple des particules de cascade (sauf remnant)      
  varntp->ntrack = nopart;          // nombre de particules pour ce tir
  if(varntp->ntrack >= VARNTPSIZE) {
    if(verboseLevel > 2) {
      G4cout <<"G4Incl error: Output data structure not big enough." << G4endl;
    }
  }
  varntp->massini = iarem;
  varntp->mzini = izrem;
  varntp->exini = esrem;
  varntp->bimpact = bimpac;
  
  //  three ways to compute the mass of the remnant: 
  //                -from the output of the cascade and the canonic mass
  //                -from energy balance (input - all emitted energies)
  //                -following the approximations of the cugnon code (esrem...)
  G4double mcorem = mprojo + calincl->bulletE() + abla->pace2(double(calincl->targetA()),double(calincl->targetZ()))
    + calincl->targetA()*uma - calincl->targetZ()*melec;

  G4double pxbil = 0.0;
  G4double pybil = 0.0;
  G4double pzbil = 0.0;         

  if(nopart > -1) {
    // Fill the projectile spectator variables
    varntp->masp = ps->a_projspec;
    varntp->mzsp = ps->z_projspec;
    varntp->exsp = ps->ex_projspec;
    varntp->spectatorP1 = ps->p1_projspec;
    varntp->spectatorP2 = ps->p2_projspec;
    varntp->spectatorP3 = ps->p3_projspec;
    varntp->spectatorT = ps->t_projspec;
    for(G4int j = 0; j <= nopart; j++) {
      if(ep[j] < 0.0) { // Workaround to avoid negative energies (and taking std::sqrt of a negative number).
	G4cout <<"G4Incl: Not registering particle with energy: " << ep[j] << G4endl;
	continue;
      }
      if(kind[j] == 0) continue; // Empty particle rows are sometimes produced by lurking indexing problems. We can simply skip those "bad" entries...
      if(gam[j] > CLHEP::pi) {
	if(verboseLevel > 2) {
	  G4cout <<"G4Incl: Just avoided floating point exception by using an ugly hack..." << G4endl;
	}
	continue; // Avoid floating point exception
      }

      varntp->itypcasc[j] = 1; // Particle was produced by the cascade
      // Spectators of composite projectiles (7/2006, AB)
      // (kind is negative in that case)
      if(kind[j] <= 0) { // Particle is a projectile spectator that comes directly from the cascade
	kind[j] *= -1;
	varntp->itypcasc[j]=-1;
	//	G4cout <<"Spectator registered!" << G4endl;
	//	continue;
      }
	
      // kind(): 1=proton, 2=neutron, 3=pi+, 4=pi0, 5=pi -      
      if(kind[j] == 1) { 
	varntp->avv[j] = 1;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+1876.5592)); // cugnon
	multp = multp + 1;
	mcorem = mcorem - ep[j] - 938.27231;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Proton produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 2) { 
	varntp->avv[j] = 1;
	varntp->zvv[j] = 0;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+1876.5592)); // cugnon
	multn = multn + 1;
	mcorem = mcorem - ep[j] - 939.56563;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Neutron produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }
      
      if(kind[j] == 3) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 139.56995;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi+ produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 4) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = 0;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 134.9764;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi0 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 5) { 
	varntp->avv[j] = -1;
	varntp->zvv[j] = -1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j]+276.0)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 139.56995;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Pi+ produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 6) { 
	varntp->avv[j] = 2;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*1874.34)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2806.359;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Deuteron produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 7) { 
	varntp->avv[j] = 3;
	varntp->zvv[j] = 1;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*2806.359)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2806.359;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Triton produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 8) { 
	varntp->avv[j] = 3;
	varntp->zvv[j] = 2;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*2807.119)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 2807.119;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: He3 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      if(kind[j] == 9) { 
	varntp->avv[j] = 4;
	varntp->zvv[j] = 2;
	varntp->plab[j] = std::sqrt(ep[j]*(ep[j] + 2.0*3724.818)); // cugnon
	numpi = numpi + 1;
	mcorem = mcorem - ep[j] - 3724.818;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: He3 produced! " << G4endl;
	  G4cout <<"G4Incl: Momentum: "<< varntp->plab[j] << G4endl;
	}
      }

      varntp->enerj[j] = ep[j];
      varntp->tetlab[j] = 180.0*std::acos(gam[j])/3.141592654;
      varntp->philab[j] = 180.0*std::atan2(beta[j],alpha[j])/3.141592654;
      pxbil = pxbil + varntp->plab[j]*alpha[j];
      pybil = pybil + varntp->plab[j]*beta[j];
      pzbil = pzbil + varntp->plab[j]*gam[j];

      if(verboseLevel > 3) {
	G4cout <<"Momentum: " << varntp->plab[j]   << G4endl;
	G4cout <<"Theta: "    << varntp->tetlab[j] << G4endl;
	G4cout <<"Phi: "      << varntp->philab[j] << G4endl;
      }
    }

    // calcul de la masse (impulsion) du remnant coherente avec la conservation d'energie:
    pcorem = std::sqrt(erecrem*(erecrem + 2.0 * 938.2796 * iarem));   // cugnon
    mcorem = 938.2796 * iarem;                               // cugnon
    varntp->pcorem = pcorem;
    varntp->mcorem = mcorem;
    // Note: Il faut negliger l'energie d'excitation (ESREM) pour que le bilan 
    // d'impulsion soit correct a la sortie de la cascade.....et prendre la
    // masse MCOREM comme ci-dessus (fausse de ~1GeV par rapport aux tables...)        
    pxrem = pcorem * alrem;
    pyrem = pcorem * berem;
    pzrem = pcorem * garem;
    varntp->pxrem = pxrem;
    varntp->pyrem = pyrem;
    varntp->pzrem = pzrem;
    pxbil = pxbil + pxrem;
    pybil = pybil + pyrem;
    pzbil = pzbil + pzrem;

    // If on purpose, add here the spectator nuclei:	
    if(calincl->bulletType() < 0 && ps->a_projspec != 0) {
      pxbil=pxbil+ps->p1_projspec;
      pybil=pybil+ps->p2_projspec;
      pzbil=pzbil+ps->p3_projspec;
    }

    if((std::fabs(pzbil - pbeam) > 5.0) || (std::sqrt(std::pow(pxbil,2) + std::pow(pybil,2)) >= 3.0)) {
      if(verboseLevel > 3) {
	G4cout <<"Bad momentum conservation after INCL:" << G4endl;
	G4cout <<"delta Pz = " << std::fabs(pzbil - pbeam) << G4endl;
	G4cout <<"      Pt = " << std::sqrt(std::pow(pxbil, 2) + std::pow(pybil, 2)) << G4endl;
      }
    }
       
    volant->iv = 0;   // init du compteur des part evaporees
    varntp->kfis = 0;  //drapeau de fission copie dans le ntuple
    varntp->estfis = 0.0;
    varntp->izfis = 0;
    varntp->iafis = 0;

    // on recopie le remnant dans le ntuple
    varntp->ntrack = varntp->ntrack + 1;
    varntp->massini = iarem;
    varntp->mzini = izrem;
    varntp->exini = esrem;
    varntp->pxrem = pxrem;
    varntp->pyrem = pyrem;
    varntp->pzrem = pzrem;
    varntp->mcorem = mcorem;
    varntp->erecrem = pcorem;
    varntp->erecrem = erecrem;

    if(verboseLevel > 2) {
      G4cout << __FILE__ << ":" << __LINE__ << "Dump varntp after cascade: " << G4endl;
      varntp->dump();
    }
    // Maximum remnant for Geant4 Fermi break-up: A = 17, Z = 8
    if(varntp->massini < 17 && varntp->exini < 1000 && varntp->massini > 0 && varntp->mzini > 0 && useFermiBreakup) { // Choose between Fermi break-up and ABLA
      varntp->needsFermiBreakup = true;
    } else {
      varntp->needsFermiBreakup = false;
      // Evaporation/fission:
      evaporationResult->ntrack = 0;
      abla->breakItUp(G4int(varntp->massini), G4int(varntp->mzini), mcorem, varntp->exini, varntp->jremn,
      		      erecrem, pxrem, pyrem, pzrem, eventnumber);

      if(verboseLevel > 2) {
	G4cout << __FILE__ << ":" << __LINE__ << "Dump evaporationResult after evaporation: " << G4endl;
	evaporationResult->dump();
      }
      for(G4int evaporatedParticle = 0; evaporatedParticle < evaporationResult->ntrack; evaporatedParticle++) {
	if(evaporationResult->avv[evaporatedParticle] == 0 && evaporationResult->zvv[evaporatedParticle] == 0) { //Fix: Skip "empty" particles with A = 0 and Z = 0
	  // continue;
	}
	varntp->kfis = evaporationResult->kfis;
	varntp->itypcasc[varntp->ntrack] = evaporationResult->itypcasc[evaporatedParticle];
	varntp->avv[varntp->ntrack] = evaporationResult->avv[evaporatedParticle];
	varntp->zvv[varntp->ntrack]= evaporationResult->zvv[evaporatedParticle];
	varntp->plab[varntp->ntrack] = evaporationResult->plab[evaporatedParticle];
	varntp->enerj[varntp->ntrack] = evaporationResult->enerj[evaporatedParticle];
	varntp->tetlab[varntp->ntrack] = evaporationResult->tetlab[evaporatedParticle];
	varntp->philab[varntp->ntrack] = evaporationResult->philab[evaporatedParticle];
	varntp->ntrack++;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Returning evaporation result"           << G4endl;
	  G4cout <<"G4Incl: A = " << varntp->avv[varntp->ntrack]    << " Z = " << varntp->zvv[varntp->ntrack] << G4endl;
	  G4cout <<"Energy: "     << varntp->enerj[varntp->ntrack]  << G4endl;
	  G4cout <<"Momentum: "   << varntp->plab[varntp->ntrack]   << G4endl;
	  G4cout <<"Theta: "      << varntp->tetlab[varntp->ntrack] << G4endl;
	  G4cout <<"Phi: "        << varntp->philab[varntp->ntrack] << G4endl;
	}
      }
      if(verboseLevel > 2) {
	G4cout <<"G4Incl: ntrack = "          << varntp->ntrack << G4endl;
	G4cout <<"G4Incl: Done extracting..." << G4endl;
      }
    }
    if(calincl->bulletType() < 0 && useProjSpect && !useFermiBreakup) { // If we use projectile spectators for carbon but no fermi breakup
      // Evaporation/fission:
      evaporationResult->ntrack = 0;
      //      G4cout <<"Warning: Using ABLA to de-excite projectile spectator..." << G4endl;
      abla->breakItUp(G4int(varntp->masp), G4int(varntp->mzsp), ps->m_projspec, varntp->exsp, 0.0,
      		      ps->t_projspec, ps->p1_projspec, ps->p1_projspec, ps->p1_projspec, eventnumber);

      if(verboseLevel > 2) {
	G4cout << __FILE__ << ":" << __LINE__ << "Dump evaporationResult after evaporation: " << G4endl;
	evaporationResult->dump();
      }
      for(G4int evaporatedParticle = 0; evaporatedParticle < evaporationResult->ntrack; evaporatedParticle++) {
	if(evaporationResult->avv[evaporatedParticle] == 0 && evaporationResult->zvv[evaporatedParticle] == 0) { //Fix: Skip "empty" particles with A = 0 and Z = 0
	}
	varntp->kfis = evaporationResult->kfis;
	varntp->itypcasc[varntp->ntrack] = evaporationResult->itypcasc[evaporatedParticle];
	varntp->avv[varntp->ntrack] = evaporationResult->avv[evaporatedParticle];
	varntp->zvv[varntp->ntrack]= evaporationResult->zvv[evaporatedParticle];
	varntp->plab[varntp->ntrack] = evaporationResult->plab[evaporatedParticle];
	varntp->enerj[varntp->ntrack] = evaporationResult->enerj[evaporatedParticle];
	varntp->tetlab[varntp->ntrack] = evaporationResult->tetlab[evaporatedParticle];
	varntp->philab[varntp->ntrack] = evaporationResult->philab[evaporatedParticle];
	varntp->ntrack++;
	if(verboseLevel > 3) {
	  G4cout <<"G4Incl: Returning evaporation result"           << G4endl;
	  G4cout <<"G4Incl: A = " << varntp->avv[varntp->ntrack]    << " Z = " << varntp->zvv[varntp->ntrack] << G4endl;
	  G4cout <<"Energy: "     << varntp->enerj[varntp->ntrack]  << G4endl;
	  G4cout <<"Momentum: "   << varntp->plab[varntp->ntrack]   << G4endl;
	  G4cout <<"Theta: "      << varntp->tetlab[varntp->ntrack] << G4endl;
	  G4cout <<"Phi: "        << varntp->philab[varntp->ntrack] << G4endl;
	}
      }
    }
#ifdef G4INCLDEBUG
    theLogger->fillHistogram1D("bimpact", varntp->bimpact);
#endif

#ifdef G4INCLDEBUG
    theLogger->fillHistogram1D("mzini", varntp->mzini);
#endif
  }
  if(nopart == -2) {
    varntp->ntrack = -2; //FIX: Error flag to remove events containing unphysical events (Ekin > Ebullet).
    evaporationResult->ntrack = -2; //FIX: Error flag to remove events containing unphysical events (Ekin > Ebullet).
  }
  else if(nopart == -1) {
    varntp->ntrack = -1;
    evaporationResult->ntrack = -1;
  }
  if(verboseLevel > 2) {
    G4cout << __FILE__ << ":" << __LINE__ << "Dump varntp after combining: " << G4endl;
    varntp->dump();
  }
}

// Initialization routines
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

  G4double xrand = 0.0;
  G4long ialdep = 0;
  G4int imat = 0;
  G4int iamat = 0, izmat = 0;

  // for the 19 secondary seeds of hazard:
  G4int nbtirhaz[IGRAINESIZE] = {38,82,76,18,39,31,41,59,26,54,
				 14,84,13,15,91,89,10,6,52};

  // specific parameters for incl:	 
  // espace de phases test (r et p) pour pauli: 
  // valeur recommandee par j.c. v-test=0.592 h**3:
  G4double rbl = 2.;

  // valeur pour avoir v-test=2 h**3 (avec pbl=200)
  rbl = 3.1848;

  // preparation of 19 other seeds (can also be initialized from outside):
  if(initRandomSeed) {
    ialdep=hazard->ial;
    for(G4int i = 0; i < IGRAINESIZE; i++) {
      for(G4int j = 0; j < nbtirhaz[i]; j++) {
	standardRandom(&xrand,&(hazard->ial));
      }

      // Zero is not accepted as random seed!
      do {
	standardRandom(&xrand,&(hazard->ial));
      } while(xrand == 0);

      xrand = xrand * 100000;

      while(xrand < 10000) {
	xrand = xrand * 10;
      }
      hazard->igraine[i] = (int) xrand;
      if(hazard->igraine[i] == ((hazard->igraine[i] / 2) * 2)) {
	hazard->igraine[i] = hazard->igraine[i] + 1;
      }
    }

    hazard->ial = ialdep;
  }

  // calculation with realistic nuclear density (saxon-wood)
  if (ws->nosurf <= 0) {
    // prepare nucleus density for nbmat nucleus defined in struct mat
    if(mat->nbmat >= 500) {
      if(verboseLevel > 2) {
	G4cout <<"You need " << mat->nbmat << " nuclei in your problem. The maximum number of nuclei is 500 " << G4endl;
      }
      return;
    }

    for(G4int i = 0; i < mat->nbmat; i++) {
      imat = i;
      izmat = int(mat->zmat[i]);
      iamat = int(mat->amat[i]);

      initMaterial(izmat, iamat, imat);
    }
  }

  // deuteron density in momentum space:
  densDeut();
}


void G4Incl::initMaterial(G4int izmat, G4int iamat, G4int imat)
{
  G4double res_dws = 0.0;
  G4double fnor = 0.0;

  G4double rcour = 0.0, geom = 0.0;
  G4int nbr = 0;

  G4double step = 0.0, f_r = 0.0;

  // rms espace r, espace p, fermi momentum and energy for light gauss nuc.      
  const G4double datarms1t[LGNSIZE] = {0.0, 0.0, 0.0, 0.0, 0.0, 2.10, 1.80, 1.80, 1.63};
  const G4double datapf1t[LGNSIZE] = {0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 110.0, 110.0, 153.0};

  for(G4int i = 0; i < LGNSIZE; i++) {
    light_gaus_nuc->rms1t[i] = datarms1t[i];
    light_gaus_nuc->pf1t[i] = datapf1t[i];
  }

  // fermi 2 param from a=19 to 28, modified harm oscil a=6 to 18
  // (h. de vries et al. at. data and nuc. data tab. 36 (1987) 495)
  const G4double datarln[LNSIZE] = {0.0,0.0,0.0,0.0,0.0,0.334,0.327,0.479,0.631,0.838,
				    0.811,1.07,1.403,1.335,1.25,1.544,1.498,1.513,
				    2.58,2.77, 2.775,2.78,2.88,2.98,3.22,3.03,2.84,
				    3.14,0.0,0.0};

  const G4double dataaln[LNSIZE] = {0.0,0.0,0.0,0.0,0.0,1.78,1.77,1.77,1.77,1.71,
				    1.69,1.69,1.635,1.730,1.81,1.833,1.798,
				    1.841,0.567,0.571, 0.560,0.549,0.550,0.551,
				    0.580,0.575,0.569,0.537,0.0,0.0};

  for(G4int i = 0; i < LNSIZE; i++) {
    light_nuc->r[i] = datarln[i];
    light_nuc->a[i] = dataaln[i];
  }

  if(verboseLevel > 3) {
    G4cout <<"Nuclear density for nucleus (z, a): " << izmat << " " << iamat << " " << imat << G4endl;
  }
  
  const G4double fmp = 938.2796; // From INCL data

  // parametres moyens de densite de la cible (fermi 2 parametres)
  if (iamat >= 28) {
    ws->r0 = (2.745e-4*iamat+1.063)*std::pow(G4double(iamat), 0.33333333);
    ws->adif = 1.63e-4*iamat+0.510;
    ws->rmaxws = ws->r0 + (ws->xfoisa)*(ws->adif);
  }
  else if(iamat >= 19) {
    ws->r0 = light_nuc->r[iamat-1];
    ws->adif = light_nuc->a[iamat-1];
    ws->rmaxws = ws->r0 + (ws->xfoisa)*(ws->adif);
  }
  else if(iamat >= 6) {
    ws->r0 = light_nuc->r[iamat-1];
    ws->adif = light_nuc->a[iamat-1];
    ws->rmaxws = 5.5 + 0.3*(double(iamat)-6.)/12.;
  }
  else if(iamat >= 2) {
    if(iamat == 2) {
      ws->r0=light_gaus_nuc->rms1t[5]; // Orig: rms1t(6)
      light_gaus_nuc->pfln[5] = light_gaus_nuc->pf1t[5]*1.291;  // Orig [6], std::sqrt(5/3)=1.291
      light_gaus_nuc->tfln[5] = std::sqrt(std::pow(light_gaus_nuc->pfln[5],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[5] = light_gaus_nuc->tfln[5] + 2.22;
      if(verboseLevel > 2) {
	G4cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[5] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[5] << " " << light_gaus_nuc->tfln[5] << G4endl;
      }
    }
    if(iamat == 3 && izmat == 1) {
      ws->r0=light_gaus_nuc->rms1t[6]; // Orig: rms1t(7)
      light_gaus_nuc->pfln[6] = light_gaus_nuc->pf1t[6]*1.291;  // Orig [7], std::sqrt(5/3)=1.291
      light_gaus_nuc->tfln[6] = std::sqrt(std::pow(light_gaus_nuc->pfln[6],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[6] = light_gaus_nuc->tfln[6] + 4.24;
      if(verboseLevel > 2) {
	G4cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[6] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[6] << " " << light_gaus_nuc->tfln[6] << G4endl;
      }
    }
    if(iamat == 3 && izmat == 2) {
      ws->r0 = light_gaus_nuc->rms1t[7]; // Orig: rms1t(8)
      light_gaus_nuc->pfln[7] = light_gaus_nuc->pf1t[7]*1.291;   //!std::sqrt(5/3)=1.291
      light_gaus_nuc->tfln[7] = std::sqrt(std::pow(light_gaus_nuc->pfln[7],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[7] = light_gaus_nuc->tfln[7] + 3.86;
      if(verboseLevel > 2) {
	G4cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[7] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[7] << " " << light_gaus_nuc->tfln[7] << G4endl;
      }
    }
    if(iamat == 4) {
      ws->r0 = light_gaus_nuc->rms1t[8]; // Orig: rms1t(9)
      light_gaus_nuc->pfln[8] = light_gaus_nuc->pf1t[8]*1.291;   // !std::sqrt(5/3)=1.291
      light_gaus_nuc->tfln[8] = std::sqrt(std::pow(light_gaus_nuc->pfln[8],2) + fmp*fmp) - fmp;
      light_gaus_nuc->vnuc[8] = light_gaus_nuc->tfln[8] + 9.43;
      if(verboseLevel > 2) {
	G4cout <<"Nuclear potential: " << light_gaus_nuc->vnuc[8] << "MeV, Fermi momentum and energy: " << light_gaus_nuc->pfln[8] << " " << light_gaus_nuc->tfln[8] << G4endl;
      }
    }
    ws->adif = 0.57735*ws->r0;
    ws->rmaxws = ws->r0 + 2.5;
  }
  ws->drws = (ws->rmaxws)/29.0;

  // bmax for sigma geom and various projectiles (p,n,pion/d/t/he3/he4/)
  G4int j = 0;
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
      //	 G4double cte = std::pow(ws->adif,3)*std::sqrt(2.*3.141592654);        
      res_dws = 3.0*(std::pow(ws->adif, 3)*std::sqrt(2.0*3.141592654))/2.0;
    }
  }
  fnor = res_dws;

  // calcul de q/pf=f(r)      
  nbr = int(std::floor((ws->rmaxws)/(ws->drws) + 1.5));
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
	f_r = std::pow(f_r,(1./3.));
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

  if(verboseLevel > 3) {
    if(iamat >= 19) {
      G4cout <<"Wood-Saxon density, r0 = " << ws->r0 << " a = " << ws->adif << G4endl;
    }
    if(iamat >= 6 && iamat <= 19) {
      G4cout <<"Modif. harm. oscil. density, alpha = " << ws->r0 << " a = " << ws->adif << G4endl;
    }
    if(iamat >= 2 && iamat <= 6) {
      G4cout <<"Gaussian density, r.m.s = " << ws->r0 << " sigma = " << ws->adif << G4endl;
    }
  }
  
  geom = 31.41592653*std::pow(ws->rmaxws,2);

  if(verboseLevel > 3) {
    G4cout <<"For incident nucleons or pions rmax = " << ws->rmaxws << " and geometrical (pi*rmaxws*rmaxws) reaction cross section (mb) is " << geom << G4endl;
    for(G4int k = 2; k < MATGEOSIZE; k++) {
      G4cout << "Rmaxws for d/t/3he/4he = " << mat->bmax_geo[k][imat] << G4endl;
    }

    G4cout <<"Exact calculation of the r(q) function for the target nucleus density q/pf  r(q/pf)" << G4endl;
  }
}

G4double G4Incl::deutv(G4int l, G4double q)
{
  G4double res = 0.0;

  if (l == 0) {
    for(G4int i = 0; i < DTONSIZE; i++) {
      res = res + dton->c[i]/(std::pow(q,2) + fm2(i+1));
    }
  }
  if(l != 0) {
    for(G4int i = 0; i < DTONSIZE; i++) { 
      res = res + dton->d[i]/(std::pow(q,2) + fm2(i+1));
    }
  }

  return res*std::sqrt(2./CLHEP::pi)*dton->fn; // See G4InclDataDefs.hh
}

G4double G4Incl::fm2(G4int j)
{
  /**
   * In implementation Returns the values of the function:
   * \f[
   * (0.23162461 + (j - 1))^2
   * \f]
   * @param j an integer parameter
   * @return a double value
   */
  
  return std::pow((0.23162461 + (j - 1)),2);
}

G4double G4Incl::interpolateFunction(G4double xv)
{
  // fonction d'interpolation au point xv ( meme hors bornes )             
  // de la fn x->y dont les derivees premieres (s) ont ete                 
  // evaluees par l'appel prealable de flin2                              
  // les indices vont de 1 a n

  saxw->k = saxw->imat;
  G4double tz = xv - saxw->x[0][saxw->imat];
  G4int j = 0;
    
  if(tz < 0) {
#ifdef G4INCLDEBUG
  theLogger->fillHistogram1D("interpolationResult", (saxw->y[0][saxw->imat] + saxw->s[0][saxw->imat]*tz));
  theLogger->fillHistogram2D("interpolationPoints", (saxw->y[0][saxw->imat] + saxw->s[0][saxw->imat]*tz), xv);
#endif
    return (saxw->y[0][saxw->imat] + saxw->s[0][saxw->imat]*tz);
  }
  else if(tz == 0) {
#ifdef G4INCLDEBUG
  theLogger->fillHistogram1D("interpolationResult", (saxw->y[0][saxw->imat]));
  theLogger->fillHistogram2D("interpolationPoints", (saxw->y[0][saxw->imat]), xv);
#endif
  return (saxw->y[0][saxw->imat]);
  }
  else { // tz > 0
    for(G4int i = 1; i < saxw->n; i++) {
      j = i;
      tz = xv - saxw->x[j][saxw->imat];
      if(tz <= 0) {
	break;
      }
    }
    if(tz >= 0) {
#ifdef G4INCLDEBUG
      theLogger->fillHistogram1D("interpolationResult", saxw->y[j][saxw->imat]);
      theLogger->fillHistogram2D("interpolationPoints", saxw->y[j][saxw->imat], xv);
#endif
      return saxw->y[j][saxw->imat];
    } else if(tz < 0.0) {
      j = j - 1;
      G4double dgx = xv - saxw->x[j][saxw->imat];
#ifdef G4INCLDEBUG
      theLogger->fillHistogram1D("interpolationResult", (saxw->y[j][saxw->imat] + saxw->s[j][saxw->imat]*dgx));
      theLogger->fillHistogram2D("interpolationPoints", (saxw->y[j][saxw->imat] + saxw->s[j][saxw->imat]*dgx), xv);
#endif
      return(saxw->y[j][saxw->imat] + saxw->s[j][saxw->imat]*dgx);
    }
  }

  return 0.0;
}

void G4Incl::firstDerivative(G4int k)
{
  for(G4int i=0; i < saxw->n-1; i++) {
    if((saxw->x[i+1][k] - saxw->x[i][k]) == 0.0) { // Safeguard to avoid division by zero
      saxw->s[i][k] = 0.0;
      continue;
    }
    saxw->s[i][k] = (saxw->y[i+1][k] - saxw->y[i][k]) / (saxw->x[i+1][k] - saxw->x[i][k]);
  }
  saxw->s[saxw->n-1][k] = saxw->s[saxw->n-2][k];
}

G4double G4Incl::wsax(G4double r) {
  return std::pow(r,2) / (1.0+std::exp(r-(ws->r0)/(ws->adif)));
}

G4double G4Incl::derivWsax(G4double r)
{
  G4double derivwsax = std::pow(r,3)*std::exp((r-(ws->r0))/(ws->adif))/std::pow((1.0+std::exp((r-(ws->r0))/(ws->adif))),2);
  return derivwsax/(ws->adif);
}

G4double G4Incl::dmho(G4double r)
{
  G4double arg=std::pow((r/(ws->adif)),2);
  return r*r*(1.+(ws->r0)*arg)*std::exp(-arg);
}

G4double G4Incl::derivMho(G4double r)
{
  G4double arg=std::pow((r/(ws->adif)),2);
  return -2.*r*r*arg*((ws->r0) -1.-(ws->r0)*arg)*std::exp(-arg);
}

G4double G4Incl::derivGaus(G4double r)
{
  G4double arg=std::pow((r/(ws->adif)),2);
  return r*r*arg*std::exp(-arg/2.);      
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
  for(G4int init_i = 0; init_i < qsize; init_i++) {
    q[init_i] = 0.0;
    f[init_i] = 0.0;
  }

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

  if(verboseLevel > 3) {
    G4cout << "deuteron density in q space from Paris potential: " << spl2->n << " Exact values from 0 to " 
	      << q[spl2->n-1] << " fm-1 " << G4endl;
  }
}

G4double G4Incl::integrate(G4double ami, G4double ama, G4double step, G4int functionChoice)
{
  G4double res = 0.0;
  G4double x1[5];
  for(G4int init_i = 0; init_i < 5; init_i++) {
    x1[init_i] = 0.0;
  }
  G4double ri = ami;
  G4double ra = ama;
  G4int nb = 0;
  G4double acont = 1.0;
  G4double dr = step;

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
  nb = int(std::floor(((ra - ri)/dr + 1.0000000001))); // 1.0000000001 -> 0.0
  dr = (ra - ri)/(double(nb - 1)); 
  res = 0.0;

  if(nb < 10) {
    if(verboseLevel > 2) {
      G4cout <<"pas assez de points d'integration" << G4endl;
    }
    return 0.0;
  }
  
  for(G4int i = 0; i < 5; i++) {
    res = res + (callFunction(functionChoice, ri) + callFunction(functionChoice, ra))*x1[i];
    ri = ri + dr;
    ra = ra - dr;
  }

  nb = nb - 10;

  if(nb == 0) {
    return (res*dr*acont);
  }

  for(G4int i = 0; i < nb; i++) {
    res = res + callFunction(functionChoice, ri);
    ri = ri + dr;
  }

  return(res*dr*acont);
}


G4double G4Incl::dens(G4double q)
{
  return q*q*(std::pow(deutv(0,q),2)+std::pow(deutv(2,q),2));
}

void G4Incl::spl2ab() 
{
  G4int i = 0, j = 0, k = 0;

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
  G4double tz = xv-spl2->x[0];
  G4int j;

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
  if(verboseLevel > 2) {
    G4cout <<"G4Incl::splineab : requested point outside defined region! Returning 0.0." << G4endl;
  }
  return 0.0; 
}

// Actual calculation

void G4Incl::pnu(G4int *ibert_p, G4int *nopart_p, G4int *izrem_p, G4int *iarem_p, G4double *esrem_p,
		 G4double *erecrem_p, G4double *alrem_p, G4double *berem_p, G4double *garem_p,
		 G4double *bimpact_p, G4int *l_p, G4double *xjrem, G4double *yjrem, G4double *zjrem)
{
  //  G4cout <<"Now in pnu..." << G4endl;
  //  G4cout <<"(clean-up)" << G4endl;
  G4int npenter = 0;
  G4int nnenter = 0;
  G4int n_enter_pot = 0;
  G4int avatarCounter = 0;
  ps->clear(); // For projectile spectators
  //  G4int i_c = 0;
  G4double p1spec=0.;
  G4double p2spec=0.;
  G4double p3spec=0.;
  //  G4double p_spec2 = 0.0;
  //  G4double s_spec = 0.0;
  // G4double e_spec = 0.0;
  G4double xl1_i, xl2_i, xl3_i;
  G4int idq = 0;
  G4double destar = 0.0;
  //  G4cout <<"(checkpoint 'setipszero0)" << G4endl;
  G4int ips = 0;
  //  G4double p1_spec = 0.0;
  //  G4double p2_spec = 0.0;
  //  G4double p3_spec = 0.0;
  G4int n_activnuc = 0;
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

  G4double minus_b1 = 0.0, minus_b2 = 0.0, minus_b3 = 0.0;
  G4double aml1 = 0.0;
  G4double aml2 = 0.0;
  G4double amlnew = 0.0;
  G4double arg = 0.0;
  G4double b1 = 0.0;
  G4double b2 = 0.0;
  G4double b3 = 0.0;
  G4double bb2 = 0.0;
  G4double bevar = 0.0;
  G4double bmass[2000]; 
  for(G4int init_i = 0; init_i < 2000; init_i++) {
    bmass[init_i] = 0.0;
  }
  G4double bmax2 = 0.0;
  G4double c1 = 0.0;
  G4double c2 = 0.0;
  G4double cb0 = 0.0;
  G4double cchi = 0.0;
  G4double ccr = 0.0;
  G4double cg = 0.0;
  G4double cif = 0.0;
  G4double cmultn = 0.0;
  G4double cobe = 0.0;
  G4double coeffb0 = 0.0;
  G4double comom = 0.0;
  G4double cstet = 0.0;
  G4double dis1 = 0.0;
  G4double dis2 = 0.0;
  G4double dis3 = 0.0;
  G4double dist = 0.0;
  G4double eb0 = 0.0;
  G4double ecoreh5 = 0.0;
  G4double efer = 0.0;
  G4double egs = 0.0;
  G4double eh5 = 0.0;
  G4double eh6 = 0.0;
  G4double eij = 0.0;
  G4double ekout = 0.0;
  G4double elead = 0.0;
  G4double energie_in = 0.0;
  G4double ener_max = 0.0;
  G4double eout = 0.0;
  //  G4double eps_c[BL1SIZE]; 
  //  G4double p3_c[BL1SIZE];
  //  for(G4int init_i = 0; init_i < BL1SIZE; init_i++) {
    //    eps_c[init_i] = 0.0;
    //    p3_c[init_i] = 0.0;
  //  }
  G4double epsv = 0.0;
  G4double erecg = 0.0;
  G4double erem = 0.0;
  G4double exi = 0.0;
  G4double expob0 = 0.0;
  G4double factemp = 0.0;
  G4double fffc = 0.0;
  G4double fm = 0.0;
  G4double g1 = 0.0;
  G4double g2 = 0.0;
  G4double ge = 0.0;
  G4double geff = 0.0;
  G4double gg = 0.0;
  G4double gl1 = 0.0;
  G4double gl2 = 0.0;
  G4double gpsg = 0.0;
  G4int i1 = 0;
  G4int i20 = 0;
  G4int ic33 = 0;
  G4int ich1 = 0;
  G4int ich2 = 0;
  G4int ich3 = 0;
  G4int ich4 = 0;
  G4int ichd = 0;
  G4int ichpion = 0;
  G4int idecf = 0;
  G4int idep = 0;
  G4int iej = 0;
  G4int iejn = 0;
  G4int iejp = 0;
  G4int i_emax = 0;
  G4int iflag = 0;
  G4int iflag20 = 0; 
  G4int iflag40 = 0; 
  G4int iflag60 = 0; 
  G4int ilm = 0; 
  G4int imin = 0;
  G4int indic[3000]; 
  G4int nb_transprojo=0;
  G4double v_proj = 0.0;
  for(G4int init_i = 0; init_i < 3000; init_i++) {
    indic[init_i] = 0;
  }
  G4int inrem = 0;
  G4int ip = 0;
  G4int ipi[2000]; 
  for(G4int init_i = 0; init_i < 2000; init_i++) {
    ipi[init_i] = 0;
  }
  G4int iqe = 0;
  G4int irem = 0;
  G4int irst_avatar = 0;
  G4int isos = 0; 
  G4int itch = 0;
  G4int iteste = 0;
  G4int itt = 0;
  G4int ixr1 = 0;
  G4int ixr2 = 0;
  G4int ixr3 = 0;
  //  G4int k; 
  G4int kcol = 0;
  G4int kd = 0;
  //  G4int klm = 0; 
//   G4int l1; 
//   G4int l2; 
  G4int ldel = 0;
  G4int lead = 0;
  G4int led = 0;
  G4int lnew = 0;
  G4int lp = 0; 
  G4int lp1 = 0;
  G4double mcdd = 0.0;
  //G4double mg; 
  G4int mg = 0;
  G4double mpaul1 = 0.0;
  G4double mpaul2 = 0.0;
  G4double mrdd = 0.0;
  G4double mrdn = 0.0;
  G4double mrdp = 0.0;
  G4double mrnd = 0.0;
  G4double mrnn = 0.0;
  G4double mrpd = 0.0;
  G4int n20 = 0;
  G4int nbalttf = 0;
  G4int nbquit = 0;
  G4int nbtest = 0;
  G4int nc[300];
  G4bool isPartOfSpectatorNucleus[300];
  G4int npproj[300];
  for(G4int init_i = 0; init_i < 300; init_i++) {
    nc[init_i] = 0;
    npproj[init_i] = 0;
  }
  G4int ncol = 0;
  G4int ncol_2c = 0;
  G4int next = 0;
  G4int nmiss = 0;
  G4int np = 0;
  G4int npidir = 0;
  G4int npion = 0;
  G4int npx = 0;
  G4int nsum_col = 0;
  G4double p1v = 0.0;
  G4double p2v = 0.0;
  G4double p3v = 0.0;
  G4double pfrem1 = 0.0;
  G4double pfrem2 = 0.0;
  G4double pfrem3 = 0.0;
  G4double pfreml = 0.0;
  G4double pfreml2 = 0.0;
  G4double phi = 0.0;
  G4double p_mod = 0.0;
  G4double pot = 0.0;
  G4double pout1 = 0.0;
  G4double pout2 = 0.0;
  G4double pout3 = 0.0;
  G4double pppp = 0.0;
  G4double prem1 = 0.0;
  G4double prem2 = 0.0;
  G4double prem3 = 0.0;
  G4double psf = 0.0;
  G4double pspr = 0.0;
  G4double ptotl = 0.0;
  G4double qdeut = 0.0;
  G4double qqq = 0.0;
  G4double r22 = 0.0;
  G4double rcm1 = 0.0;
  G4double rcm2 = 0.0;
  G4double rcm3 = 0.0;
  G4double rcorr = 0.0;
  G4double rhopi = 0.0;
  G4double rndm = 0.0;
  G4double rr = 0.0;
  G4double rrrr = 0.0;
  G4double s = 0.0;
  G4double s1t1 = 0.0;
  G4double s2t1 = 0.0;
  G4double s3t1 = 0.0;
  G4double schi = 0.0;
  G4double sepa = 0.0;
  G4double sif = 0.0;
  G4double sitet = 0.0;
  G4double sp1t1 = 0.0;
  G4double sp2t1 = 0.0;
  G4double sp3t1 = 0.0;
  G4double sq = 0.0;
  G4double sueps = 0.0;
  G4double t[50];
  for(G4int init_i = 0; init_i < 50; init_i++) {
    t[init_i] = 0.0;
  }
  G4double t0 = 0.0;
  G4double t1 = 0.0;
  G4double t2 = 0.0;
  G4double t3 = 0.0;
  G4double t33 = 0.0;
  G4double t4 = 0.0;
  G4double t5 = 0.0;
  G4double t6 = 0.0;
  G4double t7 = 0.0;
  G4double t8 = 0.0;
  G4double tau = 0.0;
  G4double tbid = 0.0;
  G4double tdel = 0.0;
  G4double temfin = 0.0;
  G4double tim = 0.0;
  G4double timi = 0.0;
  G4double tlabu = 0.0;
  G4double tp = 0.0;
  G4double tref = 0.0;
  G4double tri = 0.0;
  G4double tt31 = 0.0;
  G4double tt32 = 0.0;
  G4double tt33 = 0.0;
  G4double tt34 = 0.0;
  G4double tt35 = 0.0;
  G4double tt36 = 0.0;
  G4double tte = 0.0;
  G4double u = 0.0;
  G4double v = 0.0;
  G4double var_ab = 0.0;
  G4double x = 0.0;
  G4double x1l1 = 0.0;
  G4double x1l2 = 0.0;
  G4double x1_target = 0.0; 
  G4double x2_target = 0.0; 
  G4double x3_target = 0.0; 
  G4double x2cour = 0.0;
  G4double x2l1 = 0.0;
  G4double x2l2 = 0.0;
  G4double x3l1 = 0.0;
  G4double x3l2 = 0.0;
  G4double xapres = 0.0;
  G4double xavant = 0.0;
  G4double xbl1 = 0.0;
  G4double xbl2 = 0.0;
  G4double xc = 0.0;
  G4double xe = 0.0;
  G4double xga = 0.0;
  G4double xl1 = 0.0;
  G4double xl2 = 0.0;
  G4double xl3 = 0.0;
  G4double xlab = 0.0;
  G4double xleng = 0.0;
  G4double xlengm = 0.0;
  G4double xpb = 0.0;
  G4double xq = 0.0;
  G4double xr1 = 0.0;
  G4double xr2 = 0.0;
  G4double xr3 = 0.0;
  G4double xr4 = 0.0;
  G4double xr5 = 0.0;
  G4double xr6 = 0.0;
  G4double xr7 = 0.0;
  G4double xr8 = 0.0;
  G4double xv = 0.0;
  G4double xxx = 0.0;
  G4double xy1 = 0.0;
  G4double xy2 = 0.0;
  G4double xy3 = 0.0;
  G4double xye = 0.0;
  G4double y = 0.0;
  G4double q1[BL1SIZE];
  G4double q2[BL1SIZE];
  G4double q3[BL1SIZE];
  G4double q4[BL1SIZE];
  G4double ym[BL1SIZE];
  for(G4int init_i = 0; init_i < BL1SIZE; init_i++) {
    q1[init_i] = 0.0;
    q2[init_i] = 0.0;
    q3[init_i] = 0.0;
    q4[init_i] = 0.0;
    ym[init_i] = 0.0;
  }
  G4double y1[BL3SIZE];
  G4double y2[BL3SIZE];
  G4double y3[BL3SIZE];
  for(G4int init_i = 0; init_i < BL1SIZE; init_i++) {
    y1[init_i] = 0.0;
    y2[init_i] = 0.0;
    y3[init_i] = 0.0;
  }
  G4double z = 0.0;
  G4double za_i = 0.0;
  G4double zai2 = 0.0;
  G4double zshif = 0.0;
  G4double ztouch = 0.0; 
  G4double ztu = 0.0;

  // LIEGE INC-model as a subroutine 

  // The Liege INC model has been applied to several systems and in
  // several conditions. Refinements are still in progress.

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
  G4double fmpinct[9] = {938.2796,938.2796,138.0,138.0,138.0,1874.35,2806.8,2806.8,3727.};

  // Initialize array:
   for(G4int i = 0; i < 300; i++) {
     npproj[i] = 0;
     isPartOfSpectatorNucleus[i] = false;
     nc[i] = 0;
   }

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
  //                  Kclst=F(10) Light clusters computed if not =0
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

  // for logging
  //  std::ofstream dumpOut("inclDump.txt");
  // end for logging

  G4int ia_breakup = 0;
  G4int jparticip[BL1SIZE];
  for(G4int i = 0; i < BL1SIZE; i++) {
    jparticip[i] = 0;
  }

  G4double beproj = 0.;
  bl3->ia2 = calincl->targetA(); // f(1)->f[0] and so on..., calincl added
  G4int iz2 = calincl->targetZ();
  G4double r02 = 1.12;
  kindstruct->kindf7 = calincl->bulletType();
  if(kindstruct->kindf7 < 0) kindstruct->kindf7 = kindstruct->kindf7-1;
  G4int kclst = calincl->getClusterOption();

  // G4cout <<"Projectile kind = " << kindstruct->kindf7 << G4endl;

  G4int iz1 = 0;
  G4double fmpinc = 0.0;
  G4double rms1 = 0.0;
  G4double pf1 = 0.0;
  G4double tlab = calincl->bulletE();;

  if(kindstruct->kindf7 > 0) {
    bl3->ia1 = ia1t[G4int(kindstruct->kindf7)-1]; // In the C++ version indices start from 0
    iz1 = iz1t[G4int(kindstruct->kindf7)-1];
    fmpinc = fmpinct[G4int(kindstruct->kindf7)-1];
    rms1 = light_gaus_nuc->rms1t[G4int(kindstruct->kindf7)-1];
    pf1 = light_gaus_nuc->pf1t[G4int(kindstruct->kindf7)-1];

  } else {
   // Starting values for carbon beams:
    if(kindstruct->kindf7 == -12) { // Handle Carbon-12 as a special case
      bl3->ia1 = std::abs(kindstruct->kindf7 + 1);
      // check here that kindf7 is -a and ia1 is a of the projectile.
      //  			write(6,*) 'kindf7,ia1:',kindf7,ia1
      be->ia_be=bl3->ia1;
      iz1=be->iz_be;
      fmpinc=bl3->ia1*fmp - be->bind_be;
      rms1=be->rms_be;
      pf1=be->pms_be;
      ia_breakup = 10; // Stopping threshold for the cascade in case of light nuclei
      //      ia_breakup=calincl->getCascadeStoppingAThreshold();
    } else { // For extended projectiles
      bl3->ia1 = calincl->extendedProjectileA();
      be->ia_be = calincl->extendedProjectileA();
      iz1=be->iz_be;
      fmpinc=bl3->ia1*fmp - be->bind_be;
      rms1=be->rms_be;
      pf1=be->pms_be;
      ia_breakup = 10; // Stopping threshold for the cascade in case of light nuclei
      //      ia_breakup=calincl->getCascadeStoppingAThreshold();
    }
  }

  if(verboseLevel > 1) {
    G4cout <<"rms1 = " << rms1 << G4endl;
    G4cout <<"pf1 = " << pf1 << G4endl;
  }

  G4int k2 = 0;
  G4int k3 = 0;
  //  G4int k3 = 1; // No deltas!
  G4int k4 = 2;
  G4int k5 = 1;
  G4int k6 = 0;

  // material number:      
  saxw->imat = 0;
  // espace de phases test (r et p) pour pauli: 
  // valeur recommandee par j.c. v-test=0.589 h**3:
  //  G4double rbl = 2.0;

  // Valeur pour avoir V-test=2.38 h**3 (avec pbl=200)
  G4double rbl=3.1848;
  G4double pbl=200.0;

  paul->xrr = rbl;
  paul->xrr2 = (paul->xrr) * (paul->xrr);
  paul->pr=pbl;
  paul->pr2 = paul->pr*(paul->pr);

  G4double tem[10];
  tem[0] = 100000.0;  // tem(1) -> tem[0]
  // temfin (time at which the inc is stopped), tmax5 defined after chosing b

  G4double v0 = calincl->getPotential();
  G4double v1 = v0;
  bl8->rathr = 0.;
  bl8->ramass = 0.;

  // constants and derived data
  bl10->pf = 1.37*hc;
  G4double tf = std::sqrt(bl10->pf*(bl10->pf)+fmp*fmp)-fmp;
  G4double g0 = 115.0;
  G4double th = 0.;
  G4double pm2 = fmp*fmp;
  G4int ia = bl3->ia1 + bl3->ia2;
  G4int a2 = bl3->ia2;
  bl3->r2 = r02*std::pow(G4double(a2),0.33333333);

  // parametres moyens de densite de la cible (fermi 2 parametres)
  if (bl3->ia2 >= 28) { //then
    ws->r0 = (2.745e-4*bl3->ia2+1.063)*std::pow(G4double(bl3->ia2),0.33333333);
    ws->adif = 1.63e-4*bl3->ia2 + 0.510;
    ws->rmaxws = ws->r0 + ws->xfoisa*(ws->adif);
  }
  else if(bl3->ia2 >= 19) { //then
    ws->r0 = light_nuc->r[bl3->ia2-1];
    ws->adif = light_nuc->a[bl3->ia2-1];
    ws->rmaxws = ws->r0 + ws->xfoisa*(ws->adif);
  }
  else if(bl3->ia2>=6) { //then
    ws->r0 = 1.581*(light_nuc->a[bl3->ia2-1]) * (2.0 + 5.0 * (light_nuc->r[bl3->ia2-1]))/(2.0 + 3.0*(light_nuc->r[bl3->ia2-1]));
    ws->adif = light_nuc->a[bl3->ia2-1];
    ws->rmaxws = 5.5 + 0.3*(double(bl3->ia2) - 6.0)/12.0;
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
  if(verboseLevel > 3) {
    G4cout <<"Radius bl3->r2 = " << bl3->r2 << G4endl;
  }

  G4double tnr = tlab;
  G4double binc = std::sqrt(tnr*tnr + 2.0*tlab*fmpinc)/(tnr+fmpinc);
  G4double ginc = 1.0/std::sqrt(1.0 - binc*binc);
  G4double pinc = fmpinc*binc*ginc;

  for(G4int bli = 0; bli < BL1SIZE; bli++) {
    bl1->eps[bli] = 0.0;
    bl1->p1[bli] = 0.0;
    bl1->p2[bli] = 0.0;
    bl1->p3[bli] = 0.0;
    q1[bli] = 0.0;
    q2[bli] = 0.0;
    q3[bli] = 0.0;
    q4[bli] = 0.0;
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
  if (kindstruct->kindf7 < 6 && kindstruct->kindf7 > 0)  { //then
    ws->bmax = ws->bmax;     // comme alain
  }
  else {
    beproj = fmpinc - bl3->ia1*fmp;
    ws->bmax = ws->rmaxws + rms1;     // maximum extension of the nucleus ( w.s.)
  }

  // deutons     
  G4double al;
  standardRandom(&al, &(hazard->ial));
  G4double b = std::sqrt(al)*(ws->bmax);
  G4double bred = b/bl3->r2;
  //G4double bimpact=b;
  bimpact = b;
  G4double tnor = 0.0;

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
    if (kindstruct->kindf7 <= 2 && kindstruct->kindf7 > 0) {
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
      temfin = 1.25*cb0/amax1(1.0,0.854 + 0.438*bred)*std::pow(G4double(bl3->ia2),(eb0/amax1(1.0,0.941+0.177*bred)));
      temfin = temfin*tnor;
    }
    else {
      if (kindstruct->kindf7 < 6 && kindstruct->kindf7 > 0) {
	// here for pions:
	temfin = 30.0*std::pow((float(bl3->ia2)/208.0),0.25)*(1.0 - 0.2*bred)*(1.0 - tlab/1250.0);
	// correction for pions in the case nosurf=0 or -1 (a.b., c.v. 2/2002)
	temfin = temfin*tnor;
      }
      else {
	// deutons
	tlabu = tlab/double(bl3->ia1);
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
	temfin = 1.25*(coeffb0/xc)*std::pow(G4double(bl3->ia2),(expob0/xe));
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
      //      temfin = 25.5*std::pow(G4double(bl3->ia2),0.16);  // pb208->60fm/c

      // c modified in June 2005 (function of A and TLAB)
      temfin = 30.18*std::pow(bl3->ia2,0.17)*(1.0-5.7*std::pow(10.0,(-5))*tlab);
    }
    else {
      // here for other hadrons
      temfin = 29.8*std::pow(G4double(bl3->ia2),0.16);  // pb208->70fm/c
    }
  }

  // deutons
  // a way to change stopping time f[5] not used here
  factemp = calincl->getTimeScale();
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
    //    bl5->nesc[0] = 0; // nesc(1)->nesc[0] and so on...
    bl5->nesc[1] = 0;
    //    bl1->ind2[0] = 2*iz1 - 1;
    bl1->ind2[1] = 2*iz1 - 1;
    bl1->ind1[1] = 0;
    bl3->x1[1] = 0.0;
    bl3->x2[1] = 0.0;
    bl3->x3[1] = 0.0;
    bl1->p1[1] = 0.0;
    bl1->p2[1] = 0.0;
    bl1->p3[1] = 0.0;
    bl9->hel[0] = 0.0;
    bl9->hel[1] = 0.0;
    jparticip[0] = 0;
    jparticip[1] = 1;
  }
  else {
    npion = 1;
    if(kindstruct->kindf7 == -666) { // Extended projectile...
      ipi[1] = 8 - 2*be->ia_be;
    } else {
      ipi[1] = 8 - 2*(kindstruct->kindf7);
    }
    y1[1] = 0.0;
    y2[1] = 0.0;
    y3[1] = 0.0;
    q1[1] = 0.0;
    q2[1] = 0.0;
    q3[1] = 0.0;
    q4[1] = fmpi;
  }

  // deutons
  goto pnu9;
  //    850	
 pnu7:
  //  G4cout <<"(checkpoint 'pnu7)" << G4endl;
  if(kindstruct->kindf7 == 6) { // Deuteron
    gaussianRandom(&xga);
    bl3->x1[1] = xga*rms1*0.5775;
    gaussianRandom(&xga);
    bl3->x2[1] = xga*rms1*0.5775;
    gaussianRandom(&xga);
    bl3->x3[1] = xga*rms1*0.5775;

    bl3->x1[2] = -bl3->x1[1];
    bl3->x2[2] = -bl3->x2[1];
    bl3->x3[2] = -bl3->x3[1];

    // Deuteron density from Paris potential in q space:
    standardRandom(&xq, &(hazard->igraine[9]));
    qdeut = splineab(xq) * 197.3289;
    standardRandom(&u, &(hazard->igraine[11]));
    cstet = u * 2 - 1;
    sitet = std::sqrt(1.0 - std::pow(cstet, 2));
    standardRandom(&v, &(hazard->igraine[11]));
    phi = 2.0 * 3.141592654 * v;
    bl1->p1[1] = qdeut * sitet * std::cos(phi);
    bl1->p2[1] = qdeut * sitet * std::sin(phi);
    bl1->p3[1] = qdeut * cstet;
    bl1->eps[1] = w(bl1->p1[1], bl1->p2[1], bl1->p3[1], fmp);
    bl1->p1[2] = - bl1->p1[1];
    bl1->p2[2] = - bl1->p2[1];
    bl1->p3[2] = - bl1->p3[1]; 
    bl1->eps[2] = bl1->eps[1];
    //    bl1->eps[2] = w(bl1->p1[2], bl1->p2[2], bl1->p3[2], fmp);

    jparticip[1] = 1;
    bl9->hel[1] = 0;
    bl5->nesc[1] = 0;
    bl1->ind2[1] = 1;
    bl1->ind1[1] = 0;
    jparticip[2] = 1;
    bl9->hel[2] = 0;
    bl5->nesc[2] = 0;
    bl1->ind2[2] = -1;
    bl1->ind1[2] = 0;
  } else {
    // Composite heavier than a deuteron
    s1t1 = 0.0;
    s2t1 = 0.0;
    s3t1 = 0.0;
    sp1t1 = 0.0;
    sp2t1 = 0.0;
    sp3t1 = 0.0;
    G4double sei = 0.0;

    bl9->hel[0] = 0;
    bl5->nesc[0] = 0;
    for(G4int i = 1; i <= bl3->ia1; i++) {
      bl9->hel[i] = 0;
      bl5->nesc[i] = 0;
      bl1->ind2[i] = 1;
      if(i > iz1) {
	bl1->ind2[i] = -1;
      }
      bl1->ind1[i] = 0;
      jparticip[i] = 1;
      gaussianRandom(&xga);
      bl3->x1[i] = xga * rms1 * 0.57735;
      s1t1 = s1t1 + bl3->x1[i];
      gaussianRandom(&xga);
      bl3->x2[i] = xga * rms1 * 0.57735;
      s2t1 = s2t1 + bl3->x2[i];
      gaussianRandom(&xga);
      bl3->x3[i] = xga * rms1 * 0.57735;
      s3t1 = s3t1 + bl3->x3[i];

      G4int igaussfermi = 2;
      if(igaussfermi == 1) { // Gaussian
	// Density of composite as a gaussian in q
	// space:
	gaussianRandom(&xga);
	bl1->p1[i] = xga * pf1 * 0.57735;
	gaussianRandom(&xga);
	bl1->p2[i] = xga * pf1 * 0.57735;
	gaussianRandom(&xga);
	bl1->p3[i] = xga * pf1 * 0.57735;
      } else { // Fermi
	// Density of composite as a Fermi sphere in q space:
	for(int iloc=1; iloc <= 3; iloc++) {
	  standardRandom(&t[iloc],&(hazard->igraine[11]));
	}
	t[2]=-1+2.*t[2];
	t[3]=6.283185*t[3];
	t2=t[2];
	t3=std::sqrt(1.-t2*t2);
	t4=std::cos(t[3]);
	t5=std::sin(t[3]);
	y=fermi->pf*std::pow(t[1],0.33333333);
	bl1->p1[i]=y*t3*t4;
	bl1->p2[i]=y*t3*t5;
	bl1->p3[i]=y*t2;
      }
      bl1->eps[i] = w(bl1->p1[i], bl1->p2[i], bl1->p3[i], fmp);
      sei = sei + bl1->eps[i];

      if(verboseLevel > 1) {
	G4cout <<"projectile nucleon = " << i << G4endl;
	G4cout <<"x1 = " << bl3->x1[i] << G4endl;
	G4cout <<"x2 = " << bl3->x2[i] << G4endl;
	G4cout <<"x3 = " << bl3->x3[i] << G4endl;
	G4cout <<"p1 = " << bl1->p1[i] << G4endl;
	G4cout <<"p2 = " << bl1->p2[i] << G4endl;
	G4cout <<"p3 = " << bl1->p3[i] << G4endl;
	G4cout <<"eps = " << bl1->eps[i] << G4endl;
      }
      sp1t1 = sp1t1 + bl1->p1[i];
      sp2t1 = sp2t1 + bl1->p2[i];
      sp3t1 = sp3t1 + bl1->p3[i];
    }

    s1t1 = s1t1/double(bl3->ia1);
    s2t1 = s2t1/double(bl3->ia1);
    s3t1 = s3t1/double(bl3->ia1);
    sp1t1 = sp1t1/double(bl3->ia1);
    sp2t1 = sp2t1/double(bl3->ia1);
    sp3t1 = sp3t1/double(bl3->ia1);
    sei = 0.0;

    for(G4int i = 1; i <= bl3->ia1; i++) {
      bl3->x1[i] = bl3->x1[i] - s1t1;
      bl3->x2[i] = bl3->x2[i] - s2t1;
      bl3->x3[i] = bl3->x3[i] - s3t1;
      bl1->p1[i] = bl1->p1[i] - sp1t1;
      bl1->p2[i] = bl1->p2[i] - sp2t1;
      bl1->p3[i] = bl1->p3[i] - sp3t1;
      bl1->eps[i] = w(bl1->p1[i], bl1->p2[i], bl1->p3[i], fmp);
      sei = sei + bl1->eps[i];
    }
    v_proj = 0.0;
    if(bl3->ia1 > 4) v_proj = (sei - fmpinc)/bl3->ia1;

    //    bl1->dump(26);
    //    bl3->dump();
    //    G4cout <<"(checkpoint 'endheavyprojinit)" << G4endl;
  } // End of: Composites heavier than deuteron


 pnu9: // continue
  //  G4cout <<"(checkpoint 'pnu9)" << G4endl;
  // deutons
  // target preparation for 1 < a < 5 (with sum of momentum =0)
  if(bl3->ia2 >= 2 && bl3->ia2 <= 4) {
  pnu1633:  
    s1t1 = 0.0;
    s2t1 = 0.0;
    s3t1 = 0.0;
    sp1t1 = 0.0;
    sp2t1 = 0.0;
    sp3t1 = 0.0;
    efer = 0.0;
    for(G4int i = bl3->ia1+1; i <= ia; i++) {
    //    for(G4int i = bl3->ia1; i < ia-1; i++) {
    //    for(G4int i = bl3->ia1; i < ia; i++) {
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
      t2 = std::sqrt(1.0 - t1*t1);                                       
      t3 = std::cos(t[2]);  //t(3)->t[2]                                            
      t4 = std::sin(t[2]);   //t(3)->t[2]                                                                                      
      t5 = t[4];       // t(5)->t[4]                                                 
      t6 = std::sqrt(1.0 - t5*t5);                                         
      t7 = std::cos(t[5]);   //t(6) -> t[5]                                             
      t8 = std::sin(t[5]);   // t(6)->t[5]                                          
      if (ws->nosurf == 1) {
	x = bl3->r2*std::pow(t[0],0.33333333); // t(1)->t[0]                                      
	y = (bl10->pf)*std::pow(t[3],0.33333333); // t(4)->t[3]  
      }
      else {
	// surface..w.s.: impulsion (sphere dure), puis r(q)
	t33 = std::pow(t[6],0.33333333); // t(7)->t[6]
	y = (bl10->pf)*t33;

	rr = interpolateFunction(t33);
	x = rr*std::pow(t[3],0.33333333); // t(4)->t[3] 
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

    p_mod = std::sqrt(std::pow(bl1->p1[ia],2) + std::pow(bl1->p2[ia],2) + std::pow(bl1->p3[ia],2));
    if(p_mod > ((bl10->pf)+0.05)) {
      goto pnu1633;
    }

    bl1->eps[ia] = w(bl1->p1[ia],bl1->p2[ia],bl1->p3[ia],fmp);
    
    efer = efer + bl1->eps[ia]-fmp;                                        

  } //end if       !(bl3->ia2 >= 2 && bl3->ia2 <= 4)

  // target preparation for a > 4
  if(bl3->ia2 > 4) {
    x1_target = 0.0;
    x2_target = 0.0;
    x3_target = 0.0;
    for(G4int i = bl3->ia1+1; i <= ia; i++) { //do 1 i=bl3->ia1+1,ia
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
      t2 = std::sqrt(1.0 - t1*t1);
      t3 = std::cos(t[2]);  //t(3)->t[2]
      t4 = std::sin(t[2]); //t(3)->t[2]
      t5 = t[4]; //t(5)->t[4]
      t6 = std::sqrt(1.0 - t5*t5);
      t7 = std::cos(t[5]);  //t(6)->t[5]
      t8 = std::sin(t[5]);   // t(6)->t[5]

      if (ws->nosurf == 1) {
	x = bl3->r2*std::pow(t[0],0.33333333); // t(1)->t[0]
	y = bl10->pf*std::pow(t[3],0.33333333); // t(4)->t3
      }
      else {
	// surface..w.s.: impulsion (sphere dure), puis r(q)
	t33 = std::pow(t[6],0.33333333); // t(7)->t[6]
	y=bl10->pf*t33;
	rr=interpolateFunction(t33);
	x=rr*std::pow(t[3],0.33333333); // t(4)->t[3]       
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
    x1_target = x1_target/double(bl3->ia2);
    x2_target = x2_target/double(bl3->ia2);
    x3_target = x3_target/double(bl3->ia2);
  }

  efrun = efrun + efer;

  // location of incident particle at point (b,z)
  r22 = bl3->r2*(bl3->r2);
  z = ws->bmax * (ws->bmax) - b*b;       // for the wood-saxon density...

  if (z < 0.0) {
    z=0.0;
  }

  z = std::sqrt(z);
  // random azimuthal direction of the impact parameter (sept 99)

  if (kindstruct->kindf7 <= 2 && kindstruct->kindf7 > 0) {
    standardRandom(&tbid, &(hazard->igraine[13]));
    tbid = tbid*6.283185;
    bl3->x1[1] = bl3->x1[1] + b*std::cos(tbid);  //x1(1)->x1[1]                                        
    bl3->x2[1] = bl3->x2[1] + b*std::sin(tbid); //x2(1)->x2[1]
    bl3->x3[1] = bl3->x3[1] - z;
    // Counter of p(n) entering the potential 4/2008 AB
    if(bl1->ind2[1] == 1) { // then
      npenter=1;
      nnenter=0;
    } else {
      npenter=0;
      nnenter=1;
    } // endif

    // pour le ntuple des avatars:
    if(varavat->kveux == 1) {
      varavat->r1_in[0] = bl3->x1[1]; //r1_in(1)->r1_in[0] and x1(1)->x1[0]
      varavat->r1_in[1] = bl3->x2[1]; //r1_in(2)->r1_in[1] and x1(2)->x1[1]
      varavat->r1_in[2] = bl3->x3[1]; //r1_in(3)->r1_in[2] and x1(3)->x1[2]
    } //endif
  }
  else {
    if (kindstruct->kindf7 < 6 && kindstruct->kindf7 > 0) { //then ! pour les pions on laisse
      //call standardRandom(tbid,iy14)
      standardRandom(&tbid, &(hazard->igraine[13]));
      tbid = tbid*6.283185;
      y1[1] = y1[1] + b*std::cos(tbid); //y1(1)->y1[0]
      y2[1] = y2[1] + b*std::sin(tbid); //y2(1)->y2[0]
      y3[1] = y3[1] - z;
    }
    else {
      // deutons
      //nmiss=0.;
      nmiss = 0;
      xlengm=1000.0;
      
      for(G4int i = 1; i <= bl3->ia1; i++) {
      //      for(G4int i = 1; i < bl3->ia1; i++) {
	bl3->x3[i] = bl3->x3[i]/ginc;
	zai2 = ws->rmaxws*(ws->rmaxws) - std::pow((b+bl3->x1[i]),2) - std::pow(bl3->x2[i],2);
	if (zai2 < 0.0) {
	  goto pnu22;
	}
	ztu = -std::sqrt(zai2);
	// r22 remplace par rmaxws*rmaxws et r2 par rmaxws cv correct ?
	za_i = 2.0*(ws->rmaxws) + ztu;
	xleng = za_i - bl3->x3[i];
	if (xleng > xlengm) {
	  //  goto pnu21;
	  continue;
	}
	ilm = i;
	xlengm = xleng;
	ztouch = ztu;
	//	goto pnu21;
	continue;
      pnu22:
	nmiss = nmiss + 1; 
      }
      //    pnu21:
      if (nmiss == bl3->ia1) { //then
	//	G4cout <<"(checkpoint 'projectilemissedtarget)" << G4endl;
	nopart = -1;
	//	return;
	if(verboseLevel > 3) {
	  G4cout <<"nmiss == bl3->ia1" << G4endl;
	  G4cout <<"End of algorithm after pnu21." << G4endl;
	}
	goto pnureturn;
      }
      else {
	// Counter of p(n) entering the potential 4/2008 AB
        if(bl1->ind2[ilm]==1) { //then
	  npenter=1;
	  nnenter=0;
	} else {
	  npenter=0;
	  nnenter=1;
	} // endif

	zshif = bl3->x3[ilm] - ztouch;
	standardRandom(&tbid, &(hazard->igraine[13]));
	tbid = tbid*6.283185;
	for(G4int i = 1; i <= bl3->ia1; i++) {
	  xxx = bl3->x1[i] + b;
	  bl3->x1[i] = xxx*std::cos(tbid) - bl3->x2[i]*std::sin(tbid);
	  bl3->x2[i] = xxx*std::sin(tbid) + bl3->x2[i]*std::cos(tbid);
	  bl3->x3[i] = bl3->x3[i] - zshif;
	}
	if (std::fabs(std::pow(bl3->x1[ilm],2)+std::pow(bl3->x2[ilm],2)+std::pow(bl3->x3[ilm],2)-ws->rmaxws*(ws->rmaxws)) > 0.01) {
	  if(verboseLevel > 2) {
	    G4cout <<"wrong position" << G4endl;
	  }
	}
      }
    }
  }

  // for rho(r),rho(q) checking
#ifdef G4INCLDEBUG
  for(G4int rho_i = bl3->ia1 + 1; rho_i <= ia; rho_i++) { // DO I=IA1+1,IA
    G4double r_dist = std::sqrt(bl3->x1[rho_i] * bl3->x1[rho_i] + bl3->x2[rho_i] * bl3->x2[rho_i] + bl3->x3[rho_i] * bl3->x3[rho_i]);
    G4double p_mod = std::sqrt(std::pow(bl1->p1[rho_i],2)+std::pow(bl1->p2[rho_i],2) + std::pow(bl1->p3[rho_i],2));
    theLogger->fillHistogram1D("r_distrib", r_dist);
    theLogger->fillHistogram1D("p_distrib", p_mod);
    theLogger->fillHistogram2D("r-p_correl", r_dist, p_mod);
  }
#endif

  // initial momentum for all type of incident particles:
  xl1 = b*pinc*std::sin(tbid);                                  
  xl2 = -b*pinc*std::cos(tbid);                                          
  xl3 = 0.0;                                           

  // transcription in the general frame of reference
  // (here,=lab frame)

  bevar = 0.0;
  ge = 1.0;
  b1 = (binc - bevar)/(1.0 - bevar*binc);
  b2 = -bevar;
  g1 = 1.0/std::sqrt(1.0 - b1*b1);
  g2 = 1.0;
  // deutons
  // here for nucleons
  if (kindstruct->kindf7 <= 2 && kindstruct->kindf7 > 0) {
    bl1->eps[1] = g1*fmp + v0;
    bl1->p3[1] = std::sqrt(std::pow(bl1->eps[1],2) - std::pow(fmp,2));
  }
  else {
    // here for pions
    if (kindstruct->kindf7 < 6 && kindstruct->kindf7 > 0) { //then
      q4[1] = g1*fmpi; // q4(1)->q4[0]
      q3[1] = b1*q4[1];
    }
    else { // Composite projectiles
      // here for composite projectiles:
      // the kinetic energy is below the threshold. put all
      // fermi momentum to 0... projectile nucleons not on shell!
      energie_in = tlab + fmpinc;
      if((energie_in) <= (bl3->ia1*fmp)) {
	for(G4int i = 1; i <= bl3->ia1; i++) {
	  bl1->eps[i] = energie_in/double(bl3->ia1);
	  bl1->p1[i] = 0.0;
	  bl1->p2[i] = 0.0;
	  bl1->p3[i] = pinc/double(bl3->ia1);
	}
	goto pnu1871;
      }
      // here the composit is above threshold
      for(G4int i = 1; i <= bl3->ia1; i++) { //do i=1,bl3->ia1	!save e,p in the composit rest frame
	qvp->t_c[i] = bl1->eps[i] - fmp;
	bl1->eps[i] = bl1->eps[i] - v_proj;
	qvp->eps_c[i] = bl1->eps[i];
	qvp->p3_c[i] = bl1->p3[i];
      } //enddo

      nbtest = bl3->ia1 - 1;
      if(kindstruct->kindf7 == 6) {
	nbtest = 2;
      }
      iflag = 0;
    pnu1870:
      //      G4cout <<"(checkpoint 'pnu1870)" << G4endl;
      sueps = 0.0;

      iflag = iflag + 1;

      for(G4int i = 1; i <= bl3->ia1; i++) { //do i=1,bl3->ia1
	tte = qvp->eps_c[i];
	bl1->eps[i] = g1*(qvp->eps_c[i] + b1*qvp->p3_c[i]);
	bl1->p3[i] = g1*(b1*tte + qvp->p3_c[i]);
	sueps = sueps + bl1->eps[i];
      } //enddo

      goto pnu987;
      // With a potential for the projectile it is nonsens to search for 
      // nucleons ON shell  31/05/2010.      

      cobe = (tlab + fmpinc)/sueps;

      // off shell problem for incident clusters (a.b. 2/2002)

      if(iflag == nbtest) { // too much..all momentum to 0
	for(G4int klm = 1; klm <= bl3->ia1; klm++) { //do klm=1,bl3->ia1
	  qvp->eps_c[klm] = fmp;
	  bl1->p1[klm] = 0.0;
	  bl1->p2[klm] = 0.0;
	  qvp->p3_c[klm] = 0;
	}
	goto pnu1870;
      }
      for(G4int i = 1; i <= bl3->ia1; i++) { //do i=1,bl3->ia1
	arg = std::pow((cobe*(bl1->eps[i])),2)-pm2;
	if (arg <= 0.) { //then	! put maximum momentum to 0. 
	  i_emax = 1; //	!find maximum
	  ener_max = bl1->eps[1];
	  for(G4int klm = 2; klm <= bl3->ia1; klm++) { //do klm=2,bl3->ia1	
	    if(bl1->eps[klm] > ener_max) {
	      ener_max = bl1->eps[klm];
	      i_emax = klm;
	    }
	  }
	  qvp->eps_c[i_emax] = fmp;
	  bl1->p1[i_emax] = 0.0;
	  bl1->p2[i_emax] = 0.0;
	  qvp->p3_c[i_emax] = 0.0;

	  if(i_emax == bl3->ia1) { //    circular permut if the last one
	    epsv = qvp->eps_c[bl3->ia1]; //	 permutation circulaire
	    p1v = bl1->p1[bl3->ia1];
	    p2v = bl1->p2[bl3->ia1];
	    p3v = qvp->p3_c[bl3->ia1];
	    for(bl2->k = bl3->ia1-1; bl2->k >= 1; bl2->k = bl2->k - 1) { //do k=bl3->ia1-1,1,-1
	      qvp->eps_c[bl2->k+1] = qvp->eps_c[bl2->k];
	      bl1->p1[bl2->k+1] = bl1->p1[bl2->k];
	      bl1->p2[bl2->k+1] = bl1->p2[bl2->k];
	      qvp->p3_c[bl2->k+1] = qvp->p3_c[bl2->k];
	    }
	    qvp->eps_c[1] = epsv;
	    bl1->p1[1] = p1v;
	    bl1->p2[1] = p2v;
	    qvp->p3_c[1] = p3v; 	// fin permut.
	  }
	  sp1t1 = 0.0;   // re-compute the last one 
	  sp2t1 = 0.0;
	  sp3t1 = 0.0;
	  for(G4int j = 1; j <= bl3->ia1-1; j++) { //do j=1,bl3->ia1-1
	    sp1t1 = sp1t1 + bl1->p1[j];
	    sp2t1 = sp2t1 + bl1->p2[j];
	    sp3t1 = sp3t1 + qvp->p3_c[j];
	  }
	  bl1->p1[bl3->ia1] = -sp1t1;
	  bl1->p2[bl3->ia1] = -sp2t1;
	  qvp->p3_c[bl3->ia1] = -sp3t1;
	  qvp->eps_c[bl3->ia1] = w(bl1->p1[bl3->ia1],bl1->p2[bl3->ia1],qvp->p3_c[bl3->ia1],fmp);	

	  goto pnu1870;  // ..and boost all of them.
	}
      }

      for(G4int i = 1; i <= bl3->ia1; i++) { //do i=1,bl3->ia1
      //      for(G4int i = 1; i < bl3->ia1; i++) { //do i=1,bl3->ia1
	arg = std::pow((cobe*(bl1->eps[i])),2) - pm2;
	comom = std::sqrt(arg/(std::pow(bl1->eps[i],2) - pm2));
	bl1->p1[i] = comom*(bl1->p1[i]);
	bl1->p2[i] = comom*(bl1->p2[i]);
	bl1->p3[i] = comom*(bl1->p3[i]);
	bl1->eps[i] = bl1->eps[i]*cobe;
	if (std::fabs(am(bl1->p1[i],bl1->p2[i],bl1->p3[i],bl1->eps[i])-fmp) > 0.01) {
	  if(verboseLevel > 2) {
	    G4cout <<"wrong correction " << i << G4endl;                  
	  }
	}
      }
      bl1->eps[ilm] = bl1->eps[ilm] + v0;  

 pnu987:
      //  G4cout <<"(checkpoint 'pnu987)" << G4endl;
  // ici sauvegarde des energies et impulsions des nucleons du projectile i.l.:
  // (hors potentiel)
  for(G4int i = 1; i <= bl3->ia1; i++) {
    qvp->eps_c[i]=bl1->eps[i];
    qvp->p1_s[i]=bl1->p1[i];
    qvp->p2_s[i]=bl1->p2[i];
    qvp->p3_s[i]=bl1->p3[i];
  }

  // first particle entering in the target potential:	
  bl1->eps[ilm]=bl1->eps[ilm]+v0;
	  
  // correction en test 16/04/2009 (pourquoi pas la correction aussi de l'impulsion
  // comme pour un nucleon in ???? ....on essaye de traiter les nucleons hors couche?)
  // 29/4/2009: je pense que c'est une erreur et je transforme aussi maintenant l'impulsion:
  var_ab = std::pow(bl1->p1[ilm],2)
    + std::pow(bl1->p2[ilm],2)
    + std::pow(bl1->p3[ilm],2);
  if(var_ab > 0.) {
    gpsg=std::sqrt((std::pow(bl1->eps[ilm],2)
		    -fmp*fmp)/var_ab);
  }
  bl1->p1[ilm]=bl1->p1[ilm]*gpsg;             
  bl1->p2[ilm]=bl1->p2[ilm]*gpsg;
  bl1->p3[ilm]=bl1->p3[ilm]*gpsg;        
    } // For pion-in
  } // For nucleon-in

 pnu1871:
  //  G4cout <<"(checkpoint 'pnu1871)" << G4endl;
  //  bl1->dump(27);
  //  bl3->dump();
  //  randomGenerator->printSeeds();
  // evaluation of the times t(a,b)
  bl2->k = 0;
  kcol = 0;
  if (kindstruct->kindf7 <= 2 && kindstruct->kindf7 > 0) {
    // modif s.vuillier tient compte propagation projectile,1e collision
    // imposee pour lui (c'est une maniere de faire!!)
    G4int ioldk = 0;
    for(G4int i = 1; i <= ia; i++) { //do 40 i=1,ia
      ioldk = bl2->k;
      //   bl2->k = bl2->k + 1;
      //		tref=ref(x1(i),x2(i),x3(i),p1(i),p2(i),p3(i),eps(i),r22)          p-n04740
      tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
      if (tref > bl4->tmax5) {
	goto pnu45;
      }

      bl2->k = bl2->k + 1;
      bl2->crois[bl2->k]=tref;
      bl2->ind[bl2->k]=i;
      bl2->jnd[bl2->k]=-1;
    pnu45:
      i1=i-1;
      if (i == 1) {
	//goto pnu40;
	continue;
      }
      //   1326	c ici on ne calcule que les G4interactions nn impliquant le projectile !!! (2/02)
      if (i1 > bl3->ia1) {
	i1=bl3->ia1;
      }
      for(G4int j = 1; j <= i1; j++) { //do 41 j=1,i1
	// no collisions before the first collision of the incident particle
	time (i, j);
	if (bl1->ta < 0.0) {
	  continue;
	}
	if(bl1->ta > bl4->tmax5) {
	  continue;
	}
	eij=am(bl1->p1[i]+bl1->p1[j],bl1->p2[i]+bl1->p2[j],bl1->p3[i]+bl1->p3[j],bl1->eps[i]+bl1->eps[j]);
	if (eij < 1925.0) {
	  continue;
	}
	isos=bl1->ind2[i]+bl1->ind2[j];

	if (31.0*(bl3->rab2) > totalCrossSection(eij,0,isos)) {
	  continue;
	}

	bl2->k = bl2->k + 1;
	if (j == 1) {
	  kcol = kcol + 1;
	}
	bl2->crois[bl2->k]=bl1->ta;
	bl2->ind[bl2->k]=i;
	bl2->jnd[bl2->k]=j;
      }
      if(verboseLevel > 3) {
	if(bl2->k == (ioldk + 2)) {
	  if(verboseLevel > 2) {
	    G4cout <<"bl2->k incremented twice!" << G4endl;
	  }
	}
      }
    }
  }
  else {
    // deutons
    if (kindstruct->kindf7 < 6 && kindstruct->kindf7 > 0) { //then
      // here for incoming pions:
      for(G4int i = bl3->ia1+1; i <= ia; i++) { //do i=ia1+1,ia
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); 
	if (tref < bl4->tmax5) {
	  bl2->k = bl2->k + 1;
	  bl2->crois[bl2->k] = tref;
	  bl2->ind[bl2->k] = i;
	  bl2->jnd[bl2->k] = -1;
	}
      }
      new2(y1[1], y2[1], y3[1], q1[1], q2[1], q3[1], q4[1], 1, 0);

      //   modif a.b. 21/06/2002: should check at least one valid collision
      //   1361	c    with incoming pion.
      //   1362	c      kcol=1
      if(bl2->k != 0) {
	kcol = 1;
      }
    }
    else {
      // Counting the spectators and transparents of composit projectiles:      
      //      G4cout <<"(checkpoint 'setipszero1)" << G4endl;
      ips=0;
      nb_transprojo=0;
      for(G4int i = 1; i <= bl3->ia1; i++) { //do 38 i=1,ia1
	bl5->nesc[i] = 1;
	// Spectators of composite projectiles (7/2006, AB)
	// Pour garder trace des spectateurs du projectile:
	npproj[i]=1;

	if (i != ilm) { 
	  goto pnu36;
	}
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22);
	if(verboseLevel > 3) {
	  if(tref < 0.0) {
	    if(verboseLevel > 2) {
	      G4cout <<"G4Incl: Reflection time < 0! (line 2579)" << G4endl;
	    }
	  }
	}
	bl5->nesc[i] = 0;
	npproj[i] = 0;
	goto pnu37;
      pnu36:
	//	G4cout <<"(checkpoint 'pnu36)" << G4endl;
	t1 = bl3->x1[i]*(bl1->p1[i])+bl3->x2[i]*(bl1->p2[i])+bl3->x3[i]*(bl1->p3[i]);                      
	t2 = bl1->p1[i]*(bl1->p1[i])+bl1->p2[i]*(bl1->p2[i])+bl1->p3[i]*(bl1->p3[i]);               
	t3 = t1/t2;
	t4 = bl3->x1[i]*(bl3->x1[i])+bl3->x2[i]*(bl3->x2[i])+bl3->x3[i]*(bl3->x3[i]); 
	//   1379	c incoming nucleons enter potential at maximum radius (modif. 13/06/01)
	t5 = t3*t3 + ((ws->rmaxws)*(ws->rmaxws) - t4)/t2;
	if (t5 < 0.) {
	  // this is a projectile spectator:
	  //	  G4cout <<"(checkpoint 'incripspnu36)" << G4endl;
	  ips=ips+1;
	  isPartOfSpectatorNucleus[i] = true;
	  ps->n_projspec[ips] = i;
	  continue;
	  //goto pnu38;
	} // endif
              	
	if(verboseLevel > 3) {
	  G4cout <<"x1 = " << bl3->x1[i] <<" x2 = " << bl3->x2[i] <<" x3 = " << bl3->x3[i] << G4endl; 
	  G4cout <<"t1 = " << t1 << G4endl;
	  G4cout <<"t2 = " << t2 << G4endl;
	  G4cout <<"t3 = " << t3 << G4endl;
	  G4cout <<"t4 = " << t4 << G4endl;
	  G4cout <<"rmaxws = " << ws->rmaxws << G4endl;
	  G4cout <<"t5 = " << t5 << G4endl;
	}
	//	if (t5 < 0.) {
	//	  ips++;
	//	  ps->n_projspec[ips] = i;
	//	  continue;
	//	}
	tref = (-1.0*t3 - std::sqrt(t5))*(bl1->eps[i]);  

	// Drop reflection avatars with too large tref and the ones
	// with negative tref (spectator, outside the nucleus)
	//	if (tref > bl4->tmax5 || tref <= 0.0) {
	//	  continue;
	//	}
	if (tref > bl4->tmax5) { //then
	  // This is also a projectile spectator:
	  ips=ips+1;
	  ps->n_projspec[ips] = i;
	  continue; // go to 38
	} // endif

	// Test to avoid a nucleon entering in the past (3/06/2010 as in incl4.6 a.b.)
	if(tref <= 0.0) { // then
	  // This is also a projectile spectator:
	  ips=ips+1;
	  ps->n_projspec[ips] = i;
	  continue; // go to 38
	} // endif
      
	npproj[i] = 1;
      pnu37:
	//	G4cout <<"(checkpoint 'pnu37)" << G4endl;
	bl2->k = bl2->k + 1; 
	bl2->crois[bl2->k] = tref; 
	bl2->ind[bl2->k] = i; 
	bl2->jnd[bl2->k] = -1;
      //    pnu38:
      }

      // For composit beams, cannot decide that there is no interaction,
      // (another nucleon can enter later at RMAX).
      kcol = 1;
      n_activnuc = ia-ips;  //number of active nucleons (not missing the target)
      n_enter_pot = ia-ips; //number of projectile nucleons that can enter

      /*
      if(useProjSpect) { // Treat the projectile spectators:
	n_activnuc = ia - ips; // number of active nucleons (not missing the target)
	if(ips != 0) { //then
	  e_spec = 0.0;
	  p1_spec = 0.0;
	  p2_spec = 0.0;
	  p3_spec = 0.0;
	  z_projspec = 0;
	  for(G4int i_spec = 1; i_spec <= ips; i_spec++) { // do i_spec=1,ips
	    G4int i_c = n_projspec[i_spec];
	    e_spec = e_spec + bl1->eps[i_c];
	    p1_spec = p1_spec + bl1->p1[i_c];
	    p2_spec = p2_spec + bl1->p2[i_c];
	    p3_spec = p3_spec + bl1->p3[i_c];
	    if(bl1->ind2[i_c] == 1) z_projspec = z_projspec + 1;
	  }
	  G4double p_spec2 = std::pow(p1_spec,2) + std::pow(p2_spec,2) + std::pow(p3_spec,2);
	  G4double s_spec = std::sqrt(std::pow(e_spec,2) - p_spec2);
	  // write(6,*) 'e,p,sqs:',e_spec,sqrt(p_spec2),s_spec

	  // c masses from dresner code:
	  G4int fja = ips;                                                                   
	  G4int fjz = z_projspec;                                                                   
	  G4double ex2 = abla->eflmac(fja, fjz, 0, 0); // dnrgy(fja,fjz);
	  G4double rnmass = ex2+um*fja;

	  // c write(6,*) fja,fjz,ex2,um                                                         
	  ex_projspec = s_spec - rnmass;
	  if(ex_projspec <= 0.0) { //then
	    // c no spectator nucleus compatible with available energy:
	    ips = 0;
	    // c n_activnuc = ia-ips!number of active nucleons (not missing the target)
	    goto pnu239;
	  }
	  t_projspec = e_spec - s_spec;
	  G4double coef = std::sqrt(t_projspec*(t_projspec + 2.0*s_spec)/p_spec2);
	  p1_projspec = p1_spec*coef;
	  p2_projspec = p2_spec*coef;
	  p3_projspec = p3_spec*coef;
	  m_projspec = rnmass;

	  // c check module...      
	  // 	c      write(6,*) 'projectile spectator a,z:',a_projspec,z_projspec
	  // 	c      write(6,*) 'm,t,e*:',m_projspec,t_projspec,ex_projspec
	  // 	c      write(6,*) 'p1,p2,p3:',p1_projspec,p2_projspec,p3_projspec
	  // 	c      write(6,*) '   constituants:',ips
	  // 	c      do i_spec=1,ips
	  // 	c        i_c = n_projspec(i_spec)
	  // 	c      write(6,*) i_c,eps(i_c),p1(i_c),p2(i_c),p3(i_c)       
	  // c      enddo
	  // c ...end check module 
	} // endif
      pnu239:
	a_projspec = ips; // number of nucleons in the spectator nucleus
      } // End of projectile spectator treatment
      */

      //      for(G4int i = bl3->ia1+1; i < ia; i++) { //do  39 i=ia1+1,ia
      for(G4int i = bl3->ia1+1; i <= ia; i++) { //do  39 i=ia1+1,ia
	npproj[i] = 0;
	tref = ref(bl3->x1[i], bl3->x2[i], bl3->x3[i], bl1->p1[i], bl1->p2[i], bl1->p3[i], bl1->eps[i], r22); // line 2609 
	if(verboseLevel > 3) {
	  if(tref < 0.0) {
	    G4cout <<"G4Incl: Reflection time < 0! (line 2609)" << G4endl;
	  }
	}

	if (tref < bl4->tmax5) { //then
	  bl2->k = bl2->k + 1;
	  bl2->crois[bl2->k]=tref;
	  bl2->ind[bl2->k]=i;
	  bl2->jnd[bl2->k]=-1;
	} //endif

	time (i, ilm);
	if (bl1->ta < 0.) {
	  continue;
	}
	if (bl1->ta > bl4->tmax5) {
	  continue;
	}
	eij=am(bl1->p1[i]+bl1->p1[ilm],bl1->p2[i]+bl1->p2[ilm],bl1->p3[i]+bl1->p3[ilm],bl1->eps[i]+bl1->eps[ilm]);
	if (eij < 1925.0) {
	  continue;
	}
	isos=bl1->ind2[i]+bl1->ind2[ilm];                                             
	if (31.*(bl3->rab2) > totalCrossSection(eij,0,isos)) {
	  continue;
	}								
	bl2->k = bl2->k + 1;                                                             
	kcol=kcol+1; 
	bl2->crois[bl2->k]=bl1->ta;                                                      
	bl2->ind[bl2->k]=i; 
	bl2->jnd[bl2->k]=ilm;
      }
    }
  }

  //  dumpBl3(dumpOut);

  if(verboseLevel > 3) {
    G4cout <<"Variables after time evaluation:" << G4endl;

    G4cout <<"bl3->ia1 + bl3->ia2 = ia " << G4endl;
    G4cout << bl3->ia1 << "         " << bl3->ia2 << "          " << ia << G4endl; 
    G4cout <<"B11" << G4endl;

    for(G4int idebugtest = 0; idebugtest <= bl2->k; idebugtest++) {
      G4cout <<"index = " << idebugtest << " ind1 = " << bl1->ind1[idebugtest] << " ind2 = " << bl1->ind2[idebugtest] << " p1 = " << bl1->p1[idebugtest] << " p2 = " << bl1->p2[idebugtest] << " p3 = " << bl1->p3[idebugtest] << " bl1->eps = " << bl1->eps[idebugtest] << G4endl;
    }
    
    G4cout <<"Bl2" << G4endl;
    for(G4int idebugtest = 0; idebugtest <= bl2->k; idebugtest++) {
      G4cout <<"index = " << idebugtest << " ind = " << bl2->ind[idebugtest] << " jnd = " << bl2->jnd[idebugtest] << " crois = " << bl2->crois[idebugtest] << G4endl;
    }

    G4cout <<"jparticip[i]" << G4endl;
    for(G4int idebugtest = 0; idebugtest < 300; idebugtest++) {
      G4cout <<" " << jparticip[idebugtest];
    }
    G4cout << G4endl;
  }
  

  // deutons
  if (kcol != 0) {
    if(verboseLevel > 3) {
      G4cout <<"After time evaluation: kcol != 0!" << G4endl;
    }
    goto pnu48;
  }
  // Here transparent event (no interaction avatar found)
  //  G4cout <<"(checkpoint 'heretransparent)" << G4endl;
  nopart = -1;
  // Pour eviter renvoi des resultats du run precedent cv 7/7/98
  iarem = bl3->ia2;
  izrem = iz2;
  esrem = 0.0;
  erecrem = 0.0;

  if(verboseLevel > 3) {
    G4cout <<"kcol == 0. No collisions..." << G4endl;
    G4cout <<"End of algorithm because kcol == 0." << G4endl;
    G4cout <<"kcol = " << kcol << G4endl;
    G4cout <<"Return after pnu39." << G4endl;
  }
  // fin ajout cv
  goto pnureturn;

  // Initialization  at the beginning of the run
 pnu48:
  //  G4cout <<"(checkpoint 'pnu48)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"Beginning a run..." << G4endl;
  }
  timi = 0.0;
  tim = 0.0;
  ncol = 0;

  // compteur des collisions a deux corps (call collis acceptes par pauli)
  ncol_2c = 0;

  //C  npion=0;
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
  for(G4int i = 1; i <= ia; i++)  { //do 47 i=1,ia
    bl5->tlg[i] = 0.0;
    nc[i] = 0;
  }
  itt = 1;

  // tableau des energies a l'initialisation
  if(varavat->kveux == 1) {
    for(G4int i = bl3->ia1+1; i <= ia; i++) {
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

  if(verboseLevel > 3) {
    G4cout <<"Now arriving to label pnu449." << G4endl;
  }
 pnu449:
  //  G4cout <<"(checkpoint 'pnu449)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"Now at 449" << G4endl;
    G4cout <<"G4Incl: Now at label pnu449." << G4endl;
  }
  next = 1;
  indic[next] = 1;

 pnu44:
  //  G4cout <<"(checkpoint 'pnu44)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"Now at 44" << G4endl;
    G4cout <<"Starting a new loop at pnu44..." << G4endl;
  }
  if(next == 0) {
    if(verboseLevel > 3) {
      G4cout <<"next == 0. Returning to label pnu449." << G4endl;
      G4cout <<"next == 0" << G4endl;
    }
    goto pnu449;
  }
  
  idep = indic[next] + 1;
  tau = bl2->crois[idep-1];

  if(idep > bl2->k) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: idep > bl2->k. Going to pnu448." << G4endl;
    }
    goto pnu448;
  }
  for(G4int i = idep; i <= bl2->k; i++) { //do 42 i=idep,k
    if (bl2->crois[i] > tau) {
      continue;
    }
    tau = bl2->crois[i]; 
    next = next + 1;
    indic[next] = i;
  }

  if(verboseLevel > 3) {
    G4cout <<"next = " << next << G4endl;
  }
 pnu448:
  //  G4cout <<"(checkpoint 'pnu448)" << G4endl;
  //  bl2->dump();
  imin = indic[next];
  bl9->l1 = bl2->ind[imin]; //NOTE: l1 changed to bl9->l1.
  bl9->l2 = bl2->jnd[imin]; //NOTE: l2 changed to bl9->l2.

#ifdef G4INCL_DEBUG_LOG
  // Print the debug log
  print_log_start_step();
  print_log_entry(1, 1, 1, imin);
  print_log_end_step();
#endif

  // test le 20/3/2003: tue sinon le dernier avatar?
  if (bl2->k == 0) {
    if(verboseLevel > -2) {
      G4cout <<"k == 0. Going to the end of the avatar." << G4endl;
    }
    goto pnu230;
  }
  bl2->k = bl2->k - 1; // bugfix k belongs to struct bl2!
  next = next - 1;

  // correction s.vuillier 25/1/96 decalage temps correct
  if (imin > bl2->k) {
    goto pnu46;
  }

  for(G4int i = imin; i <= bl2->k; i++) { //do 43 i=imin,k
    bl2->crois[i] = bl2->crois[i+1];
    bl2->ind[i] = bl2->ind[i+1];
    bl2->jnd[i] = bl2->jnd[i+1];
  }
 pnu46:
  //  G4cout <<"(checkpoint 'pnu46)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu46." << G4endl;
  }
  tim = timi + tau;
  avatarCounter++;
  //  G4cout <<"FOUND AVATAR " <<  avatarCounter << " : tau = " << tau;
  //  G4cout <<"  l1 = " << bl9->l1 << "   l2 = " << bl9->l2 << G4endl;

  // tableau des energies a t=20,40,60 fm/c
  if(varavat->kveux == 1) {
    if(iflag20 == 0 && tim >= 20.0) {
      iflag20 = 1;
      for(G4int i = 1; i <= ia; i++) {
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
      for(G4int i = 1; i <= ia; i++) {
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
      for(G4int i = 1; i <= ia; i++) {
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
  //  G4cout <<"(checkpoint 'pnu645)" << G4endl;
  // Here, stoping condition for heavy (light) ions:
  if(kindstruct->kindf7 < 0) {
    if(n_activnuc-nbquit <= ia_breakup) {
      //      G4cout <<"stopped cascade for fermi breakup " << G4endl;
      goto pnu255;
    }
  }

  // modif: pas de reflexions avant au moins un avatar du (des) nucleon incident
  // celui-ci ne peut etre qu'une collision nn (ou pin)

  /*
  if((irst_avatar == 0) && (bl9->l2 == -1)) {
    if(verboseLevel > 3) {
      G4cout <<"Interaction type: reflection (l2 = " << bl9->l2 << "). No first interaction with a participant yet." << G4endl;
    }
    goto pnu44;
  }

  irst_avatar = irst_avatar+1;
  */
  //  if((bl3->ia1 + bl3->ia2 - nbquit) <= 10) goto pnu255; // Stop the cascade if only very few nucleons are left.

  if (tim < temfin) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: tim < temfin. Going to pnu49 (line 2886)" << G4endl;
    }
    goto pnu49; // line 2886
  }

  goto pnu255;
 pnu49:
  //  G4cout <<"(checkpoint 'pnu49)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu49. " << G4endl;
  }
  
  if (bl2->k == 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: bl2->k == 0. Going to pnu255." << G4endl;
    }
    goto pnu255;
  }
  // l1 va a la surface du noyau:
  if (bl9->l2 == -1) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: l2 == -1. Going to pnu220." << G4endl;
    }
    goto pnu220;
  }

  if((k4-1) <= 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: (k4 - 1) <= 0. Going to pnu803." << G4endl;
    }
    goto pnu803;
  }
  if((k4-1) > 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: (k4 - 1) > 0. Going to pnu830." << G4endl;
    }
    goto pnu830;
  }

  // l1 est un delta: 
 pnu830:
  //  G4cout <<"(checkpoint 'pnu830)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu830." << G4endl;
  }
  
  if(bl9->l2 == 0) {
    goto pnu220;
  }
  // interaction pi(l1-ia)-nucleon(l2)
  if(bl9->l1 > ia) { // line 2916
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: l1 > ia (line 2916)" << G4endl;
    }
    goto pnu801;
  }
 pnu803:
  //  G4cout <<"(checkpoint 'pnu803)" << G4endl;
  // pas de collision entre 2 non participants:
  if(jparticip[bl9->l1] == 0 && jparticip[bl9->l2] == 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: both particles are spectators. No collision. " << G4endl;
    }
    goto pnu44;
  }

  // parameters for the next colliding pair
  t[9] = bl1->eps[bl9->l1] + bl1->eps[bl9->l2]; //t(10)->t[9] 
  t0 = 1.0/t[9]; // t(10)->t[9]
  b1 = (bl1->p1[bl9->l1] + bl1->p1[bl9->l2])*t0;
  b2 = (bl1->p2[bl9->l1] + bl1->p2[bl9->l2])*t0;
  b3 = (bl1->p3[bl9->l1] + bl1->p3[bl9->l2])*t0;
  s = (1.0 - b1*b1 - b2*b2 - b3*b3)*t[9]*t[9]; //t(10)->t[9]
  sq = std::sqrt(s);

  if(sq < 1925.5) {
    if(verboseLevel > 3) {
      G4cout <<"sq < 1925.5" << G4endl;
      G4cout <<"Particles: l1 = " << bl9->l1 << " l2 = " << bl9->l2 << G4endl;
      G4cout <<"eps[bl9->l1] = " << bl1->eps[bl9->l1] << " eps[bl9->l2] = " << bl1->eps[bl9->l2] << G4endl;
      G4cout <<"p1[bl9->l1] = " << bl1->p1[bl9->l1] << " p1[bl9->l2] = " << bl1->p1[bl9->l2] << G4endl;
      G4cout <<"p2[bl9->l1] = " << bl1->p2[bl9->l1] << " p2[bl9->l2] = " << bl1->p2[bl9->l2] << G4endl;
      G4cout <<"p3[bl9->l1] = " << bl1->p3[bl9->l1] << " p3[bl9->l2] = " << bl1->p3[bl9->l2] << G4endl;
      G4cout <<"sq = " << sq << G4endl;
    }
    goto pnu44;
  }

  bl1->ta = tau/bl1->eps[bl9->l1];
  x1l1 = bl3->x1[bl9->l1] + bl1->p1[bl9->l1]*(bl1->ta);
  x2l1 = bl3->x2[bl9->l1] + bl1->p2[bl9->l1]*(bl1->ta);
  x3l1 = bl3->x3[bl9->l1] + bl1->p3[bl9->l1]*(bl1->ta);

  bl1->ta = tau/bl1->eps[bl9->l2];
  x1l2 = bl3->x1[bl9->l2] + bl1->p1[bl9->l2]*(bl1->ta);
  x2l2 = bl3->x2[bl9->l2] + bl1->p2[bl9->l2]*(bl1->ta);
  x3l2 = bl3->x3[bl9->l2] + bl1->p3[bl9->l2]*(bl1->ta);

  // test on the minimum distance of approach
  t[10] = x1l1 - x1l2; //t(11)->t[10]
  t[11] = x2l1 - x2l2; //t(12)->t[11]
  t[12] = x3l1 - x3l2; //t(13)->t[12]
  t[13] = t[10]*t[10] + t[11]*t[11] + t[12]*t[12]; //t(N)->t[N-1]
  t[14] = b1*t[10] + b2*t[11] + b3*t[12]; //t(N)->t[N-1]
  t[15] = b1*b1 + b2*b2 + b3*b3; //t(16)->t[15]
  bb2 = t[13] + t[14]*t[14]/(1.0 - t[15]); //t(N)->t[N-1]

  if(verboseLevel > 3) {
    G4cout <<"Minimum dist. of approach tested..." << G4endl;
  }
  
  if (k3 == 1) goto pnu260;
  if (k4 == 0) goto pnu260;
  mg = bl1->ind1[bl9->l1] + bl1->ind1[bl9->l2];
  isos = bl1->ind2[bl9->l1] + bl1->ind2[bl9->l2];
  if (mg != 1) goto pnu260;
  ldel = bl9->l2;
  if(mg - bl1->ind1[bl9->l1] == 0) ldel = bl9->l1;
  bl6->xx10 = std::sqrt(std::pow(bl1->eps[ldel],2) - std::pow(bl1->p1[ldel], 2) - std::pow(bl1->p2[ldel], 2) - std::pow(bl1->p3[ldel], 2));
  bl6->isa = bl1->ind2[ldel];
  bmax2 = totalCrossSection(sq,mg,isos)/31.415926;
  if (k5 == 0 && mg != 0) bmax2 = bmax2 - lowEnergy(sq,mg,isos)/31.415926;
  goto pnu261;
 pnu260:
  //  G4cout <<"(checkpoint 'pnu260)" << G4endl;
  bmax2 = totalCrossSection(sq,mg,isos)/31.41592;
 pnu261:
  //  G4cout <<"(checkpoint 'pnu261)" << G4endl;
  if (bb2 < bmax2) {
    goto pnu220;
  }
  if (bl2->k == 0) {
    goto pnu230;
  }

  if(verboseLevel > 3) {
    G4cout <<"bb2 >= bmax2 or bl2->k == 0" << G4endl;
  }
  goto pnu44;
  // loop while((bb2 >= bmax2) && (k != 0)) (PK)
  // evaluation of the positions at time = tim
 pnu220:
  //  G4cout <<"(checkpoint 'pnu220)" << G4endl;
  timi = tim;
   if(verboseLevel > 3) {
     G4cout <<"Evaluating positions at time = tim" << G4endl;
     G4cout <<"tim = " << tim << G4endl;
   }

//   dumpSaxw(dumpOut);
//   dumpBl1(dumpOut);
//   dumpBl2(dumpOut);
//   dumpBl3(dumpOut);

  if(varavat->kveux == 1) { //then
    iavat = iavat + 1;
    varavat->timeavat[iavat] = tim;
    varavat->l1avat[iavat] = bl9->l1;
    varavat->l2avat[iavat] = bl9->l2;
    varavat->energyavat[iavat] = sq;

    if(bl9->l1 <= ia) {
      varavat->jpartl1[iavat] = jparticip[bl9->l1];
    }
    else {
      varavat->jpartl1[iavat] = 0;
    }

    if(bl9->l2 > 0) {
      varavat->jpartl2[iavat] = jparticip[bl9->l2];
    }
    else {
      varavat->jpartl2[iavat] = 0;
    }
  }

  // gel des nucleons non participants sur le premier avatar (nn)=(l1,1)      
  /*
  if (irst_avatar == 1) {
    for(G4int i = 1; i <= bl9->l1; i = i + bl9->l1 - 1) { // bugfix!
      bl1->ta = tau/bl1->eps[i];                                                
      bl3->x1[i] = bl3->x1[i] + bl1->p1[i]*(bl1->ta);                                     
      bl3->x2[i] = bl3->x2[i] + bl1->p2[i]*(bl1->ta);                                    
      bl3->x3[i] = bl3->x3[i] + bl1->p3[i]*(bl1->ta);
      if(verboseLevel > 3) {
	G4cout <<"G4Incl: i = " << G4endl;
      }
    }
    for(G4int i = 1; i <= bl2->k; i++) {
      bl2->crois[i] = bl2->crois[i] + tau;
    }
  }
  else {
  */
    for(G4int i = 1; i <= ia; i++) {
      bl1->ta = tau/bl1->eps[i];
      bl3->x1[i] = bl3->x1[i] + bl1->p1[i]*(bl1->ta);
      bl3->x2[i] = bl3->x2[i] + bl1->p2[i]*(bl1->ta);
      bl3->x3[i] = bl3->x3[i] + bl1->p3[i]*(bl1->ta);
    }
    //  }

  //  if(npion != 0) {
  if(npion > 0) {
    for(G4int i = 1; i <= npion; i++) {
      bl1->ta = tau/q4[i];
      y1[i] = y1[i] + q1[i]*(bl1->ta);
      y2[i] = y2[i] + q2[i]*(bl1->ta);
      y3[i] = y3[i] + q3[i]*(bl1->ta);
    }
  }
  //  G4cout <<"(checkpoint 'pnu840)" << G4endl;

  if(bl9->l2 == 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: l2 == 0. Going to pnu805." << G4endl;
    }
    goto pnu805;
  }
  // Candidate: if(l2!=0)...

  // reflexions sur le potentiel, sortie eventuelle de la particule:
  if (bl9->l2 == -1) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: l2 == -1. Going to pnu600." << G4endl;
    }
    goto pnu600;
  }

  if(bl9->l1 > ia) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: l1 > ia. Going to pnu831." << G4endl;
    }
    goto pnu831;
  }

  // collision of particles l1 and l2
  if(verboseLevel > 3) {
    G4cout <<"Particles l1 and l2 collide!" << G4endl;
  }
  
  ich1 = bl1->ind1[bl9->l1];
  ich2 = bl1->ind1[bl9->l2];
  ich3 = bl1->ind2[bl9->l1];
  ich4 = bl1->ind2[bl9->l2];
  aml1 = std::sqrt(std::pow(bl1->eps[bl9->l1],2) - std::pow(bl1->p1[bl9->l1],2) - std::pow(bl1->p2[bl9->l1],2) - std::pow(bl1->p3[bl9->l1],2));
  aml2 = std::sqrt(std::pow(bl1->eps[bl9->l2],2) - std::pow(bl1->p1[bl9->l2],2) - std::pow(bl1->p2[bl9->l2],2) - std::pow(bl1->p3[bl9->l2],2));
  gl1 = bl1->eps[bl9->l1]/aml1;
  gl2 = bl1->eps[bl9->l2]/aml2;
  // l-conservation
  if (k6 == 1) {
    t[30] = (aml1*(bl3->x1[bl9->l1]) + aml2*(bl3->x1[bl9->l2]))/(aml1 + aml2); //t(31)->t[30]
    t[31] = (aml1*(bl3->x2[bl9->l1]) + aml2*(bl3->x2[bl9->l2]))/(aml1 + aml2); //t(32)->t[31]
    t[32] = (aml1*(bl3->x3[bl9->l1]) + aml2*(bl3->x3[bl9->l2]))/(aml1 + aml2); //t(33)->t[32]
    tt31 = bl3->x1[bl9->l1] - bl3->x1[bl9->l2];
    tt32 = bl3->x2[bl9->l1] - bl3->x2[bl9->l2];
    tt33 = bl3->x3[bl9->l1] - bl3->x3[bl9->l2];
    t[33] = (aml2*(bl1->p1[bl9->l1]) - aml1*(bl1->p1[bl9->l2]))/(aml1 + aml2); //t(34)->t[33]
    t[34] = (aml2*(bl1->p2[bl9->l1]) - aml1*(bl1->p2[bl9->l2]))/(aml1 + aml2); //t(35)->t[34]
    t[35] = (aml2*(bl1->p3[bl9->l1]) - aml1*(bl1->p3[bl9->l2]))/(aml1 + aml2); //t(36)->t[35]
    tt34 = bl1->p1[bl9->l1] + bl1->p1[bl9->l2];
    tt35 = bl1->p2[bl9->l1] + bl1->p2[bl9->l2];
    tt36 = bl1->p3[bl9->l1] + bl1->p3[bl9->l2];
  }

  // l-conservation
  t[20] = bl1->p1[bl9->l1]; //t(21)->t[20]
  t[21] = bl1->p2[bl9->l1]; //t(22)->t[21]
  t[22] = bl1->p3[bl9->l1]; //t(23)->t[22]
  t[23] = bl1->eps[bl9->l1]; //t(24)->t[23]
  t[24] = bl1->p1[bl9->l2]; //t(25)->t[24]
  t[25] = bl1->p2[bl9->l2]; //t(26)->t[25]
  t[26] = bl1->p3[bl9->l2]; //t(27)->t[26]
  t[27] = bl1->eps[bl9->l2]; //t(28)->t[27]

  // info delta ou nucleon:
  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[bl9->l1];
    varavat->del2avat[iavat] = bl1->ind1[bl9->l2];                                               
  }

  minus_b1 = -1.0*b1;
  minus_b2 = -1.0*b2;
  minus_b3 = -1.0*b3;
  loren(&(bl1->p1[bl9->l1]), &(bl1->p2[bl9->l1]), &(bl1->p3[bl9->l1]), &minus_b1, &minus_b2, &minus_b3, &(bl1->eps[bl9->l1]));
  loren(&(bl1->p1[bl9->l2]), &(bl1->p2[bl9->l2]), &(bl1->p3[bl9->l2]), &minus_b1, &minus_b2, &minus_b3, &(bl1->eps[bl9->l2]));

  if(verboseLevel > 3) {
    G4cout <<"Calling collis..." << G4endl;
    G4cout <<"Energy eps[bl9->l1] = " << bl1->eps[bl9->l1] << G4endl;
    G4cout <<"Momentum: p1 = " << bl1->p1[bl9->l1] << " p2 = " << bl1->p2[bl9->l1] << " p3 = " << bl1->p3[bl9->l1] << G4endl;
    G4cout <<"Energy eps[bl9->l2] = " << bl1->eps[bl9->l2] << G4endl;
    G4cout <<"Momentum: p1 = " << bl1->p1[bl9->l2] << " p2 = " << bl1->p2[bl9->l2] << " p3 = " << bl1->p3[bl9->l2] << G4endl;
  }
//   bl9->l1 = l1;
//   bl9->l2 = l2;
  collis(&(bl1->p1[bl9->l1]), &(bl1->p2[bl9->l1]), &(bl1->p3[bl9->l1]),
	 &(bl1->eps[bl9->l1]), &(bl1->p1[bl9->l2]), &(bl1->p2[bl9->l2]), &(bl1->p3[bl9->l2]), &(bl1->eps[bl9->l2]),
	 &(t[11]), &(t[12]), &(t[13]), &(t[14]), &np, &ip, &k2, &k3, &k4, &k5,
	 &(bl1->ind1[bl9->l1]), &(bl1->ind1[bl9->l2]), &(bl1->ind2[bl9->l1]), &(bl1->ind2[bl9->l2]));
//   l1 = bl9->l1;
//   l2 = bl9->l2;
  if(verboseLevel > 3) {
    G4cout <<"End of collis call" << G4endl;
    G4cout <<"Energy eps[bl9->l1] = " << bl1->eps[bl9->l1] << G4endl;
    G4cout <<"Energy eps[bl9->l2] = " << bl1->eps[bl9->l2] << G4endl;
  }
  
  if(verboseLevel > 3) {
    G4cout <<"Variables after collis call: "<< G4endl;
    G4cout <<"bl1->p1[" << bl9->l1 <<"] = " << bl1->p1[bl9->l1] <<" bl1->p2[" << bl9->l1 <<"] = " << bl1->p2[bl9->l1] <<" bl1->p3[" << bl9->l1 <<"] = " << bl1->p3[bl9->l1] <<" bl1->eps[" << bl9->l1 <<"] = " << bl1->eps[bl9->l1] << G4endl;
    G4cout <<"bl1->p1[" << bl9->l2 <<"] = " << bl1->p1[bl9->l2] <<" bl2->p2[" << bl9->l2 <<"] = " << bl1->p2[bl9->l2] <<" bl2->p3[" << bl9->l2 <<"] = " << bl1->p3[bl9->l2] <<" bl1->eps[" << bl9->l2 <<"] = " << bl1->eps[bl9->l2] << G4endl;
  }
  
  loren(&(bl1->p1[bl9->l1]), &(bl1->p2[bl9->l1]), &(bl1->p3[bl9->l1]), &b1, &b2, &b3, &(bl1->eps[bl9->l1]));
  loren(&(bl1->p1[bl9->l2]), &(bl1->p2[bl9->l2]), &(bl1->p3[bl9->l2]), &b1, &b2, &b3, &(bl1->eps[bl9->l2]));

  if (bl1->ind1[bl9->l1] == 1) { // bugfix 1 -> 0
    goto pnu243;
  }
  xbl1 = pauliBlocking(bl9->l1, rbl, pbl);
  standardRandom(&rndm,&(hazard->igraine[10])); 
  if (rndm > (1.0 - xbl1)) {
    goto pnu248;
  }
 pnu243:
  //  G4cout <<"(checkpoint 'pnu243)" << G4endl;
  if (bl1->ind1[bl9->l2] == 1) { // bugfix 1 -> 0
    goto pnu241;
  }
  xbl2 = pauliBlocking(bl9->l2, rbl, pbl);
  standardRandom(&rndm,&(hazard->igraine[10]));
  if (rndm > (1.0 - xbl2)) {
    goto pnu248;
  }
  goto pnu241;
 pnu248:
  //  G4cout <<"(checkpoint 'pnu248)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"Pauli blocked transition!" << G4endl;
  }
  
  mpaul1 = mpaul1 + 1;
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 1;
  }
  // restitution de l1 et l2 si rejet de la col. par pauli:
  bl1->p1[bl9->l1] = t[20]; //t(21)->t[20]
  bl1->p2[bl9->l1] = t[21]; //t(22)->t[21]
  bl1->p3[bl9->l1] = t[22]; //t(23)->t[22]
  bl1->eps[bl9->l1] = t[23]; //t(24)->t[23]
  bl1->p1[bl9->l2] = t[24]; //t(25)->t[24]
  bl1->p2[bl9->l2] = t[25]; //t(26)->t[25]
  bl1->p3[bl9->l2] = t[26]; //t(27)->t[26]
  bl1->eps[bl9->l2] = t[27]; //t(28)->t[27]
  bl1->ind1[bl9->l1] = ich1;
  bl1->ind1[bl9->l2] = ich2;
  bl1->ind2[bl9->l1] = ich3;
  bl1->ind2[bl9->l2] = ich4;

  if (bl2->k == 0) {
    goto pnu230;
  }
  for(G4int i = 1; i <= bl2->k; i++) {
    bl2->crois[i] = bl2->crois[i] - tau;
  }

  // pour le temps de calcul (a.b. 02/2002)
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: bl2->k != 0" << G4endl;
  }
  goto pnu44;

 pnu241:
  //  G4cout <<"(checkpoint 'pnu241)" << G4endl;
  // la premiere collision a deux corps ne peut pas baisser l'energie
  // du nucleon de recul (bloque par pauli dans un noyau cible froid).
  // (ici, toujours l2 < l1)
  ncol_2c = ncol_2c + 1;
  if(ncol_2c == 1) {
    for (G4int icomp = 1; icomp <= bl3->ia1; icomp++) {
      // test on the first collision modified 4/07/2001 for direct and exchange.
      if(icomp == bl9->l1 || icomp == bl9->l2) {
	xavant = min(t[23],t[27]); //t(24)->t[23], t(28)->t[27]
	xapres = min(bl1->eps[bl9->l1],bl1->eps[bl9->l2]);
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
      varavat->r1_first_avat[0] = bl3->x1[1]; //(N)->[N-1]
      varavat->r1_first_avat[1] = bl3->x2[1]; //(N)->[N-1]
      varavat->r1_first_avat[2] = bl3->x3[1]; //(N)->[N-1]
    } //endif
  }
  else {
    // les collisions suivantes ne penvent conduire a un noyau de a nucleons
    // sous l'energie de fermi et dans une config. d'energie inferieure a
    // efer-(ia2-nbalttf)*tf).
    egs = 0.0;
    nbalttf = 0;
    //    for(G4int i = 1; i <= ia-1; i++) {
    for(G4int i = 1; i <= ia; i++) {
      if(bl5->nesc[i] == 0) {
	if(std::sqrt(std::pow(bl1->p1[i],2)+std::pow(bl1->p2[i],2)+std::pow(bl1->p3[i],2)) < bl10->pf) {
	  nbalttf = nbalttf + 1;
	  egs = egs + bl1->eps[i] - fmp;
	}
      }
    }

    if(egs < (efer- double(bl3->ia2-nbalttf)*tf)) {
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

  jparticip[bl9->l1] = 1;
  jparticip[bl9->l2] = 1;
  if(verboseLevel > 3) {
    G4cout <<"Particle " << bl9->l1 << " is now participant." << G4endl;
    G4cout <<"Particle " << bl9->l2 << " is now participant." << G4endl;
  }
  
  if (ws->nosurf <= 0) {
    // surface
    pppp = std::sqrt(std::pow(bl1->p1[bl9->l1],2) + std::pow(bl1->p2[bl9->l1],2) + std::pow(bl1->p3[bl9->l1],2));
    rrrr = std::sqrt(std::pow(bl3->x1[bl9->l1],2) + std::pow(bl3->x2[bl9->l1],2) + std::pow(bl3->x3[bl9->l1],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[bl9->l1] = bl3->x1[bl9->l1]*rcorr/rrrr;
	bl3->x2[bl9->l1] = bl3->x2[bl9->l1]*rcorr/rrrr;
	bl3->x3[bl9->l1] = bl3->x3[bl9->l1]*rcorr/rrrr;
      }
    }
    pppp = std::sqrt(std::pow(bl1->p1[bl9->l2],2) + std::pow(bl1->p2[bl9->l2],2) + std::pow(bl1->p3[bl9->l2],2));
    rrrr = std::sqrt(std::pow(bl3->x1[bl9->l2],2) + std::pow(bl3->x2[bl9->l2],2) + std::pow(bl3->x3[bl9->l2],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[bl9->l2] = bl3->x1[bl9->l2]*rcorr/rrrr;
	bl3->x2[bl9->l2] = bl3->x2[bl9->l2]*rcorr/rrrr;
	bl3->x3[bl9->l2] = bl3->x3[bl9->l2]*rcorr/rrrr;
      }
    }
  }
  
  if (np == 0) {
    goto pnu240;
  }

  npion = npion + 1;
  loren(&(t[11]), &(t[12]), &(t[13]), &b1, &b2, &b3, &(t[14])); //t(N)->t[N-1]
  q1[npion] = t[11]; //t(12)->t[11]
  q2[npion] = t[12]; //t(13)->t[12]
  q3[npion] = t[13]; //t(14)->t[13]
  q4[npion] = t[14]; //t(15)->t[14]
 pnu240:
  //  G4cout <<"(checkpoint 'pnu240)" << G4endl;
  ncol = ncol + 1;
  if (bl9->l2 != 1) {
    goto pnu870;
  }

  // critere pour la leading particle: avant impulsion longitudinale max l=1
  // change en fevrier 2002: leading part. = energie totale max (l=1)
  if (bl1->p3[bl9->l2] > bl1->p3[bl9->l1]) { 
    goto pnu870;
  }

  // attention, il faut mieux faire et selectionner la plus grande energie
  // des particules participantes (jparticip()=1) et dans le noyau (nesc()=0)!

  xr1 = bl1->p1[bl9->l1];
  xr2 = bl1->p2[bl9->l1];
  xr3 = bl1->p3[bl9->l1];
  xr4 = bl1->eps[bl9->l1];
  xr5 = bl3->x1[bl9->l1];
  xr6 = bl3->x2[bl9->l1];
  xr7 = bl3->x3[bl9->l1];
  xr8 = gl1;
  ixr1 = bl1->ind1[bl9->l1];
  ixr2 = bl1->ind2[bl9->l1];
  ixr3 = ich1;
  bl1->p1[bl9->l1] = bl1->p1[bl9->l2];
  bl1->p2[bl9->l1] = bl1->p2[bl9->l2];
  bl1->p3[bl9->l1] = bl1->p3[bl9->l2];
  bl1->eps[bl9->l1] = bl1->eps[bl9->l2];

  bl3->x1[bl9->l1] = bl3->x1[bl9->l2];
  bl3->x2[bl9->l1] = bl3->x2[bl9->l2];
  bl3->x3[bl9->l1] = bl3->x3[bl9->l2];
  gl1 = gl2;
  bl1->ind1[bl9->l1] = bl1->ind1[bl9->l2];
  bl1->ind2[bl9->l1] = bl1->ind2[bl9->l2];
  ich1 = ich2;
  bl1->p1[bl9->l2] = xr1;
  bl1->p2[bl9->l2] = xr2;
  bl1->p3[bl9->l2] = xr3;
  bl1->eps[bl9->l2] = xr4;
  
  bl3->x1[bl9->l2] = xr5;
  bl3->x2[bl9->l2] = xr6;
  bl3->x3[bl9->l2] = xr7;
  gl2 = xr8;
  bl1->ind1[bl9->l2] = ixr1;
  bl1->ind2[bl9->l2] = ixr2;
  ich2 = ixr3;

  if(ich1 + ich2 - bl1->ind1[bl9->l1] - bl1->ind1[bl9->l2] != 0 || (ich1 + ich2) != 1) {
    goto pnu870;
  }
  if (bl2->k == 0) {
    goto pnu870;
  }

  for(G4int i = 1; i <= bl2->k; i++) {
    if((bl2->ind[i] != 1) || (bl2->jnd[i] != 0)) {
      if((bl2->ind[i] != bl9->l1) || (bl2->jnd[i] != 0)) {
	continue;
      }
      bl2->ind[i] = 1;
      break;
    }

    bl2->ind[i] = bl9->l1;
    break;
  }

  pnu870:
  //  G4cout <<"(checkpoint 'pnu870)" << G4endl;
  bl5->tlg[bl9->l1] = th*(bl1->eps[bl9->l1])/std::sqrt(std::pow(bl1->eps[bl9->l1],2)-std::pow(bl1->p1[bl9->l1],2)-std::pow(bl1->p2[bl9->l1],2)-std::pow(bl1->p3[bl9->l1],2));
  bl5->tlg[bl9->l2] = th*(bl1->eps[bl9->l2])/std::sqrt(std::pow(bl1->eps[bl9->l2],2)-std::pow(bl1->p1[bl9->l2],2)-std::pow(bl1->p2[bl9->l2],2)-std::pow(bl1->p3[bl9->l2],2));
  nc[bl9->l1] = nc[bl9->l1] + 1;
  nc[bl9->l2] = nc[bl9->l2] + 1;
  led = 0;

  if((ich1+ich2-bl1->ind1[bl9->l1]-bl1->ind1[bl9->l2]) < 0) {
    mrnd = mrnd + 1;
  }
  if((ich1+ich2-bl1->ind1[bl9->l1]-bl1->ind1[bl9->l2]) == 0) {
    if((ich1+ich2-1) < 0) {
      mrnn = mrnn + 1;
    }
    if((ich1+ich2-1) == 0) {
      mrdd = mrdd + 1;
      led = 1;
    }
    if((ich1+ich2-1) > 0) {
      mcdd = mcdd + 1;
      led = 1;
    }
  }
  if((ich1+ich2-bl1->ind1[bl9->l1]-bl1->ind1[bl9->l2]) > 0) {
    mrdn = mrdn + 1;
  }


  // reevaluation of the times t(a,b) for (a or b)=(l1 or l2)

  // l-conservation 
  if (k6 == 1) {
    aml1 = am(bl1->p1[bl9->l1],bl1->p2[bl9->l1],bl1->p3[bl9->l1],bl1->eps[bl9->l1]);
    aml2 = am(bl1->p1[bl9->l2],bl1->p2[bl9->l2],bl1->p3[bl9->l2],bl1->eps[bl9->l2]);
    
    t[36] = (aml2*(bl1->p1[bl9->l1]) - aml1*(bl1->p1[bl9->l2]))/(aml1+aml2); //t(37)->t[36]
    t[37] = (aml2*(bl1->p2[bl9->l1]) - aml1*(bl1->p2[bl9->l2]))/(aml1+aml2); //t(38)->t[37]
    t[38] = (aml2*(bl1->p3[bl9->l1]) - aml1*(bl1->p3[bl9->l2]))/(aml1+aml2); //t(39)->t[38]
    t[39] = std::sqrt(t[33]*t[33] + t[34]*t[34] + t[35]*t[35]); //t(N)->t[N-1]
    t[40] = std::sqrt(t[36]*t[36] + t[37]*t[37] + t[38]*t[38]); //t(N)->t[N-1]
    rhopi = tt31*t[33] + tt32*t[34] + tt33*t[35]; //t(N)->t[N-1]
    t[42] = tt31 - rhopi*t[33]/std::pow(t[39],2); //t(N)->t[N-1]
    t[43] = tt32 - rhopi*t[34]/std::pow(t[39],2); //t(N)->t[N-1]
    t[44] = tt33 - rhopi*t[35]/std::pow(t[39],2); //t(N)->t[N-1]
    t[45] = std::sqrt(t[42]*t[42] + t[43]*t[43] + t[44]*t[44]); //t(N)->t[N-1]
    t[42] = t[42]/t[45]; //t(N)->t[N-1]
    t[43] = t[43]/t[45]; //t(N)->t[N-1]
    t[44] = t[44]/t[45];
    cif = (t[33]*t[36] + t[34]*t[37] + t[35]*t[38])/t[39]/t[40]; //t(N)->t[N-1]
   
    // trouble with forward scattering 22/3/95
    if(std::fabs(cif) > 1.) {
      cif = sign(1.0,cif);
    }
    sif = std::sqrt(1.0 - cif*cif);
    t[36] = (t[33]*cif/t[39] + t[42]*sif)*t[40]; //t(N)->t[N-1]
    t[37] = (t[34]*cif/t[39] + t[43]*sif)*t[40]; //t(N)->t[N-1]
    t[38] = (t[35]*cif/t[39] + t[44]*sif)*t[40]; //t(N)->t[N-1]
    tri = std::sqrt(tt31*tt31 + tt32*tt32 + tt33*tt33);
    cchi = rhopi/tri/t[39]; //t(40)->t[39]
    schi = std::sqrt(1.0 - cchi*cchi);
    c1 = cif*cchi - sif*schi;
    c2 = sif*cchi + cif*schi;
    tt31 = (c1*t[33]/t[39] + c2*t[42])*tri*t[39]/t[40]; //t(N)->t[N-1]
    tt32 = (c1*t[34]/t[39] + c2*t[43])*tri*t[39]/t[40]; //t(N)->t[N-1]
    tt33 = (c1*t[35]/t[39] + c2*t[44])*tri*t[39]/t[40]; //t(N)->t[N-1]
    bl3->x1[bl9->l1] = t[30] + aml2*tt31/(aml1 + aml2); //t(31)->t[30]
    bl3->x2[bl9->l1] = t[31] + aml2*tt32/(aml1 + aml2); //t(32)->t[30]
    bl3->x3[bl9->l1] = t[32] + aml2*tt33/(aml1 + aml2); //t(33)->t[32]
    bl3->x1[bl9->l2] = t[30] - aml1*tt31/(aml1 + aml2); //t(31)->t[30]
    bl3->x2[bl9->l2] = t[31] - aml1*tt32/(aml1 + aml2); //t(32)->t[31]
    bl3->x3[bl9->l2] = t[32] - aml1*tt33/(aml1 + aml2); //t(33)->t[32]
    bl1->p1[bl9->l1] = aml1*tt34/(aml1 + aml2) + t[36]; //t(37)->t[36]
    bl1->p2[bl9->l1] = aml1*tt35/(aml1 + aml2) + t[37]; //t(38)->t[37]
    bl1->p3[bl9->l1] = aml1*tt36/(aml1 + aml2) + t[38]; //t(39)->t[38]
    bl1->eps[bl9->l1] = w(bl1->p1[bl9->l1],bl1->p2[bl9->l1],bl1->p3[bl9->l1],aml1);
    bl1->p1[bl9->l2] = aml2*tt34/(aml1 + aml2) - t[36]; //t(37)->t[36]
    bl1->p2[bl9->l2] = aml2*tt35/(aml1 + aml2) - t[37]; //t(38)->t[37]
    bl1->p3[bl9->l2] = aml2*tt36/(aml1 + aml2) - t[38]; //t(39)->t[38]
    bl1->eps[bl9->l2] = w(bl1->p1[bl9->l2],bl1->p2[bl9->l2],bl1->p3[bl9->l2],aml2);
  }
  // l-conservation

  if(bl2->k != 0) {
    kd = 0;
    ccr = tau;
//     for(G4int i = 1; i <= bl2->k; i++) {
//       i20 = i - kd;
//       if (k4 != 2 || led != 1) {
// 	if((bl2->ind[i] == l1) || (bl2->ind[i] == l2) || (bl2->jnd[i] == l1) || (bl2->jnd[i] == l2)) {
// 	  kd = kd + 1;
// 	  continue;
// 	}
//       }
//       else {
// 	if(bl2->jnd[i] == 0) {
// 	  if (bl2->ind[i] == l1 && bl1->ind1[l1] == 1) {
// 	    bl2->crois[i]=(bl2->crois[i]-ccr)*(bl1->eps[l1])/std::sqrt(std::pow(bl1->eps[l1],2)-std::pow(bl1->p1[l1],2)-std::pow(bl1->p2[l1],2)-std::pow(bl1->p3[l1],2))/gl1+ccr;
// 	  }
// 	  if (bl2->ind[i] == l2 && bl1->ind1[l2] == 1) {
// 	    bl2->crois[i]=(bl2->crois[i]-ccr)*(bl1->eps[l2])/std::sqrt(std::pow(bl1->eps[l2],2)-std::pow(bl1->p1[l2],2)-std::pow(bl1->p2[l2],2)-std::pow(bl1->p3[l2],2))/gl2+ccr;
// 	  }
// 	}
//       }
    
//       bl2->crois[i20]=bl2->crois[i]-ccr;
//       bl2->ind[i20]=bl2->ind[i];
//       bl2->jnd[i20]=bl2->jnd[i];
//     }

    for(G4int i = 1; i <= bl2->k; i++) { //do 50 i=1,k                                                       p-n09070
      i20 = i - kd;
      if(k4 != 2 || led != 1) { //if (k4.ne.2.or.led.ne.1) go to 512                                p-n09090
	goto pnu512;
      }
      if(bl2->jnd[i] == 0) {//if (jnd(i).eq.0) go to 511                                        p-n09100
	goto pnu511;
      }
    pnu512:
      if(bl2->ind[i] == bl9->l1) {//if (ind(i).eq.l1) go to 52                                        p-n09120
	goto pnu52; 

      }
      if(bl2->ind[i] == bl9->l2) {//if (ind(i).eq.l2) go to 52                                        p-n09130
	goto pnu52;
	}
      if(bl2->jnd[i] == bl9->l2) {//if (jnd(i).eq.l2) go to 52                                        p-n09140
	goto pnu52;
      }
      if(bl2->jnd[i] == bl9->l1) {//if (jnd(i).eq.l1) go to 52                                        p-n09150
	goto pnu52;
      }
      goto pnu513;
    pnu511:
      //if (ind(i).eq.l1.and.ind1(l1).eq.1) crois(i)=(crois(i)-ccr)*eps(l1p-n09170
      //-)/std::sqrt(eps(l1)**2-p1(l1)**2-p2(l1)**2-p3(l1)**2)/gl1+ccr          p-n09180
      if(bl2->ind[i] == bl9->l1 && bl1->ind1[bl9->l1] == 1) {
	bl2->crois[i]=(bl2->crois[i]-ccr)*bl1->eps[bl9->l1]/std::sqrt(std::pow(bl1->eps[bl9->l1],2)-std::pow(bl1->p1[bl9->l1],2)-std::pow(bl1->p2[bl9->l1],2)-std::pow(bl1->p3[bl9->l1],2))/gl1+ccr;
      }
      //if (ind(i).eq.l2.and.ind1(l2).eq.1) crois(i)=(crois(i)-ccr)*eps(l2p-n09190
      //-)/std::sqrt(eps(l2)**2-p1(l2)**2-p2(l2)**2-p3(l2)**2)/gl2+ccr          p-n09200
      if(bl2->ind[i] == bl9->l2 && bl1->ind1[bl9->l2] == 1) {
	bl2->crois[i]=(bl2->crois[i]-ccr)*bl1->eps[bl9->l2]/std::sqrt(std::pow(bl1->eps[bl9->l2],2)-std::pow(bl1->p1[bl9->l2],2)-std::pow(bl1->p2[bl9->l2],2)-std::pow(bl1->p3[bl9->l2],2))/gl1+ccr;
      }
 pnu513:
      bl2->crois[i20]=bl2->crois[i]-ccr;
      bl2->ind[i20]=bl2->ind[i];
      bl2->jnd[i20]=bl2->jnd[i];
      continue; // goto pnu50
    pnu52:
      kd=kd+1;
    }
  
    bl2->k = bl2->k - kd;
  }

  newt(bl9->l1,bl9->l2);

  tref=ref(bl3->x1[bl9->l1], bl3->x2[bl9->l1], bl3->x3[bl9->l1], bl1->p1[bl9->l1], bl1->p2[bl9->l1], bl1->p3[bl9->l1], bl1->eps[bl9->l1],r22); // line 3502
  
  if(verboseLevel > 3) {
    if(tref < 0.0) {
      G4cout <<"G4Incl: Reflection time < 0! (line 3502)" << G4endl;
    }
  }

  if(tref <= bl4->tmax5) {
    bl2->k = bl2->k + 1;
    bl2->crois[bl2->k] = tref;
    bl2->ind[bl2->k] = bl9->l1;
    bl2->jnd[bl2->k] = -1;
  }

  tref=ref(bl3->x1[bl9->l2], bl3->x2[bl9->l2], bl3->x3[bl9->l2], bl1->p1[bl9->l2], bl1->p2[bl9->l2], bl1->p3[bl9->l2], bl1->eps[bl9->l2],r22); // line 3516
  
  if(verboseLevel > 3) {
    if(tref < 0.0) {
      G4cout <<"G4Incl: Reflection time < 0! (line 3516)" << G4endl;
    }
  }

  if(tref <= bl4->tmax5) {
    bl2->k = bl2->k + 1;
    bl2->crois[bl2->k] = tref;
    bl2->ind[bl2->k] = bl9->l2;
    bl2->jnd[bl2->k] = -1;
  }

  if (k4 == 2) {
    goto pnu848;
  }

  if(bl2->k < 0) {
    goto pnu230;
  }
  if(bl2->k == 0) {
    goto pnu230;
  }
  if(bl2->k > 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: bl2->k > 0 at line 3488. Going back to label pnu449." << G4endl;
    }
    goto pnu449;
  }

 pnu848:
  //  G4cout <<"(checkpoint 'pnu848)" << G4endl;
  if (npion == 0) {
    goto pnu844;
  }
  if (bl1->ind1[bl9->l1] == 1) {
    goto pnu843;
  }
  for(G4int k20 = 1; k20 <= npion; k20++) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: calling G4Incl::new3" << G4endl;
    }
    new3((y1[k20]), (y2[k20]), (y3[k20]), (q1[k20]), (q2[k20]), (q3[k20]), (q4[k20]), k20, bl9->l1);
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: After new3:" << G4endl;
      G4cout <<"y1[" << k20 << "] = " << y1[k20] <<" y2[" << k20 << "] = " << y2[k20] <<" y3[" << k20 << "] = " << y3[k20] << G4endl;
      G4cout <<"q1[" << k20 << "] = " << q1[k20] <<" q2[" << k20 << "] = " << q2[k20] <<" q3[" << k20 << "] = " << q3[k20] << G4endl;
    }
  }

 pnu843:
  //  G4cout <<"(checkpoint 'pnu843)" << G4endl;
  if(bl1->ind1[bl9->l2] != 1) {
    for(G4int k20 = 1; k20 <= npion; k20++) {
      //new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], -k20, l2);
      new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], k20, bl9->l2);
    }
  }
 pnu844:
  //  G4cout <<"(checkpoint 'pnu844)" << G4endl;
  if(bl1->ind1[bl9->l1]+bl1->ind1[bl9->l2] <= ich1+ich2) {
    goto pnu849;
  }
  if(bl1->ind1[bl9->l1]-ich1 != 1) {
    goto pnu820;
  }
  lnew = bl9->l1;
  goto pnu821;
 pnu820:
  //  G4cout <<"(checkpoint 'pnu820)" << G4endl;
  if(bl1->ind1[bl9->l2]-ich2 != 1) {
    goto pnu849;
  }
  lnew = bl9->l2;

 pnu821:
  //  G4cout <<"(checkpoint 'pnu821)" << G4endl;
  standardRandom(&rndm,&(hazard->igraine[16]));
  // largeur variable du delta (phase space factor G4introduced 4/2001)
  amlnew = std::sqrt(std::pow(bl1->eps[lnew],2)-std::pow(bl1->p1[lnew],2)-std::pow(bl1->p2[lnew],2)-std::pow(bl1->p3[lnew],2)); 
  
  geff = bl1->eps[lnew]/amlnew;
  qqq = std::sqrt((std::pow(amlnew,2) - std::pow((fmp+fmpi),2))*(std::pow(amlnew,2) - std::pow((fmp-fmpi),2)))/(2.0*amlnew);
  
  psf = std::pow(qqq,3)/(std::pow(qqq,3) + 5832000.0);
  tdel = -hc/(g0*psf)*std::log(rndm)*geff;                                  
  
  if(tdel <= bl4->tmax5) {
    bl2->k = bl2->k + 1;
    bl2->crois[bl2->k] = tdel;
    bl2->ind[bl2->k] = lnew;
    bl2->jnd[bl2->k] = 0;
  }

 pnu849:
  if (bl2->k == 0) { 
    goto pnu230;
  }
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: bl2->k != 0. Going back to label pnu449." << G4endl;
  }
  goto pnu449;

  // decay of the delta particle                                       p-n09780
 pnu805:
  //  G4cout <<"(checkpoint 'pnu805)" << G4endl;
  npion = npion + 1;
  ichd = bl1->ind2[bl9->l1];
  t[30] = bl1->p1[bl9->l1];   //t(31)->t[30]
  t[31] = bl1->p2[bl9->l1]; //t(32)->t[31]
  t[32] = bl1->p3[bl9->l1]; //t(33)->t[32]
  t[33] = bl1->eps[bl9->l1]; //t(34)->t[33]
  var_ab = std::pow(bl1->eps[bl9->l1],2) - std::pow(bl1->p1[bl9->l1],2) - std::pow(bl1->p2[bl9->l1],2) - std::pow(bl1->p3[bl9->l1],2);
  ym[npion] = 0.0;

  if(var_ab > 0.0) {
    ym[npion] = std::sqrt(var_ab);
  }

  //PK: This workaround avoids a NaN problem encountered with
  //geant4. The problem does not seem to exist if we run in standalone
  //mode.
  //  if(((std::pow(ym[npion],2)-std::pow(fmp+fmpi,2))*(std::pow(ym[npion],2)-std::pow(fmp-fmpi,2))) < 0) {
  //    ym[npion] = ym[npion]+fmpi+100.0;
  //  }
  // PK
  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[bl9->l1];
    varavat->energyavat[iavat] = ym[npion];
  }

  if(verboseLevel > 3) {
    G4cout <<"Calling decay2 from pnu" << G4endl;
    G4cout <<"npion = " << npion << G4endl;
    G4cout <<"q1 = " << q1[npion] << " q2 = " << q2[npion] << " q3 = " << q3[npion] << " q4 = " << q4[npion] << G4endl;
  }
  decay2(&(bl1->p1[bl9->l1]), &(bl1->p2[bl9->l1]), &(bl1->p3[bl9->l1]), &(bl1->eps[bl9->l1]), &(q1[npion]), &(q2[npion]), &(q3[npion]),
 	 &(q4[npion]), &(ym[npion]), &fmp, &fmpi, &(bl9->hel[bl9->l1]));

  if(verboseLevel > 3) {
    G4cout <<"Quantities after decay2: " << G4endl;
    G4cout <<"l1 = " << bl9->l1 << " bl1->p1[l1] = " << bl1->p1[bl9->l1] << " bl1->p2[l1] = " << bl1->p2[bl9->l1] << " bl1->p3[l1] = " << bl1->p3[bl9->l1] << " bl1->eps[l1] = " << bl1->eps[bl9->l1] << G4endl;
  }
  
  // decay
  if (bl1->ind2[bl9->l1]*(bl1->ind2[bl9->l1]) == 9) {
    goto pnu806;
  }
  
  standardRandom(&rndm, &(hazard->ial));
  if (rndm < 0.333333333) {
    goto pnu837;
  }
  
  ipi[npion]=0;
  goto pnu809;

 pnu837:
  //  G4cout <<"(checkpoint 'pnu837)" << G4endl;
  ipi[npion]=bl1->ind2[bl9->l1]*2;
  bl1->ind2[bl9->l1]=-1*(bl1->ind2[bl9->l1]);
  goto pnu809;
 pnu806:
  //  G4cout <<"(checkpoint 'pnu806)" << G4endl;
  bl1->ind2[bl9->l1]=bl1->ind2[bl9->l1]/3;
  ipi[npion]=2*(bl1->ind2[bl9->l1]);
 pnu809: // continue
  //  G4cout <<"(checkpoint 'pnu809)" << G4endl;
  bl1->ind1[bl9->l1]=0;
  bl5->tlg[bl9->l1]=0.;

  // escape ?
  if (bl5->nesc[bl9->l1] > 0) {
    goto pnu850;
  }

  iteste = 0;
  xpb = pauliBlocking(bl9->l1, rbl, pbl);
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
  for(G4int i = 1; i <= ia; i++) {
    if(bl5->nesc[i] == 0) {
      if(std::sqrt(std::pow(bl1->p1[i],2)+std::pow(bl1->p2[i],2)+std::pow(bl1->p3[i],2)) < bl10->pf) {
	nbalttf = nbalttf + 1;
	egs = egs + bl1->eps[i] - fmp;
      }
    }
  }
  if(egs >= (efer - double(bl3->ia2-nbalttf)*tf)) {
    goto pnu850;
  }

  // attention, logique negative!!! liberer le goto si on veut supprimer la
  // sequence precedente (cdpp sur delta-> pi n)

  if(varavat->kveux == 1) {
    varavat->bloc_cdpp[iavat] = 1;
  }

  // it is blocked!      
 pnu1848:
  //  G4cout <<"(checkpoint 'pnu1848)" << G4endl;
  mpaul2 = mpaul2 + 1;
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 1;
  }

  // largeur variable du delta (phase space factor G4introduced 4/2001)
  // (180.**3 = 5832000.)
  qqq = std::sqrt((std::pow(ym[npion],2) - std::pow((fmp+fmpi),2))*(std::pow(ym[npion],2) - std::pow((fmp-fmpi),2)))/(2.*ym[npion]);
  psf = std::pow(qqq,3)/(std::pow(qqq,3)+5832000.0);
  tdel = hc*t[33]/(g0*psf*ym[npion]); //t(34)->t[33]                                 

  if (iteste == 0) {
    tdel = tdel*xpb/(1.000001-xpb);
  }

  if(tdel <= bl4->tmax5) {
    bl2->k = bl2->k + 1;
    bl2->crois[bl2->k] = tdel;
    bl2->ind[bl2->k] = bl9->l1;
    bl2->jnd[bl2->k] = 0;
  }

  bl1->p1[bl9->l1] = t[30]; //t(31)->t[30]
  bl1->p2[bl9->l1] = t[31]; //t(32)->t[31]
  bl1->p3[bl9->l1] = t[32]; //t(33)->t[32]
  bl1->eps[bl9->l1] = t[33]; //t(34)->t[33]
  bl1->ind1[bl9->l1] = 1;
  bl1->ind2[bl9->l1] = ichd;
  npion = npion - 1;

  if (bl2->k == 0) {
    goto pnu230;
  }

  for(G4int i = 1; i <= bl2->k; i++) {
    bl2->crois[i] = bl2->crois[i] - tau;
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Going back to label pnu449 after a loop from 1 to bl2->k" << G4endl;
  }
  goto pnu449;

  // valid decay of the delta
 pnu850:
  //  G4cout <<"(checkpoint 'pnu850)" << G4endl;
  if(varavat->kveux == 1) {
    varavat->bloc_paul[iavat] = 0;
    varavat->bloc_cdpp[iavat] = 0;
  }

  if (ws->nosurf <= 0) {
    // surface
    pppp = std::sqrt(std::pow(bl1->p1[bl9->l1],2) + std::pow(bl1->p2[bl9->l1],2) + std::pow(bl1->p3[bl9->l1],2));
    rrrr = std::sqrt(std::pow(bl3->x1[bl9->l1],2) + std::pow(bl3->x2[bl9->l1],2) + std::pow(bl3->x3[bl9->l1],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[bl9->l1] = bl3->x1[bl9->l1]*rcorr/rrrr;
	bl3->x2[bl9->l1] = bl3->x2[bl9->l1]*rcorr/rrrr;
	bl3->x3[bl9->l1] = bl3->x3[bl9->l1]*rcorr/rrrr;
      }
    }
  }

  ncol = ncol + 1;
  mrdp = mrdp + 1;
  y1[npion] = bl3->x1[bl9->l1];
  y2[npion] = bl3->x2[bl9->l1];
  y3[npion] = bl3->x3[bl9->l1];

  if (bl2->k == 0) {
    goto pnu4047;
  }
  
  kd = 0;
  ccr = tau;
  for(G4int i = 1; i <= bl2->k; i++) {
    i20 = i - kd;

    if((bl2->ind[i] == bl9->l1) || (bl2->jnd[i] == bl9->l1)) {
      kd = kd + 1;
    }
    else {
      bl2->crois[i20] = bl2->crois[i] - ccr;
      bl2->ind[i20] = bl2->ind[i];
      bl2->jnd[i20] = bl2->jnd[i];
    }
  }
  bl2->k = bl2->k - kd;
  
 pnu4047:
  //  G4cout <<"(checkpoint 'pnu4047)" << G4endl;
  if (bl5->nesc[bl9->l1] != 0) {
    goto pnu845;
  }

  new1(bl9->l1);
  bl2->k = bl2->k + 1;
  bl2->crois[bl2->k] = ref(bl3->x1[bl9->l1], bl3->x2[bl9->l1], bl3->x3[bl9->l1], bl1->p1[bl9->l1], bl1->p2[bl9->l1], bl1->p3[bl9->l1], bl1->eps[bl9->l1],r22);
  bl2->ind[bl2->k] = bl9->l1;
  bl2->jnd[bl2->k] = -1;
  if(verboseLevel > 3) {
    if(bl2->crois[bl2->k] < 0.0) {
      G4cout <<"G4Incl: Reflection time < 0! (line 3797)" << G4endl;
    }
  }


  if(npion > 1) {
    n20 = npion - 1;
    for(G4int k20 = 1; k20 <= n20; k20++) {
      new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], k20, bl9->l1);
    }
  }

 pnu845: 
  //  G4cout <<"(checkpoint 'pnu845)" << G4endl;
  new2(y1[npion], y2[npion], y3[npion], q1[npion], q2[npion], q3[npion], q4[npion], npion, bl9->l1);
  if(bl2->k == 0) {
    goto pnu230;
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: bl2->k == 0 after a call to new2. Going back to label pnu449." << G4endl;
  }
  goto pnu449;

  // pion-nucleon collision
 pnu801:
  //  G4cout <<"(checkpoint 'pnu801)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"Pion-nucleon collision!" << G4endl;
  }
  lp = bl9->l1 - ia;
  dis1 = bl3->x1[bl9->l2]-y1[lp] + (bl1->p1[bl9->l2]/bl1->eps[bl9->l2] - q1[lp]/q4[lp])*tau;
  dis2 = bl3->x2[bl9->l2]-y2[lp] + (bl1->p2[bl9->l2]/bl1->eps[bl9->l2] - q2[lp]/q4[lp])*tau;
  dis3 = bl3->x3[bl9->l2]-y3[lp] + (bl1->p3[bl9->l2]/bl1->eps[bl9->l2] - q3[lp]/q4[lp])*tau;
  dist = dis1*dis1 + dis2*dis2 + dis3*dis3;
  t[9] = bl1->eps[bl9->l2] + q4[lp]; //t(10)->t[9]
  t0 = 1.0/t[9]; //t(10)->t[9]
  b1 = (bl1->p1[bl9->l2] + q1[lp])*t0;
  b2 = (bl1->p2[bl9->l2] + q2[lp])*t0;
  b3 = (bl1->p3[bl9->l2] + q3[lp])*t0;
  s = (1.0 - b1*b1 - b2*b2 - b3*b3)*t[9]*t[9]; //t(10)->t[9]
  sq = std::sqrt(s);
  cg = 4+bl1->ind2[bl9->l2]*ipi[lp];

  if(verboseLevel > 3) {
    G4cout <<"Pion-Nucleon collision done! " << G4endl;
  }
  
  if(sq > 3000.0) {
    if(verboseLevel > 3) {
      G4cout <<"sq = " << sq << G4endl;
    }
    goto pnu832;
  }
  if(31.41592*dist > pionNucleonCrossSection(sq)*cg/6.0) {
    goto pnu832;
  }

  if(verboseLevel > 3) {
    G4cout <<"Going to pnu220!" << G4endl;
  }
  
  goto pnu220;

 pnu832:
  //  G4cout <<"(checkpoint 'pnu832)" << G4endl;
  if (bl2->k == 0) {
    goto pnu230;
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: bl2->k != 0. Going back to pnu44." << G4endl;
  }
  goto pnu44;
 pnu831:
  //  G4cout <<"(checkpoint 'pnu831)" << G4endl;
  standardRandom(&rndm, &(hazard->igraine[18]));
  geff = t[9]/sq; //t(10)->t[9]
  gg = g0;
  if (sq > 1500.0) {
    gg=200.0;
  }

  // largeur variable du delta (phase space factor G4introduced 4/2001)
  // (180.**3 = 5832000.)
  qqq = std::sqrt((std::pow(sq,2) - std::pow((fmp+fmpi),2))*(std::pow(sq,2) - std::pow((fmp-fmpi),2)))/(2.0*sq);
  psf = std::pow(qqq,3)/(std::pow(qqq,3)+5832000.);
  tdel = -hc/(gg*psf)*std::log(rndm)*geff;         

  bl1->ind1[bl9->l2] = 1;
  bl1->ind2[bl9->l2] = bl1->ind2[bl9->l2] + ipi[lp];
  nc[bl9->l2] = nc[bl9->l2] + 1;
  bl1->eps[bl9->l2] = t[9]; //t(10)->t[9]
  bl1->p1[bl9->l2] = bl1->p1[bl9->l2] + q1[lp];
  bl1->p2[bl9->l2] = bl1->p2[bl9->l2] + q2[lp];
  bl1->p3[bl9->l2] = bl1->p3[bl9->l2] + q3[lp]; 

  // ce nucleon (ici delta) devient un participant:
  jparticip[bl9->l2] = 1;
  if(verboseLevel > 3) {
    G4cout <<"Particle " << bl9->l2 << " is now participant." << G4endl;
  }
  
  if (ws->nosurf <= 0) {
    // surface
    pppp = std::sqrt(std::pow(bl1->p1[bl9->l2],2) + std::pow(bl1->p2[bl9->l2],2) + std::pow(bl1->p3[bl9->l2],2));
    rrrr = std::sqrt(std::pow(bl3->x1[bl9->l2],2) + std::pow(bl3->x2[bl9->l2],2) + std::pow(bl3->x3[bl9->l2],2));
    if (pppp <= bl10->pf) {
      xv = pppp/bl10->pf;
      rcorr = interpolateFunction(xv);
      if (rrrr > rcorr) {
	bl3->x1[bl9->l2] = bl3->x1[bl9->l2]*rcorr/rrrr;
	bl3->x2[bl9->l2] = bl3->x2[bl9->l2]*rcorr/rrrr;
	bl3->x3[bl9->l2] = bl3->x3[bl9->l2]*rcorr/rrrr;
      }
    }
    // fin surface
  }

  // difference with standard cascade :
  // the delta is located at the nucleon site to avoid problems
  // with the reflexion on the wall
  if(lp != npion) {
    lp1 = lp + 1;
    //    for(G4int i10 = lp1; i10 <= npion; lp1++) {
    for(G4int i10 = lp1; i10 <= npion; i10++) {
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

  if(bl2->k != 0) {
    kd = 0;
    ccr = tau;
    for(G4int i = 1; i <= bl2->k; i++) {
      i20 = i - kd;
      if((bl2->ind[i] == bl9->l1) || (bl2->ind[i] == bl9->l2) || (bl2->jnd[i] == bl9->l1) || (bl2->jnd[i] == bl9->l2)) {
	kd = kd + 1;
      }
      else {
	bl2->crois[i20] = bl2->crois[i] - ccr;
	bl2->ind[i20] = bl2->ind[i];
	bl2->jnd[i20] = bl2->jnd[i];
      }
    }

    bl2->k = bl2->k - kd;
    for(G4int i10 = 1; i10 <= bl2->k; i10++) {
      if(bl2->ind[i10] <= bl9->l1) {
	continue;
      }
      else {
	bl2->ind[i10] = bl2->ind[i10] - 1;
      }
    }
  }

  new1(bl9->l2);

  if(tdel <= bl4->tmax5) {
    bl2->k = bl2->k + 1;
    bl2->crois[bl2->k] = tdel;
    bl2->ind[bl2->k] = bl9->l2;
    bl2->jnd[bl2->k] = 0;
  }

  bl2->k = bl2->k + 1;
  bl2->crois[bl2->k] = ref(bl3->x1[bl9->l2], bl3->x2[bl9->l2], bl3->x3[bl9->l2], bl1->p1[bl9->l2], bl1->p2[bl9->l2], bl1->p3[bl9->l2], bl1->eps[bl9->l2],r22);
  if(verboseLevel > 3) {
    if(bl2->crois[bl2->k] < 0.0) {
      G4cout <<"G4Incl: Reflection time < 0! (line 3955)" << G4endl;
    }
  }
  bl2->ind[bl2->k] = bl9->l2;
  bl2->jnd[bl2->k] = -1;

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Going back to pnu449 at line 3917." << G4endl;
  }
  goto pnu449; // Line 3917

  // reflection on or transmission through the potential wall
 pnu600:
  //  G4cout <<"(checkpoint 'pnu600)" << G4endl;
  // deutons pas bien compris ici cv ?
  if (npproj[bl9->l1] == 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: npproj[l1] == 0. Going to pnu608." << G4endl;
    }
    goto pnu608;
  }
  
  if (bl1->ind1[bl9->l1] != 0) {
    if(verboseLevel > 3) {
      G4cout <<"wrong reentering particle (ind1[l1] != 0)" << G4endl;
      G4cout <<"ind1[" << bl9->l1 << "] = " << bl1->ind1[bl9->l1] << G4endl;
    }
  }
  
  if (bl3->x1[bl9->l1]*(bl1->p1[bl9->l1])+bl3->x2[bl9->l1]*(bl1->p2[bl9->l1])+bl3->x3[bl9->l1]*(bl1->p3[bl9->l1]) > 0.0) {
    if(verboseLevel > 3) {
      G4cout <<"wrong reentering particle" << G4endl;
      G4cout <<"particle: l1 = " << bl9->l1 << G4endl;
    }
  }
  
  if (bl1->ind2[bl9->l1] == 1) { //then
    npenter=npenter+1;
  } else {
    nnenter=nnenter+1;
  } // endif

  var_ab = std::pow(bl1->p1[bl9->l1],2) + std::pow(bl1->p2[bl9->l1],2) + std::pow(bl1->p3[bl9->l1],2);
  gpsg = 0.0;
  if (var_ab > 0.0) {
    gpsg = std::sqrt((std::pow(bl1->eps[bl9->l1]+v0,2)-pm2)/var_ab); 
  }
  
  bl1->p1[bl9->l1] = gpsg*(bl1->p1[bl9->l1]);                                                
  bl1->p2[bl9->l1] = gpsg*(bl1->p2[bl9->l1]);                                               
  bl1->p3[bl9->l1] = gpsg*(bl1->p3[bl9->l1]);                                               
  bl1->eps[bl9->l1] = bl1->eps[bl9->l1] + v0;
  npproj[bl9->l1] = 0;
  bl5->nesc[bl9->l1] = 0;

  // reevaluation of the times tab after entrance of 2nd,..nucleon 
  // of the projectile (goto 602 instead of 607 modif. 13/06/01)
  goto pnu602;

  // deutons
  // pour un non participant la transmission est impossible:
 pnu608:
  //  G4cout <<"(checkpoint 'pnu608)" << G4endl;
  if(varavat->kveux == 1) {
    varavat->del1avat[iavat] = bl1->ind1[bl9->l1];
    varavat->energyavat[iavat] = bl1->eps[bl9->l1] - fmp;
  }
  if(jparticip[bl9->l1] == 0) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: jparticip[l1] == 0. Going to pnu601." << G4endl;
    }
    goto pnu601;
  }
  if(varavat->kveux == 1) {
    varavat->go_out[iavat]=1;
  }
  if (bl1->ind1[bl9->l1] == 0) {
    goto pnu605;
  }
  fm = am(bl1->p1[bl9->l1],bl1->p2[bl9->l1],bl1->p3[bl9->l1],bl1->eps[bl9->l1]);
  pot = v1;
  goto pnu606;

 pnu605:
  //  G4cout <<"(checkpoint 'pnu605)" << G4endl;
  fm = fmp;
  pot = v0;

 pnu606:
  //  G4cout <<"(checkpoint 'pnu606)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu606. Calculating transmission probability." << G4endl;
  }
  tp = transmissionProb(bl1->eps[bl9->l1]-fm,bl1->ind2[bl9->l1],1,itch,bl3->r2,v0);
  if(varavat->kveux == 1) {
    varavat->energyavat[iavat] = bl1->eps[bl9->l1] - fm;
  }
  standardRandom(&rndm,&(hazard->igraine[10]));

  if (rndm > tp) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl:" << G4endl;
    }
    goto pnu601;
  }
  // Commande clusters ou non (1/0)! 
  if(kclst == 0) goto pnu1610;
  // Cluster goes here!
 pnu1610:
  //  G4cout <<"(checkpoint 'pnu1610)" << G4endl;
  // ici la particule l1 s'échappe du noyau:
  bl5->nesc[bl9->l1] = 1;
  nbquit = nbquit + 1;
  itch = itch - (1 + bl1->ind2[bl9->l1])/2;
  var_ab = std::pow(bl1->p1[bl9->l1],2) + std::pow(bl1->p2[bl9->l1],2) + std::pow(bl1->p3[bl9->l1],2);
  gpsg = 0.0;
  if(var_ab > 0.0) {
    gpsg = std::sqrt((std::pow(bl1->eps[bl9->l1]-pot,2) - fm*fm)/(var_ab));
  }
  bl1->p1[bl9->l1] = gpsg*(bl1->p1[bl9->l1]);
  bl1->p2[bl9->l1] = gpsg*(bl1->p2[bl9->l1]);
  bl1->p3[bl9->l1] = gpsg*(bl1->p3[bl9->l1]);
  bl1->eps[bl9->l1] = bl1->eps[bl9->l1] - pot;

  // Test for a transparent nucleon for composit beams and stored for FermiBreakup
  if(bl3->ia1 > 1 && nc[bl9->l1] == 0) {
    //    G4cout <<"(checkpoint 'incripspnu1610)" << G4endl;
    ips=ips+1;
    isPartOfSpectatorNucleus[bl9->l1] = true;
    ps->n_projspec[ips]=bl9->l1;
    nb_transprojo = nb_transprojo + 1;
  } // endif

  // pnu610: goes here

  // comptage des particules hors du noyau (7/6/2002):
  // (remnant minimum=1 nucleon)
  if(nbquit >= (ia-1)) {
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: nbquit >= (ia - 1). Going to pnu255." << G4endl;
    }
    goto pnu255;
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Going to pnu 602." << G4endl;
  }
  goto pnu602;

  // here no transmission possible
 pnu601:
  //  G4cout <<"(checkpoint 'pnu601)" << G4endl;
  pspr=bl3->x1[bl9->l1]*(bl1->p1[bl9->l1])+bl3->x2[bl9->l1]*(bl1->p2[bl9->l1])+bl3->x3[bl9->l1]*(bl1->p3[bl9->l1]);
  if(varavat->kveux == 1) {
    varavat->go_out[iavat]=0;
  }

  // surface: modif a.b. pour tenir compte du rayon variable du noyau.
  // (x2cour remplace r22 le rayon**2 fixe du noyau)
  x2cour = std::pow(bl3->x1[bl9->l1],2) + std::pow(bl3->x2[bl9->l1],2) + std::pow(bl3->x3[bl9->l1],2);
  bl1->p1[bl9->l1] = bl1->p1[bl9->l1] - 2.0*(bl3->x1[bl9->l1])*pspr/x2cour;
  bl1->p2[bl9->l1] = bl1->p2[bl9->l1] - 2.0*(bl3->x2[bl9->l1])*pspr/x2cour;
  bl1->p3[bl9->l1] = bl1->p3[bl9->l1] - 2.0*(bl3->x3[bl9->l1])*pspr/x2cour;
  // fin modif surface a.b.             

 pnu602:
  //  G4cout <<"(checkpoint 'pnu602)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu602." << G4endl;
  }
  
  if(bl2->k != 0) { // GOTO elimination?
    kd = 0;
    ccr = tau;
    for(G4int i = 1; i <= bl2->k; i++) {
      i20 = i - kd;

      if((bl2->jnd[i] == bl9->l1) || ((bl2->ind[i] == bl9->l1) && (bl2->jnd[i] != 0))) {
	kd = kd + 1;
	continue;
      }
      
      bl2->crois[i20] = bl2->crois[i] - ccr;
      bl2->ind[i20] = bl2->ind[i];
      bl2->jnd[i20] = bl2->jnd[i];
    }
    bl2->k = bl2->k - kd;

    if (bl5->nesc[bl9->l1] != 0) { // modif 10/02 pour logique des clusters (nesc()=2!)
      goto pnu613;
    }
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now calling new1(l1) (new1(" << bl9->l1 << "))" << G4endl;
  }
  new1(bl9->l1);

  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now calling ref." << G4endl;
    G4cout <<"x1 = " << bl3->x1[bl9->l1] <<" x2 = " << bl3->x2[bl9->l1] <<" x3 = " << bl3->x3[bl9->l1] << G4endl;
    G4cout <<"p1 = " << bl1->p1[bl9->l1] <<" p2 = " << bl1->p2[bl9->l1] <<" p3 = " << bl1->p3[bl9->l1] <<" eps = " << bl1->eps[bl9->l1] << G4endl;
  }
  tref = ref(bl3->x1[bl9->l1],bl3->x2[bl9->l1],bl3->x3[bl9->l1],bl1->p1[bl9->l1],bl1->p2[bl9->l1],bl1->p3[bl9->l1],bl1->eps[bl9->l1],r22); // line 4101
  if(verboseLevel > 3) {
    G4cout <<"Returned from function ref. tref = " << tref << G4endl;
  }
  
  if(verboseLevel > 3) {
    if(tref < 0.0) {
      G4cout <<"G4Incl: Reflection time < 0 (line 4101)!" << G4endl;
      G4cout <<"G4Incl: bl1->eps[" << bl9->l1 << "] = " << bl1->eps[bl9->l1] << G4endl;
    }
  }
  
  if (tref > bl4->tmax5) {
    goto pnu615;
  }
  bl2->k = bl2->k + 1;
  bl2->crois[bl2->k] = tref;
  bl2->ind[bl2->k] = bl9->l1;
  bl2->jnd[bl2->k] = -1;

 pnu615:
  //  G4cout <<"(checkpoint 'pnu615)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu615." << G4endl;
  }
  
  if (npion == 0) {
    goto pnu613;
  }
  if (bl1->ind1[bl9->l1] == 1) {
    goto pnu613;
  }
  for(G4int k20 = 1; k20 <= npion; k20++) { //do 614 k20=1,npion
    new3(y1[k20], y2[k20], y3[k20], q1[k20], q2[k20], q3[k20], q4[k20], k20, bl9->l1);
  }
 pnu613:
  //  G4cout <<"(checkpoint 'pnu613)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu613." << G4endl;
  }

  if (bl2->k == 0) {
    goto pnu230;
  }

  if(verboseLevel > 3) {
    G4cout <<"G4Incl. bl2->k != 0. Going back to pnu449 from line 4077." << G4endl;
    G4cout <<"G4Incl: bl2->k = " << bl2->k << G4endl;
    for(G4int myindex = 0; myindex <= bl2->k; myindex++) {
      G4cout <<"index = " << myindex << " ind = " << bl2->ind[myindex] << " jnd = " << bl2->jnd[myindex] << " crois = " << bl2->crois[myindex] << G4endl;
      if(bl2->crois[myindex] < 0.0) {
	G4cout <<"Negative time!!! Dumping information on collision: " << G4endl;
	G4int part1 = bl2->ind[myindex];
	G4int part2 = bl2->jnd[myindex];
	G4cout <<"particle 1 index = " << bl2->ind[myindex] << " \t particle 2 index = " << bl2->jnd[myindex] << G4endl;
	if(part1 >= 0) {
	  G4cout <<"Particle 1: " << G4endl;
	  G4cout <<"p1 = " << bl1->p1[part1] <<"p2 = " << bl1->p2[part1] <<"p3 = " << bl1->p3[part1] <<" eps = " << bl1->eps[part1] << G4endl;
	  G4cout <<"x1 = " << bl3->x1[part1] <<"x2 = " << bl3->x2[part1] <<"x3 = " << bl3->x3[part1] << G4endl;
	}
	if(part2 >= 0) {
	  G4cout <<"Particle 2: " << G4endl;
	  G4cout <<"p1 = " << bl1->p1[part2] <<"p2 = " << bl1->p2[part2] <<"p3 = " << bl1->p3[part2] <<" eps = " << bl1->eps[part2] << G4endl;
	  G4cout <<"x1 = " << bl3->x1[part2] <<"x2 = " << bl3->x2[part2] <<"x3 = " << bl3->x3[part2] << G4endl;
	}
      }
    }
  }

  goto pnu449; // Line 4077

  // decay of the surviving deltas
 pnu230:
  //  G4cout <<"(checkpoint 'pnu230)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu230." << G4endl;
  }
  // decay of the surviving deltas
pnu255:
  //  G4cout <<"(checkpoint 'pnu255)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu6255." << G4endl;
  }

  if (k3 == 1) {
    goto pnu256;
  }
  if (k4 == 0) {
    goto pnu256;
  }

  npidir = npion;

  idq = 0;
  destar = 0.;

  for(G4int i = 1; i <= ia; i++) {
    G4int iblcdpp=0;
    if (bl1->ind1[i] == 0) {
      continue;
    }
    npion = npion + 1;
    var_ab = std::pow(bl1->eps[i],2) - std::pow(bl1->p1[i],2) - std::pow(bl1->p2[i],2) - std::pow(bl1->p3[i],2);
    ym[npion] = 0.0;

    if(var_ab > 0.0) {
      ym[npion] = std::sqrt(var_ab);
    }
    xy1 = bl1->p1[i];
    xy2 = bl1->p2[i];
    xy3 = bl1->p3[i];
    xye = bl1->eps[i];
    ichd=bl1->ind2[i];
    if(varavat->kveux == 1) {
      iavat = iavat + 1;
      varavat->timeavat[iavat] = tim;
      varavat->l1avat[iavat] = i;
      varavat->l2avat[iavat] = -2;
      varavat->energyavat[iavat] = ym[npion];
      varavat->bloc_paul[iavat] = 0;
      varavat->bloc_cdpp[iavat] = 0;
      varavat->del1avat[iavat] = bl1->ind1[bl9->l1];
      varavat->jpartl1[iavat] = 1;
      varavat->jpartl2[iavat] = 0;
    }

    decay2(&(bl1->p1[i]), &(bl1->p2[i]), &(bl1->p3[i]), &(bl1->eps[i]), &(q1[npion]), &(q2[npion]), &(q3[npion]),
	   &(q4[npion]), &(ym[npion]), &fmp, &fmpi, &(bl9->hel[i]));

    if(verboseLevel > 3) {
      G4cout <<"Quantities after decay2: " << G4endl;
      G4cout <<"i = " << i << " bl1->p1[i] = " << bl1->p1[i] << " bl1->p2[i] = " << bl1->p2[i] << " bl1->p3[i] = " << bl1->p3[i] << " bl1->eps[i] = " << bl1->eps[i] << G4endl;
    }

// C 6/08/09; Don't have to check CDPP for the forced decay of a DELTA OUTSIDE!
//       IF(NESC(I).EQ.1) GO TO 1850
      
// C------
// C On teste CDPP apres la decroissance forcee	June/2005 (AB-CV)
// C	Si rejet de la decroissance -> delta devient nucleon sans pion 
// C       et toute l'energie en exces est mise en energie d'excitation du noyau
// C       On respecte aussi la conservation de charge (Z du remnant!)
 
// C Le decay ne peut conduire a un noyau de A nucleons
// C  sous l'energie de Fermi et dans une config. d'energie inferieure a
// C  EFER-(IA2-NBALTTF)*TF).
//        EGS=0.
//        NBALTTF=0
//        DO II=1,IA
//          IF(NESC(II).EQ.0) THEN
// 	    IF(SQRT(P1(II)**2+P2(II)**2+P3(II)**2).LT.PF) THEN
// 	       NBALTTF=NBALTTF+1
// 	       EGS=EGS+EPS(II)-FMP
// 	    ENDIF
// 	 ENDIF
//        END DO
//        IF(EGS.GE.(EFER-(IA2-NBALTTF)*TF)) GO TO 1850
// C ATTENTION, logique negative!!! Liberer le goto si on veut supprimer la
// C   sequence precedente (CDPP sur Delta-> pi N)
// C        GO TO 1850

// C Ici ce decay est interdit par CDPP.
// C Restitution du Delta: un neutron au repos et correction conservee pour 
// C la charge IDQ et l'energie d'excitation DESTAR du remnant. On viole donc
// C un peu la conservation d'impulsion, mais semble moins grave que de convertir
// C la masse d'un pion en impulsion d'un nucleon unique....
//       IBLCDPP=1		!flag de blocage par CDPP
//       IND1(i) = 0	!identif delta (nucleon)
//       IND2(i) = -1	!isospin du neutron
// 	IDQ = IDQ + (ICHD+1)/2
//       p1(i)=0.		!impulsion energie
//       p2(i)=0.
//       p3(i)=0.
//       eps(i)=FMP
// 	DESTAR = DESTAR + xye - FMP
//       NPION=NPION-1	!pas le dernier pion

// 1850  CONTINUE
// C------    
    if(bl5->nesc[i] == 0) {
      idecf = 1;
    }

    if (ws->nosurf <= 0) {
      // surface
      if (bl5->nesc[i] == 0.0) {
	pppp = std::sqrt(std::pow(bl1->p1[i],2) + std::pow(bl1->p2[i],2) + std::pow(bl1->p3[i],2));
	rrrr = std::sqrt(std::pow(bl3->x1[i],2) + std::pow(bl3->x2[i],2) + std::pow(bl3->x3[i],2));
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
    if(iblcdpp == 1) continue; // goto pnu257;

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
    //    G4cout <<"(checkpoint 'pnu283)" << G4endl;
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: Now at pnu283." << G4endl;
    }

    ipi[npion] = bl1->ind2[i]*2;
    bl1->ind2[i] = -1*(bl1->ind2[i]);
    goto pnu285;

  pnu280:
    //    G4cout <<"(checkpoint 'pnu285)" << G4endl;
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: Now at pnu280." << G4endl;
    }

    bl1->ind2[i] = bl1->ind2[i]/3;
    ipi[npion] = 2*(bl1->ind2[i]);

  pnu285:
    //    G4cout <<"(checkpoint 'pnu285)" << G4endl;
    if(verboseLevel > 3) {
      G4cout <<"G4Incl: Now at pnu285." << G4endl;
    }

    y1[npion] = bl3->x1[i];
    y2[npion] = bl3->x2[i];
    y3[npion] = bl3->x3[i];
  }

  if(verboseLevel > 3) {
    G4cout <<"out of the loop..." << G4endl;
  }
  
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
    for(G4int i = 1; i <= ia; i++) {
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
  //  G4cout <<"(checkpoint 'pnu256)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu256." << G4endl;
  }

  elead = 0.0;
  lead = 0;
  npx = 0;
  erem = 0.;
  // Excitation energy of the final delta of blocked decay: (AB CV 6/2005)
  erem = erem + destar;
  izrem = 0;
  // Correct the charge of the remnant if the final delta is blocked: (A.B., C.V. 6/2005)
  izrem = izrem + idq;
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
  p1spec=0.;
  p2spec=0.;
  p3spec=0.;
  cmultn = 0.0;

  if (kindstruct->kindf7 <= 2 && kindstruct->kindf7 > 0) {
    if (ncol == 0 || nc[1] == 0) { // then nc(1)->nc[0]
      if(verboseLevel > 3) {
	G4cout <<"no collisioms" << G4endl;
      }
      goto pnu9100;
    }
  }
  else {
    if (kindstruct->kindf7 <= 5 && kindstruct->kindf7 > 0) {
      if (ncol == 0) {
	if(verboseLevel > 3) {
	  G4cout <<"no collisioms" << G4endl;
	}
	goto pnu9100;
      }
    }
    else {
      // ici faisceau composite: modif a.b. 2/2002 pour tous les composites:
      nsum_col = 0;
      for(G4int i = 1; i <= bl3->ia1; i++) {
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
  //  G4cout <<"(checkpoint 'pnu9100)" << G4endl;
  iarem = bl3->ia2;
  izrem = iz2;
  esrem = 0.0;
  erecrem = 0.0;
  nopart = -1;
  // fin ajout cv
  if(verboseLevel > 3) {
    G4cout <<"End of algorithm after pnu9100." << G4endl;
  }
  goto pnureturn;

 pnu9101:
  //  G4cout <<"(checkpoint 'pnu9101)" << G4endl;
  if(verboseLevel > 3) {
    G4cout <<"G4Incl: Now at pnu9101." << G4endl;
  }

  nopart = 0;
  ekout = 0.0;

  if(kclst != 0) { // Clusters go here!
  }

// Treat here the projectile spectators (including transparents):
  projo_spec(bl3->ia1,ips,fmpinc,pinc,tlab);

  for(G4int i = 1; i <= ia; i++) { // do 258 i=1,ia
    //  pnu259:
    if (bl5->nesc[i] == 0) goto pnu254;
    if (bl5->nesc[i] > 1) continue; // on evite les nucleons des clusters

    // here nesc=1 or negative value (emitted nucleons, spectators included)      
    xl1_i=bl3->x2[i]*bl1->p3[i]-bl3->x3[i]*bl1->p2[i];	// moment angulaire emporte
    xl2_i=bl3->x3[i]*bl1->p1[i]-bl3->x1[i]*bl1->p3[i];
    xl3_i=bl3->x1[i]*bl1->p2[i]-bl3->x2[i]*bl1->p1[i];

    xl1=xl1-xl1_i;		                                      
    xl2=xl2-xl2_i;		                                       
    xl3=xl3-xl3_i;		                                       

    //   ici ajout de pout cv le 5/7/95
    pout1=pout1+bl1->p1[i];
    pout2=pout2+bl1->p2[i];
    pout3=pout3+bl1->p3[i];
    //      write(6,*)'eout nucleon',eout
    eout = eout + bl1->eps[i]-fmp;

    ic33=(bl1->ind2[i]+3)/2;

    // here we don't copy the spectators of the projectile bound in a
    //   spectator nucleus (treated by the deexcitation module):
    if(kindstruct->kindf7 < 0 && i <= bl3->ia1 && ps->a_projspec != 0) {
      for(G4int iloc =  1; iloc <= ps->a_projspec; iloc++) { //do iloc=1,a_projspec
	if(i == ps->n_projspec[iloc]) {
	  p1spec=p1spec+bl1->p1[i];
	  p2spec=p2spec+bl1->p2[i];
	  p3spec=p3spec+bl1->p3[i];
	  // for balance, a transparent from projectile is out:
	  if(i > ips-nb_transprojo) ekout=ekout+bl1->eps[i]-fmp;
	  goto pnu258_end_of_outer_loop;
	  //	  continue;
	} // endif
      } // enddo
    } // endif
    //    if (isPartOfSpectatorNucleus[i]) goto pnu258_end_of_outer_loop;

    nopart=nopart+1;

    ekout=ekout+bl1->eps[i]-fmp;

    kind[nopart]=3-ic33;
    // spectators of composite projectiles (7/2006, ab)
    if(npproj[i] == 1) kind[nopart] = -kind[nopart];

    ep[nopart]=bl1->eps[i]-fmp;
    bmass[nopart]=fmp;
    ptotl=std::sqrt(std::pow(bl1->p1[i],2) + std::pow(bl1->p2[i],2)+std::pow(bl1->p3[i],2));
    alpha[nopart]=bl1->p1[i]/ptotl;
    beta[nopart]=bl1->p2[i]/ptotl;
    gam[nopart]=bl1->p3[i]/ptotl;
 
    continue;
  pnu254:
    t[4]=bl3->x1[i]*bl3->x1[i]+bl3->x2[i]*bl3->x2[i]+bl3->x3[i]*bl3->x3[i];
    erem=erem+bl1->eps[i]-fmp;
    rcm1=rcm1+bl3->x1[i];
    rcm2=rcm2+bl3->x2[i];
    rcm3=rcm3+bl3->x3[i];
    prem1=prem1+bl1->p1[i];
    prem2=prem2+bl1->p2[i];
    prem3=prem3+bl1->p3[i];
    izrem=izrem+(1+bl1->ind2[i])/2;
    iarem=iarem+1;
  pnu258_end_of_outer_loop:
    continue;
  }
  //  258 continue  

  //pnu258:
  //  G4cout <<"(checkpoint 'pnu258)" << G4endl;
  //  correction pions 21/3/95 jc
  ichpion = 0;
  if(npion != 0) {
    for(G4int ipion = 1; ipion <= npion; ipion++) {
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
      ptotl = std::sqrt(std::pow(q1[ipion],2) + std::pow(q2[ipion],2) + std::pow(q3[ipion],2));
      ep[nopart] = q4[ipion] - fmpi;
      bmass[nopart] = fmpi;
      ekout = ekout + ep[nopart];
      alpha[nopart] = q1[ipion]/ptotl;
      beta[nopart] = q2[ipion]/ptotl;
      gam[nopart] = q3[ipion]/ptotl;
      //      nopart++;
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
  l = int(std::sqrt(xl1*xl1 + xl2*xl2 + xl3*xl3)/hc + 0.5);

  iej = bl3->ia2 - irem;

  eh5 = erem - std::pow(double(irem)/double(a2),1.666667)*efer;
  if(verboseLevel > 3) {
    G4cout <<"erem used for excitation energy calculation = " << erem << G4endl;
    G4cout <<"irem used for excitation energy calculation = " << irem << G4endl;
    G4cout <<"a2 used for excitation energy calculation = " << a2 << G4endl;
    G4cout <<"eh5 used for excitation energy calculation = " << eh5 << G4endl;
  }
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
	if(verboseLevel > 3) {
	  G4cout <<"npion == 0" << G4endl;
	  G4cout <<"End of algorithm because npion == 0" << G4endl;
	}
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
  pfreml2 = std::pow(pfrem1,2) + std::pow(pfrem2,2) + std::pow(pfrem3,2);
  if (pfreml2 > 1.0e-12) {
    pfreml = std::sqrt(pfreml2);
    alrem = pfrem1/pfreml;
    berem = pfrem2/pfreml;
    garem = pfrem3/pfreml;
  }
  else {
    alrem = 0.0;
    berem = 0.0;
    garem = 1.0;
  }

  erecrem = pfreml2/(std::sqrt(pfreml2 + std::pow((fmp*iarem),2)) + fmp*iarem);

  if(iarem == 1) {
    erecrem = erecrem + eh5;
  }

  // correction recul
  erecg = erecrem + ecoreh5;
  // correction energie d'excitation pour une absorption (a.b., c.v. 2/2002)
  esrem = eh5;

  if (ekout < 0.001) {
    if(verboseLevel > 3) {
      G4cout <<"ekout < 0.001" << G4endl;
      G4cout <<"End of algorithm because kout < 0.001" << G4endl;
    }
    goto pnureturn;
  }
 
  // on ote l'instruction precedente car esrem toujours nulle 14/9/99

  if (erecg > 0.25) {
    fffc = (ekout - erecg)/ekout;
    if (fffc < 0.0) {
      fffc = 0.0;
    }

    for(G4int ipart = 0; ipart <= nopart; ipart++) {
      ep[ipart] = ep[ipart]*fffc;
    }
  }

  // modif boudard juillet 99 (il faut tenir compte de la renormalisation
  // des energies pour les impulsions.)
  pfrem1= -p1spec;
  pfrem2= -p2spec;
  pfrem3=pinc-p3spec;
//   if(useProjSpect) {
//     pfrem1 = -p1_spec;
//     pfrem2 = -p2_spec;
//     pfrem3 = pinc - p3_spec;
//   } else {
//     pfrem1 = 0.0;
//     pfrem2 = 0.0;
//     pfrem3 = pinc;
//   }
  for(G4int ipart = 0; ipart <= nopart; ipart++) {
    G4double xmodp = 0.0;
    if(ep[ipart] < 0.0) { // Safeguard against particles with negative energy
      xmodp = 0.0;
      G4cout <<"The energy of particle " << ipart << " is negative (" << ep[ipart] << "). Forcing xmodp = " << xmodp << "." << G4endl;
    } else {
      xmodp = std::sqrt(ep[ipart]*(2.0*bmass[ipart] + ep[ipart]));
    }
    pfrem1 = pfrem1 - alpha[ipart]*xmodp;
    pfrem2 = pfrem2 - beta[ipart]*xmodp;
    pfrem3 = pfrem3 - gam[ipart]*xmodp;
  }
  // fin modif a.b.

  pfreml2 = std::pow(pfrem1,2) + std::pow(pfrem2,2) + std::pow(pfrem3,2);
  erecrem = pfreml2/(std::sqrt(pfreml2 + std::pow((fmp*iarem),2)) + fmp*iarem);

  if (pfreml2 > 1.0e-12) {
    pfreml = std::sqrt(pfreml2);
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
  goto pnureturn;

  if(verboseLevel > 3) {
    G4cout <<"Reached end of routine..." << G4endl;
  }
  goto pnureturn;

  if(verboseLevel > 3) {
    G4cout <<"ia1 > 1 ! " << G4endl;
  }
 pnureturn:
  //  G4cout <<"(checkpoint 'pnureturn)" << G4endl;
//   if(useProjSpect) {
//     varntp->spectatorA = a_projspec;
//     varntp->spectatorZ = z_projspec;
//     varntp->spectatorEx = ex_projspec;
//     varntp->spectatorT = t_projspec;
//     varntp->spectatorM = m_projspec;
//     varntp->spectatorP1 = p1_projspec;
//     varntp->spectatorP2 = p1_projspec;
//     varntp->spectatorP3 = p1_projspec;
//   }
  (*ibert_p) = ibert;
  (*nopart_p) = nopart;
  (*izrem_p) = izrem;
  (*iarem_p) = iarem;
  (*esrem_p) = esrem; 
  (*erecrem_p) = erecrem;
  (*alrem_p) = alrem;
  (*berem_p) = berem;
  (*garem_p) = garem;
  (*bimpact_p) = bimpact;
  (*l_p) = l;
  (*xjrem) = 0;
  (*yjrem) = 0;
  (*zjrem) = 0;
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

  G4double debugOutput = 0.0;
  debugOutput = am(p1,p2,p3,e1);

  G4double pout11 = (*pout11_p);
  G4double pout12 = (*pout12_p);
  G4double pout13 = (*pout13_p);
  G4double eout1 = (*eout1_p);

  G4double q1 = (*q1_p);
  G4double q2 = (*q2_p);
  G4double q3 = (*q3_p);
  //  G4double q4 = -1*(*q4_p);
  G4double q4 = (*q4_p);

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
  G4double a = 0.0;
  G4double aaa = 0.0;
  G4double aac = 0.0;
  //  G4double alog;
  G4double alphac = 0.0;
  //  G4double amax1;
  G4double apt = 0.0;
  G4double argu = 0.0;
  G4double b = 0.0;
  G4double btmax = 0.0;
  G4double cfi = 0.0;
  G4double cpt = 0.0;
  G4double ctet = 0.0;
  G4double e3 = 0.0;
  G4double ecm = 0.0;
  G4double ex[3];
  G4double ey[3];
  G4double ez[3];
  G4double qq[3];
  for(G4int init_i = 0; init_i < 3; init_i++) {
    ex[init_i] = 0.0;
    ey[init_i] = 0.0;
    ez[init_i] = 0.0;
    qq[init_i] = 0.0;
  }
  G4double f3 = 0.0;
  G4double f3max = 0.0;
  G4double fi = 0.0;
  G4double fracpn = 0.0;
  G4double heli = 0.0;
  G4int iexpi = 0;
  G4int ii = 0;
  G4int index = 0;
  G4int index2 = 0;
  G4int isi = 0;
  G4double pin = 0.0;
  G4double pl = 0.0;
  G4double pnorm = 0.0;
  //  G4double pq = 0.0;
  G4double psq = 0.0;
  G4double qq4 = 0.0;
  G4double ranres = 0.0;
  G4double rndm = 0.0;
  G4double s = 0.0;
  G4double s1 = 0.0;
  //  G4double sel;
  G4double sfi = 0.0;
  G4double stet = 0.0;
  G4double t = 0.0;
  G4double x = 0.0;
  G4double xkh = 0.0;
  G4double xp1 = 0.0;
  G4double xp2 = 0.0;
  G4double xp3 = 0.0;
  G4double xx = 0.0;
  G4double y = 0.0;
  G4double yn = 0.0;
  G4double z = 0.0;
  G4double zn = 0.0;
  G4double zz = 0.0;

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
  G4double xmdel = 1232.0;
  G4double xpi = 138.0;

  //   2848	c      data iy1,iy2,iy3,iy4,iy5,iy6,iy7,iy8,iy10,iy11,iy12,iy13/         p-n17840
  //   2849	c     1 12345,22345,32345,42345,52345,62345,72345,82345,34567,47059,21033p-n17850
  //   2850	c     1 12345,22345,32345,42345,52345,62345,72345,82341,34567,47059,21033p-n17850
  //   2851	c     2,32835/                                                           p-n17860
  //   2852	c      data iy9/15637/
  //  2853	      pcm(e,a,c)=0.5*std::sqrt((e**2-(a+c)**2)*(e**2-(a-c)**2))/e            p-n17870
  G4int iso = is1 + is2;
  np = 0;
  psq = p1*p1 + p2*p2 + p3*p3;
  pnorm = std::sqrt(psq);
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
    bl6->xx10=std::sqrt(e1*e1-psq);
    bl6->isa=is1;
  }
  else {
    bl6->xx10 = std::sqrt(std::pow(eout1,2)-psq);
    bl6->isa = is2;
  }

  s = s1 + srec(ecm,bl6->xx10, iso,int(bl6->isa));
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
  pl = 0.5*ecm*std::sqrt(std::pow(ecm,2) - 4.0*xm2)/xm;
  x = 0.001*pl;
  
  if (iso == 0) {
    goto collis80;
  }
  if (pl > 2000.0) {
    goto collis81;
  }
  x = x*x;
  x = std::pow(x,4);
  b = 5.5e-6*x/(7.7+x);
  goto collis82;

 collis81:
  b = (5.34 + 0.67*(x-2.0))*1.e-6;
  goto collis82;

 collis80: 
  if (pl < 800.) {
    b = (7.16 - 1.63*x)*1.e-6;
    b = b/(1.0 + std::exp(-(x - 0.45)/0.05));
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
  debugOutput = am(p1,p2,p3,e1);
  
  btmax = 4.0*psq*b;
  z = std::exp(-btmax);
  standardRandom(&rndm, &(hazard->igraine[1]));
  ranres = rndm; 
  y = 1.0 - rndm*(1.0 - z);
  t = std::log(y)/b;
  iexpi = 0;

  if (((m1+m2) == 0) && (iso == 0)) {
    apt = 1.0;                                                       
    if (pl > 800.) {
      apt = std::pow((800.0/pl),2);                                          
      cpt = amax1(6.23*std::exp(-1.79*x),0.3);
      alphac = 100.0*1.e-6;
      aaa = (1 + apt)*(1 - std::exp(-btmax))/b;
      argu = psq*alphac;

      if (argu >= 8) {
	argu = 0.0;
      }
      else {
	argu = std::exp(-4.0*argu);
      }

      aac = cpt*(1.0 - argu)/alphac;
      fracpn = aaa/(aac + aaa);
      standardRandom(&rndm, &(hazard->igraine[7]));
      if (rndm > fracpn) {
	z = std::exp(-4.0*psq*alphac);
	iexpi = 1;
	//	y = 1.0 - ranres*(10.0 - z);
	y = 1.0 - ranres*(1.0 - z);
	t = std::log(y)/alphac;
      }
    }
  }

  ctet = 1.0 + 0.5*t/psq;
  if(std::fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }      

  stet = std::sqrt(1.0 - std::pow(ctet,2));
  standardRandom(&rndm, &(hazard->igraine[2]));
  fi = 6.2832*rndm;
  cfi = std::cos(fi);
  sfi = std::sin(fi);
  xx = p1*p1 + p2*p2;
  zz = p3*p3;    

  debugOutput = am(p1,p2,p3,e1);

  if(xx >= (zz*1.0e-8)) {
    yn=std::sqrt(xx);
    zn=yn*pnorm;
    ez[0] = p1/pnorm; // ez(1) -> ez[0] and so on...
    ez[1] = p2/pnorm;
    ez[2] = p3/pnorm;
    ex[0] = p2/yn;
    ex[1] = -p1/yn;
    ex[2] = 0.0;     
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
  debugOutput = am(p1,p2,p3,e1);

  // backward scattering according the parametrization of ref
  // prc56(1997)1
  if (m1+m2 == 1) goto collis133;
  if (iso != 0) goto collis133;
  standardRandom(&rndm,&(hazard->igraine[7]));
  apt = 1.0;
  if (pl > 800.0) {
    apt = std::pow(800.0/pl,2);
  } //endif
  if (iexpi == 1 || rndm > 1.0/(1.0+apt)) { // then
    ii = is1;
    is1 = is2;
    is2 = ii;
  } // endif
 collis133:

  debugOutput = am(p1,p2,p3,e1);
  goto exitRoutine;

  // delta production
  // the production is not isotropic in this version
  // it has the same std::exp(b*t) structure as the nn elastic scatteringp-n19170
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

  y = std::tan(3.1415926*(rndm - 0.5));
  x = 1232.0 + 0.5*130.0*y + bl8->ramass;
  if (x < (xm+xpi+2.0)) {
    goto collis101;
  }
  if (ecm < (x+xm+1.)) {
    goto collis101;
  }

  // generation of the delta mass with the penetration factor
  // (see prc56(1997)2431)
  y = std::pow(ecm,2);
  q2 = (y - std::pow(1076.0,2))*(y - std::pow(800.0,2))/y/4.0;                                 
  q3 = std::pow((std::sqrt(q2)),3);                                                  
  f3max = q3/(q3 + std::pow(180.0,3));                                      
  y = std::pow(x,2);
  q2 = (y - std::pow(1076.0,2))*(y - std::pow(800.0,2))/y/4.0;                                 
  q3 = std::pow((std::sqrt(q2)),3);                                                 
  f3 = q3/(q3 + std::pow(180.0,3));

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
  x = 0.001*0.5*ecm*std::sqrt(std::pow(ecm,2) - 4.0*xm2)/xm;
  if(x < 1.4) {
    b = (5.287/(1.0 + std::exp((1.3 - x)/0.05)))*1.e-6;
  }
  else {
    b = (4.65 + 0.706*(x - 1.4))*1.e-6;
  }

  xkh = 2.0*b*pin*pnorm;
  ctet=1.0 + std::log(1.0 - rndm*(1.0 - std::exp(-2.0*xkh)))/xkh;
  if(std::fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }
  
  stet = std::sqrt(1.0 - std::pow(ctet,2));
  standardRandom(&rndm, &(hazard->igraine[6]));
  fi = 6.2832*rndm;
  cfi = std::cos(fi);
  sfi = std::sin(fi);

  // delta production: correction of the angular distribution 02/09/02
  xx = p1*p1 + p2*p2;
  zz = p3*p3;
  if(xx >= (zz*1.0e-8)) {
    yn = std::sqrt(xx);
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

  e3 = std::sqrt(xp1*xp1 + xp2*xp2 + xp3*xp3 + xm*xm);
  if(k4 != 0) {
    goto collis161;
  }
  
  // decay of the delta particle (k4=0)
  np = 1;
  ip = 0;
  qq[0] = xp1; //qq(1) -> qq[0]
  qq[1] = xp2;
  qq[2] = xp3;
  qq4 = std::sqrt(xp1*xp1 + xp2*xp2 + xp3*xp3 + xmdel*xmdel);
  heli = std::pow(ctet,2);

  if(verboseLevel > 3) {
    G4cout <<"Caling decay2 from collis" << G4endl;
  }
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
  debugOutput = am(p1,p2,p3,e1);

  if (iso == 0) {
    goto collis150;
  }
  if (rndm > 0.333333) {
    debugOutput = am(p1,p2,p3,e1);
    goto exitRoutine;
  }

  is1 = -is1;
  ip = -2*is1;

 collis150:
  if (index == 1) {
    debugOutput = am(p1,p2,p3,e1);
    goto exitRoutine;
  }
  if (rndm < 0.5) {
    goto collis152;
  }
  is1 = 1;
  is2 = 1;
  ip = -2;
  debugOutput = am(p1,p2,p3,e1);
  goto exitRoutine;

 collis152: 
  is1 = -1;
  is2 = -1;
  ip = 2;
  debugOutput = am(p1,p2,p3,e1);
  goto exitRoutine;

 collis160:
  pout11 = -p1;
  pout12 = -p2;
  pout13 = -p3;
  debugOutput = am(p1,p2,p3,e1);
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
  debugOutput = am(p1,p2,p3,e1);

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

  bl9->hel[bl9->l1] = std::pow(ctet,2);
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
  if(std::fabs(ctet) > 1.0) {
    ctet = sign(1.0,ctet);
  }
  stet = std::sqrt(1.0 - ctet*ctet);
  standardRandom(&rndm, &(hazard->igraine[12]));
  fi = 6.2832*rndm;
  cfi = std::cos(fi);
  sfi = std::sin(fi);
  p1 = pnorm*stet*cfi;
  p2 = pnorm*stet*sfi;
  p3 = pnorm*ctet;
  m1 = 0;
  m2 = 0;
  e1 = std::sqrt(p1*p1 + p2*p2 + p3*p3 + xm*xm);
  eout1 = ecm - e1;
  debugOutput = am(p1,p2,p3,e1);

  if (iso == 0) {
    goto collis160;
  }
  is1=iso/2;
  is2=iso/2;
  goto collis160;
  
 exitRoutine:
  debugOutput = am(p1,p2,p3,e1);
  (*p1_p) = p1;// Was pq
  (*p2_p) = p2;
  (*p3_p) = p3;

  (*e1_p) = e1;

  debugOutput = am(pout11,pout12,pout13,eout1);
  (*pout11_p) = pout11;
  (*pout12_p) = pout12;
  (*pout13_p) = pout13;
  (*eout1_p) = eout1;

  (*q1_p) = q1;
  (*q2_p) = q2;
  (*q3_p) = q3;
  (*q4_p) = q4;

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
  // the anisotropy is supposed to follow a 1+3*hel*(std::cos(theta))**2
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

  G4double rndm = 0.0;

  G4double xe = wp;
  G4double b1 = p1/xe;
  G4double b2 = p2/xe;
  G4double b3 = p3/xe;                                                          
  // PK: NaN workaround
  // if(((std::pow(xi,2)-std::pow(x1+x2,2))*(std::pow(xi,2)-std::pow(x1-x2,2))) < 0) {
  //   xi = xi+x2+100.0;
  // }
  // PK
  G4double xq = pcm(xi,x1,x2);                                                  
  G4double ctet = 0.0, stet = 0.0;

  G4double fi = 0.0, cfi = 0.0, sfi = 0.0;
  G4double sal = 0.0, cal = 0.0;
  G4double t1 = 0.0, t2 = 0.0;
  G4double w1 = 0.0;
  G4double beta = 0.0;

  if(verboseLevel > 3) {
    G4cout <<"Delta decay in progress: " << G4endl;
    G4cout <<"Starting values: " << G4endl;
    G4cout <<"p1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << " wp = " << wp << G4endl;
    G4cout <<"q1 = " << q1 << " q2 = " << q2 << " q3 = " << q3 << " wq = " << wq << G4endl;
  }

 decay2100:
  standardRandom(&rndm,&(hazard->igraine[7]));
  ctet = -1.0 + 2.0*rndm;
  if(std::abs(ctet) > 1.0) ctet = sign(1.0,ctet);
  stet = std::sqrt(1.0 - std::pow(ctet, 2));
  standardRandom(&rndm,&(hazard->igraine[9]));
  if (rndm > ((1.0 + 3.0 * hel * std::pow(ctet,2))/(1.0 + 3.0*hel))) goto decay2100;
  standardRandom(&rndm,&(hazard->igraine[8]));
  fi = 6.2832*rndm;
  cfi = std::cos(fi);
  sfi = std::sin(fi);
  beta = std::sqrt(b1*b1 + b2*b2 + b3*b3);
  if (beta < 1.0e-10) goto decay2101;
  sal = std::sqrt(std::pow(b1, 2) + std::pow(b2, 2))/beta;
  cal = b3/beta;
  if (sal < 1.0e-6) goto decay2101;
  t1 = ctet + cal*stet*sfi/sal;
  t2 = stet/sal;
  q1 = xq*(b1*t1 + b2*t2*cfi)/beta;
  q2 = xq*(b2*t1 - b1*t2*cfi)/beta;
  q3 = xq*(b3*t1/beta - t2*sfi);
  goto decay2102;
 decay2101:
  q1 = xq * stet*cfi;
  q2 = xq * stet*sfi;
  q3 = xq * ctet;
 decay2102:
  hel = 0.0;                                                       
  w1 = q1*q1 + q2*q2 + q3*q3;
  wq = std::sqrt(w1 + x2*x2);
  p1 = -q1;
  p2 = -q2;
  p3 = -q3;
  wp = std::sqrt(w1+x1*x1);
  loren(&q1, &q2, &q3, &b1, &b2, &b3, &wq);
  loren(&p1, &p2, &p3, &b1, &b2, &b3, &wp);

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
  G4double t[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
    bl1->ta = -1.0*t[6]/t[9];
  }

  bl3->rab2 = t[3]*t[3] + t[4]*t[4] + t[5]*t[5] + bl1->ta*t[6];
}

void G4Incl::newt(G4int l1, G4int l2)
{
  G4int ig = 0, id = 0, kg = 0, kd = 0;
  G4int iy = 0, ix = 0;
  G4double E = 0.0;

  G4int ia = bl3->ia1 + bl3->ia2;
  for(G4int i = 1; i <= ia; i++) { // do 52 i=1,ia
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
    if (bl1->ta < 0.0) {
      goto newt50;
    }
    if(bl1->ta > bl4->tmax5) {
      goto newt50;
    }
    if (bl1->ta < bl5->tlg[l2]) { // tlg(12)->tlg[11]
      goto newt50;
    }
    if ((bl1->ind1[ig]+bl1->ind1[id]) > 0) {
      goto newt60;
    }

    E=am(bl1->p1[ig]+bl1->p1[id],bl1->p2[ig]+bl1->p2[id],bl1->p3[ig]+bl1->p3[id],bl1->eps[ig]+bl1->eps[id]);
    if (E < 1925.0) {
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
  G4int ia = 0, iy = 0, ix = 0;
  G4double E = 0.0;

  ia=bl3->ia1+bl3->ia2;
  for(G4int i = 1; i <= ia; i++) {
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
    if(bl1->ta < 0.0) {
      continue;
    }
    if(bl1->ta > bl4->tmax5) {
      continue;
    }
    if (bl1->ind1[i]+bl1->ind1[l1] > 0) {
      goto new160;
    }
    E=am(bl1->p1[i]+bl1->p1[l1],bl1->p2[i]+bl1->p2[l1],bl1->p3[i]+bl1->p3[l1],bl1->eps[i]+bl1->eps[l1]);
    if (E < 1925.0) {
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
    if (E < 1925.0) {
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
    if ((31.0*(bl3->rab2)) > totalCrossSection(E,iy,bl1->ind2[i]+bl1->ind2[l1])) {
      continue;
    }
  new170: 
    bl2->k=bl2->k+1;
    bl2->crois[bl2->k]=bl1->ta;
    bl2->ind[bl2->k]=i;
    bl2->jnd[bl2->k]=l1;
  }
  //  std::ofstream newout("new1Dump.out");
  //  dumpBl2(newout);
  //  dumpBl5(newout);
  //  newout.close();
  //  exit(0);
}

void G4Incl::new2(G4double y1, G4double y2, G4double y3, G4double q1, G4double q2, G4double q3, 
		  G4double q4, G4int npion, G4int l1)
{
  G4double t[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  G4int ia = bl3->ia1 + bl3->ia2;
  for(G4int i = 1; i <= ia; i++) {
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
    G4double E = std::sqrt(std::pow((bl1->eps[i]+q4),2) - std::pow((bl1->p1[i]+q1),2) - std::pow((bl1->p2[i]+q2),2) - std::pow((bl1->p3[i]+q3),2));
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
  G4double t[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double E = 0.0, xx2 = 0.0;
  G4int ia = 0;

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
  E = std::sqrt(std::pow((bl1->eps[l1]+q4),2) - std::pow((bl1->p1[l1]+q1),2) - std::pow((bl1->p2[l1]+q2),2) - std::pow((bl1->p3[l1]+q3),2));
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
  G4double gam = std::sqrt(gam2);
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

  G4double pmod = 0.0, pr2 = 0.0;
  G4double xr2 = 0.0, rdeq = 0.0, dx2 = 0.0, dp2 = 0.0;
  G4double rs = 0.0, vol = 0.0;
  G4int nl = 0;
  G4int ia = 0;

  if (ws->npaulstr == 2) {
    return 0.0;
  }

  if (ws->npaulstr == 1) {
    // pauli strict
    pmod = std::sqrt(std::pow(bl1->p1[l],2) + std::pow(bl1->p2[l],2) + std::pow(bl1->p3[l],2));
    if (pmod < fermi->pf) {
      return 1.0;
    } else {
      return 0.0;
    }
  } else {
    // Statistic Pauli blocking
    xr2 = xr*xr;
    pr2 = pr*pr;
    vol = std::pow((40.0*3.1415926/3.0),2) * (std::pow((xr*pr)/(2.0*3.1415926*197.33),3));
    rs = std::sqrt(bl3->x1[l]*bl3->x1[l] + bl3->x2[l]*bl3->x2[l] + bl3->x3[l]*bl3->x3[l]);
    if (ws->nosurf <= 0) {
      // modifs a.b.: r2 -> rmaxws pour la densite en w.s.
      rdeq = ws->rmaxws;
    } else {
      rdeq = ws->r0;
    }

    if ((rs - xr) <= rdeq) {
      if ((rs + xr) > rdeq) {
	vol = vol*0.5*(rdeq-rs+xr)/xr;
      }

      ia = bl3->ia1 + bl3->ia2;
      nl = 0;

      for(G4int i = 1; i <= ia; i++) {
	dx2 = std::pow((bl3->x1[l]-bl3->x1[i]),2) + std::pow((bl3->x2[l]-bl3->x2[i]),2) + std::pow((bl3->x3[l]-bl3->x3[i]),2);
	dp2 = std::pow((bl1->p1[l]-bl1->p1[i]),2) + std::pow((bl1->p2[l]-bl1->p2[i]),2) + std::pow((bl1->p3[l]-bl1->p3[i]),2);
      
	if((bl5->nesc[i] > 0) || (bl1->ind1[i] > 0) || (bl1->ind2[i] != bl1->ind2[l]) || (dx2 > xr2) || (dp2 > pr2)) {
	  if(((nl - 1)/vol/2.0) > 1.0) {
	    return 1.0;
	  } else {
	    return ((nl - 1)/vol/2.0);
	  }
	}
	nl = nl + 1;
      }
    } else {
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
  G4double plab = E*std::sqrt(E*E-3.52e6)/1876.6;
  G4double p1 = 0.001*plab;
  G4double alp = 0.0;

  if(plab > 2000.0) {
    // goto sel13;
    // sel13:
    return ((77.0/(p1 + 1.5))*scale);
  }

  if((m-1) < 0) {
    if (i == 0) {
      if (plab < 800.0) {
	if (plab < 450.0) {
	  alp = std::log(p1);
	  return(6.3555*std::exp(-3.2481*alp - 0.377*alp*alp));
	}
	else {
	  return((33.0 + 196.0*std::sqrt(std::fabs(std::pow((p1 - 0.95),5))))*scale);
	}
      }
      else {
	return(31.0/std::sqrt(p1)*scale);
      }
    }
  }

  if (plab < 800.0) {
    if (plab < 440.0) {
      return(34.0*std::pow((p1/0.4),(-2.104)));
    }
    else {
      return((23.5 + 1000.*std::pow((p1 - 0.7),4))*scale);
    }
  }
  else if(plab > 2000.0) {
    return ((77.0/(p1 + 1.5))*scale);
  }
  else {
    return((1250.0/(50.0 + p1) - 4.0*std::pow((p1 - 1.3),2))*scale);
  }
}

G4double G4Incl::totalCrossSection(G4double E, G4int m, G4int i)
{
  // total cross-sections
  // i=2,0,-2  for pp,pn,nn
  // m=0,1,2 for nucleon-nucleon,nucleon-delta,delta,delta

  G4double stotResult = 0.0;
  G4double sine = 0.0;

  if((m-1) < 0) { // Nucleon - Nucleon
    sine = deltaProductionCrossSection(E,int(i));
  }

  if((m-1) == 0) { // Nucleon - Delta
    sine = srec(E,(bl6->xx10),i,int(bl6->isa));
  }

  if((m-1) > 0) { // Delta - Delta
    sine = 0.0;
  }

  stotResult = sine + lowEnergy(E,m,i);
#ifdef G4INCLDEBUG
  theLogger->fillHistogram1D("totalCXResult", stotResult);
  theLogger->fillHistogram1D("totalLowECX", lowEnergy(E, m, i));
  if(m == 0) { // Nucleon - Nucleon
    if(i == 2) { // pp
      theLogger->fillHistogram2D("ppCX", E, stotResult);
    } else if(i == 0) { // pn
      theLogger->fillHistogram2D("pnCX", E, stotResult);
    } else if(i == -2) { // nn
      theLogger->fillHistogram2D("nnCX", E, stotResult);
    }
  } else if(m == 1) { // Nucleon - Delta
    theLogger->fillHistogram2D("NDCX", E, stotResult);
  } else if(m == 2) { // Delta - Delta
    theLogger->fillHistogram2D("DDCX", E, stotResult);
  }
#endif
  return stotResult;
}

G4double G4Incl::srec(G4double Ein, G4double d, G4int i, G4int isa)
{
  G4double E = Ein;
  G4double s = 0.0;
  G4double x = 0.0, y = 0.0;
  G4double srecResult = 0.0;

  if (i*i == 16) {
    return 0.0;
  }

  if(E <= (938.3 + d)) {
    return 0.0;
  } else {
    if(E < (938.3 + d + 2.0)) {
      E = 938.3 + d + 2.0;
    }
    s = E*E;
    x = (s - 3.523e6)/(s - std::pow((938.3 + d),2));
    y = s/(s - std::pow((d - 938.3),2));
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

//   G4double scali = 1.0;
//   G4double plab;
//   G4double p1;
//   G4double sproResult;

//   G4double EE = E - bl8->rathr;

//   if(EE*EE-3.53e6 < 0) {
//     return 0.0;
//   }
//   else {
//     plab = EE*std::sqrt(EE*EE-3.52e6)/1876.6;
//     p1 = 0.001*plab;
//     if (plab > 800.0) {
//       //goto spro1;
//       //  spro1:
//       if (i*i == 4) {
// 	//	goto spro10;
// 	//spro10:
// 	if (plab < 2000.0) {
// 	  //    goto spro11;
// 	  //  spro11:
// 	  if (plab < 1500.0) {
// 	    //      goto spro12;
// 	    // spro12:
// 	    sproResult = 23.5 + 24.6/(1.0 + std::exp(-10.0*p1 + 12.0)) - 1250.0/(p1+50.0)+4.0*std::pow((p1-1.3),2);
// 	    return (sproResult*scali);
// 	  }
// 	  else {
// 	    sproResult = 41.0 + 60.0*(p1 - 0.9)*std::exp(-1.2*p1) - 1250.0/(p1+50.0) + 4.*std::pow((p1 - 1.3),2);
// 	    return (sproResult*scali);
// 	  }
// 	}
// 	else {
// 	  return ((41.0 + (60.0*p1 - 54.0)*std::exp(-1.2*p1) - 77.0/(p1 + 1.5))*scali);
// 	}
//       }
//       else {
// 	if (plab < 2000.0) {
// 	  //goto spro2;
// 	  //  spro2:
// 	  if (plab < 1000.0) {
// 	    //goto spro3;
// 	    //  spro3:
// 	    return ((33.0 + 196.0*std::sqrt(std::pow(std::fabs(p1 - 0.95),5)) - 31.1/std::sqrt(p1))*scali);
// 	  }
// 	  else {
// 	    return ((24.2 + 8.9*p1 - 31.1/std::sqrt(p1))*scali);
// 	  }
// 	}
// 	else {
// 	  return ((42.-77./(p1+1.5))*scali);
// 	}
//       }
//     }
//     // plab <= 800.0
//     else {
//       return 0.0;
//     }
//   }
  double scali=1.0;
  double plab = 0.0;
  double p1 = 0.0;
  double sproResult = 0.0;

  double EE = E -(bl8->rathr);

  if(EE*EE-3.53e6 < 0) {
    goto spro22;
  }
  plab = EE * std::sqrt(EE*EE - 3.52e6)/1876.6;
  p1 = 0.001*plab;
  if (plab > 800.) {
    goto spro1;
  }
  spro22:
  sproResult = 0.0;
  return sproResult;

  spro1:
  if (i*i == 4) {
    goto spro10;
  }
  if (plab < 2000.) {
    goto spro2;
  }
  sproResult=(42.-77./(p1+1.5))*scali;
  return sproResult;
  spro2: if (plab < 1000.) {
    goto spro3;
  }
  sproResult=(24.2+8.9*p1-31.1/std::sqrt(p1))*scali;
  return sproResult;
  spro3: sproResult = (33.0 + 196.0*std::sqrt(std::pow(std::fabs(p1 - 0.95),5)) - 31.1/std::sqrt(p1))*scali;
  return sproResult;
  spro10: if (plab < 2000.) {
    goto spro11;
  }
  sproResult = (41.0 + (60.0*p1 - 54.0)*std::exp(-1.2*p1) - 77.0/(p1 + 1.5))*scali;
  return sproResult;
  spro11: if (plab < 1500.) {
    goto spro12;
  }
  sproResult = 41.0 + 60.0*(p1 - 0.9)*std::exp(-1.2*p1) - 1250.0/(p1 + 50.0)+4.0*std::pow((p1-1.3),2);
  sproResult = sproResult*scali;
  return sproResult;
 spro12: sproResult=23.5 + 24.6/(1.0 + std::exp(-10.0*p1 + 12.0)) - 1250.0/(p1 + 50.0) + 4.0*std::pow((p1 - 1.3),2);
  sproResult=sproResult*scali;
  return sproResult;
}

G4double G4Incl::pionNucleonCrossSection(G4double x)
{
  //   sigma(pi+ + p) in the (3,3) region
  //   new fit by j.vandermeulen + constant value above the (3,3)
  //   resonance

  G4double y = x*x;
  G4double q2 = (y-std::pow(1076.0,2))*(y-std::pow(800.0,2))/y/4.0;
  G4double q3 = 0.0, f3 = 0.0;
  G4double spn = 0.0;

  if(q2 <= 0) {
    return 0.0;
  } else {
    q3 = std::pow((std::sqrt(q2)),3);
    f3 = q3/(q3+std::pow(180.0,3));
    spn = 326.5/(std::pow(((x - 1215.0 - bl8->ramass)*2.0/110.0),2)+1.0);
    spn = spn*(1.0 - 5.0 * (bl8->ramass/1215.0));
    return (spn*f3);
  }
}

//      FUNCTION BARR(E,IZ,IA,IZN,R,V0)
//  tp =          transmissionProb(bl1->eps[bl9->l1]-fm,bl1->ind2[bl9->l1],1,itch,bl3->r2,v0);
G4double G4Incl::transmissionProb(G4double E, G4int iz, G4int ia, G4int izn,G4double R, G4double v0)
{
// C
// C
// C       NOUVEAU BARR RELATIVISTE ATTENTION AUX ENTREES SORTIES!
// C	(06/2005)
// C                                                                       P-N24450
// CCC   BARR=TRANSMISSION PROBABILITY FOR A PARTICLE (Nucleon, cluster or 
// CCC   pion) of kinetic energy E on the edge of the (attractive) well of 
// CCC   depth V0 felt by the particle.
// CCC   	IZ is the isospin of the particle,
// CCC   	IZN the instantaneous charge of the nucleus and R the radius of
// CCC       the well.
// CCC   	IA is the mass number of the particle, and MRE its mass energy.             
// C                                                                       P-N24500
// C     Modified 9/10/2002 for clusters (d,t,3He,4He) (IZ=isospin,IA=A)
// C     IZ must be correct so that charge Q=(IA+IZ)/2.
// C     Modified 4/2004 for relativistic expressions and pions.

  G4double mre;
  G4double barr = 0.0;
  G4double v0ia = v0*ia;
  mre = ia*938.3;

  if (E > v0ia && izn > 0 && R != 0) goto tp1;
  barr = 0.0;
  return barr;

 tp1:
  G4double x=std::sqrt((2.*mre*E+E*E)*(2.*mre*(E-v0ia)+std::pow((E-v0ia),2)));                                          
  barr=4.*x/(2.*mre*(2.*E-v0ia)+std::pow(E,2)+std::pow(E-v0ia,2)+2.*x);
  G4double iq=(ia+iz)/2.;
  if (iq > 0 && E > v0ia) goto tp2;
  return barr;

 tp2:
  G4double b=iq*izn*1.44/R;                                             
  G4double px=std::sqrt((E-v0ia)/b);                                             
  if (px < 1.) goto tp3;
  return barr;

 tp3:
  G4double g=iq*izn/137.03*std::sqrt(2.*mre/(E-v0ia)/(1.+(E-v0ia)/2./mre))
    *(std::acos(px)-px*std::sqrt(1.-px*px));
  if (g > 35.) {
    barr=0.0;
  } else {
    barr=barr*std::exp(-2.*g);
  } // endif
  return barr;
}

G4double G4Incl::ref(G4double &x1, G4double &x2, G4double &x3, G4double p1, G4double p2, G4double p3, G4double E, G4double r2)
{
  // Surface : modif de REF
  // REF=TIME NECESSARY FOR A NUCLEON TO REACH THE SURFACE

  const G4double pf = 270.33936, pf2 = 73083.4;
  G4double ref = 0.0;
  G4double t2 = p1*p1 +p2*p2 + p3*p3;
  G4double p = std::sqrt(t2);
  G4double r = r2;
  G4double xv = 0.0;
  G4double s_l = 0.0;
  G4double t1 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
   if (ws->nosurf <= 0) { // modif pour w.s.:
     xv = p/pf;
     if(t2 <= pf2) {
       r = interpolateFunction(xv);
#ifdef G4INCLDEBUG
     theLogger->fillHistogram1D("refInterpolXV", xv);
     theLogger->fillHistogram1D("refInterpolRes", r);
     if(xv > 1.0) theLogger->fillHistogram1D("refInterpolResWhenXVgt1", r);
#endif
     } else {
       r = ws->rmaxws;
     }
     r = r*r;
   }
/*
  if (ws->nosurf <= 0) {
    // modif pour w.s.:
    xv=p/pf;
    r=interpolateFunction(xv);
    r=r*r;
    if (t2 > pf2) r=std::pow(ws->rmaxws,2);
  } //endif
*/
 ref21:
  t4 = x1*x1 + x2*x2 + x3*x3;
  if (t4 > r) goto ref2;
  t1 = x1*p1 + x2*p2 + x3*p3;
  t3 = t1/t2;
  t5 = t3*t3 + (r-t4)/t2;
  if (t5 > 0) goto ref1;
  ref = 10000.0;
  return ref;
 ref1:
  ref = (-t3 + std::sqrt(t5))*E;
  return ref;
 ref2:
  s_l = std::sqrt(r*0.99/t4);
  x1 = x1*s_l;
  x2 = x2*s_l;
  x3 = x3*s_l;
  goto ref21;

  return 0.0;
}

// void G4Incl::forceAbsor(G4int nopart, G4double iarem, G4double izrem, G4double esrem, G4double erecrem,
// 			G4double alrem, G4double berem, G4double garem, G4double jrem)
void G4Incl::forceAbsor(G4int *nopart, G4int *iarem, G4int *izrem, G4double *esrem, G4double *erecrem,
			G4double *alrem, G4double *berem, G4double *garem, G4int *jrem)
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
  G4double iz1 = 0.0;
  G4double del = 0.0;
  G4double bmaxt = 0.0;
  G4double proba, proba_trans;
  G4double alea;

  if((*nopart) != -1) {
    return;
  }

  //  bl3->ia2 = int(std::floor(calincl->targetA() + 0.1)); // f(1) -> f[0]
  bl3->ia2 = calincl->targetA();
  sep = 6.8309;

  if(bl3->ia2 <= 4) {
    if(bl3->ia2 == 2) {
      itg = 6 - 1;
    }
    if(bl3->ia2 == 3 && calincl->targetZ() == 1) {
      itg = 7 - 1;
    }
    if(bl3->ia2 == 3 && calincl->targetZ() == 2) {
      itg = 8 - 1;
    }
    if(bl3->ia2 == 4) {
      itg = 9 - 1;
    }
    sep = light_gaus_nuc->vnuc[itg] - light_gaus_nuc->tfln[itg]; // :::BUG::: Off-by-one!!!
  }

  if((calincl->bulletE() >= 10.0) && (calincl->bulletE() <= 100.0)) {
    if(calincl->bulletType() == 1) {
      bl3->ia1 = int(1.0);
      iz1 = 1.0;
      G4double fmpinc = 938.2796;
      G4double pbeam2 = calincl->bulletE()*(calincl->bulletE() + 2.0*fmpinc);
      bmaxt = ws->bmax;
      proba_trans = coulombTransm(calincl->bulletE(),bl3->ia1,iz1,calincl->targetA(),calincl->targetZ());
      proba = forceAbs(1,calincl->targetA(),calincl->targetZ(),calincl->bulletE(),bmaxt,proba_trans);
      
      standardRandom(&alea,&(hazard->igraine[4]));
      if(alea > proba) {
        return;
      }

      (*iarem) = int(calincl->targetA()) + bl3->ia1;
      (*izrem) = int(calincl->targetZ()) + int(iz1);
      
      del = std::sqrt(std::pow(((calincl->targetA() + 1.0)*fmpinc + calincl->bulletE()),2) - pbeam2);
      
      (*erecrem) = pbeam2/((calincl->targetA() + 1.0)*fmpinc+calincl->bulletE() + del);

      (*esrem) = calincl->bulletE() + sep - (*erecrem);

      (*alrem) = 0.00001;
      (*berem) = 0.0;
      (*garem) = 0.99999;
      (*jrem) = 0;
      (*nopart) = 0;
      return;
    }
    else if((calincl->bulletType() == 2) && (calincl->bulletE() >= 20.0)) {
      bl3->ia1 = int(1.0);
      iz1 = 0.0;
      G4double fmpinc = 938.2796;
      G4double pbeam2 = calincl->bulletE()*(calincl->bulletE() + 2.0*fmpinc);
      bmaxt = ws->bmax;
      proba_trans = coulombTransm(calincl->bulletE(),bl3->ia1,iz1,calincl->targetA(),calincl->targetZ());
      proba = forceAbs(1,calincl->targetA(),calincl->targetZ(),calincl->bulletE(),bmaxt,proba_trans);
      
      standardRandom(&alea,&(hazard->igraine[4]));
      if(alea > proba) {
        return;
      }

      (*iarem) = int(calincl->targetA()) + bl3->ia1;
      (*izrem) = int(calincl->targetZ()) + int(iz1);
      
      del = std::sqrt(std::pow(((calincl->targetA()+1.)*fmpinc+calincl->bulletE()),2)-pbeam2);
      
      (*erecrem) = pbeam2/((calincl->targetA() + 1.0)*fmpinc + calincl->bulletE() + del);

      (*esrem) = calincl->bulletE() + sep - (*erecrem);

      (*alrem) = 0.00001;
      (*berem) = 0.0;
      (*garem) = 0.99999;
      (*jrem) = 0;
      (*nopart) = 0;
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
  double sig_g = 31.41592654*bmax*bmax;
  if(iprojo == 1) {
    ap = 1.0;
    zp = 1.0;
  }
  else {
    ap=1.0;
    zp=0.0;
  }

  sig_exp = xabs2(zp, ap, zt, at, ep);
  sig_incl = crossSection(int(iprojo), ep, at);
  
  proba = (sig_exp-pt*sig_incl)/(pt*(sig_g - sig_incl));
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
  
  G4double Const, xzt, xat, Const1 = 0.0, t1, gcm, bcm, plab, ecmp, ecmt, rela, ecm, rm, bigr, bigb;
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
  G4double rp = radius(G4int(ap));                                                     
  G4double rt = radius(G4int(at));                                                     
  G4double vp = (dp1 + dpth)*dppi*std::pow(rp,3);                                          
  G4double vt = (dp1 + dpth)*dppi*std::pow(rt,3);                                          
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

  gcm = (ap*(dp1 + E/938.0) + at)/(std::pow((std::pow(ap,2) + std::pow(at,2) + dp2*ap*(E + 938.0)*at/938.e0),dph));
  bcm = std::sqrt(dp1-dp1/(std::pow(gcm,2)));
  
  plab = ap*std::sqrt(dp2*938.0*E + E*E);                                    
  ecmp = gcm*(E+938.0)*ap - bcm*gcm*plab - ap*938.0;                     
  ecmt = gcm*938.0*at - at*938.0;
  rela = ecmp + ecmt;                                                    
  ecm = rela;                                                          
  if (ecm < (0.1*rela)) {
    ecm = 0.1*rela;
  }                             
  rm = (197.32/137.01)*zp*zt/ecm;                                      
  bigr = rp + rt + 1.2*(std::pow(ap,dpth) + std::pow(at,dpth))/(std::pow(ecm,dpth));                   

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
      xm = (1 - x1*std::exp(-E/(sl*x1)));                                     
      if (E < 20) {
	//	cout <<"e,xm= " << e << " " << xm << endl;
      }
    }                            
    else {                                                            
      xm = (1-0.3*std::exp(-(E-1)/15))*(1 - std::exp(-(E-0.9)));
    }
  }
  if ((nint(zp) == 2) && (nint(ap) == 4)) {
    Const = Const1 - 0.8/(1 + std::exp((250.0-E)/75.0));                      
  }
  if ((nint(zp) == 1) && (nint(ap) == 1)) {
    if (nint(at) > 45) {
      t1 = 40.0 + at/dp3;                                               
    }
    if (nint(at) < 4) {
      t1 = 55;                                                         
    }
    Const = 2.05 - 0.05/(1 + std::exp((250.0-E)/75.0));                     
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
      Const = Const - 0.25/(1 + std::exp((250.0 - E)/75.0));                  
    }
    if (nint(zt) >= 35) {
      phst = (nint(zt)-35.e0)/260.e0;	                                 
      Const = Const - phst/(1 + std::exp((250.0 - E)/75.0));                    
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
  
  ce = Const * (1.0 - std::exp(-E/t1)) - 0.292*std::exp(-E / 792) * std::cos(0.229*std::pow(E,0.453));    
  term1 = std::pow((at*ap),dpth)/(std::pow(at,dpth) + std::pow(ap,dpth));                           
  delta = 1.615*term1 - 0.873*ce;                                      
  delta = delta + 0.140*term1/(std::pow(ecm,dpth));                                 
  delta = delta + 0.794*(at - dp2*zt)*zp/(at*ap);                        
  delta = -delta;                   
  beta = 1.0;                                                           
  twxsec = dp10*dppi*1.26e0*1.26e0*beta*std::pow((0.873e0*std::pow(ap,dpth) + 0.873e0*std::pow(at,dpth)-delta),2);   
  sig = twxsec*(dp1-bigb/ecm)*xm;
  if (sig < dp0) {
    sig = dp0;
  }

  return sig;
}

void G4Incl::standardRandom(G4double *rndm, G4long*)
{
  // Use Geant4 G4UniformRand
  //  (*rndm) = G4UniformRand();
  (*rndm) = randomGenerator->getRandom();
}

void G4Incl::gaussianRandom(G4double *rndm)
{
  // Gaussian random number generator
  G4double xg, xxg, xshuf;
 gaussianRandom1:
  xg=0.;
  for(G4int i = 1; i <= 12; i++) { //do 6 j=1,12
    standardRandom(&xxg,&(hazard->ial));
    //  gaussianRandom6:
    xg=xg+xxg;
  }
  xg=xg-6.;
  if (xg*xg > 9.0) goto gaussianRandom1;
  standardRandom(&xshuf,&(hazard->igraine[10]));
  if (xshuf > 0.5) standardRandom(&xxg,&(hazard->ial));
  // gaussianRandom2:
  (*rndm) = xg;
  return;
    /*  
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
    */
}

G4double G4Incl::safeExp(G4double x)
{                                                
  if (x < -80.0) {
    x = -80.0;                                         
  }
  if (x > 80.0) {
    x = 80.0;
  }                                           

  return std::exp(x);
}

G4double G4Incl::radius(G4int A)
{                                               
  const G4double dp1 = 1.0, dp3 = 3.0; 
  const G4double dp5 = 5.0, dpth = dp1/dp3;
  const G4int naSize = 23;
  const G4int rmsSize = naSize;

  const G4int na[naSize] = {1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26};
  
  const G4double rms[rmsSize] = {0.85,2.095,1.976,1.671,2.57,2.41,2.519,2.45,
				 2.42,2.471,2.440,2.58,2.611,2.730,2.662,2.727,
				 2.900,3.040,2.969,2.94,3.075,3.11,3.06};
             
  G4double fact = std::sqrt(dp5/dp3);                                                
  //  G4int ia = int(std::floor(A+0.4));                                                        
  G4int ia = A;
  G4double result = fact * (0.84 * std::pow(A,dpth) + 0.55);                               
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

  G4double apow[3] = {0.0, 0.0, 0.0}, epow[5] = {0.0, 0.0, 0.0};
  G4int ii = 0, jj = 0;

  if(A >= 27.0) {
    ii = 3;        
    jj = 4;
  } else {
    ii = 3;
    jj = 5;
  }

  for(int i = 1; i <= ii; i++) {
    apow[i-1] = std::pow(A,(ii-i));
  }

  for(int j = 1; j <= jj; j++) {
    epow[j-1] = std::pow(E,(jj-j));
  }

  double result = 0.0;

  if(A >= 27.0) { 
    if(projectile == 1) {
      for(G4int i = 0; i < ii; i++) {
	for(G4int j = 0; j < jj; j++) {
	  result = result + coefp[j][i]*apow[i]*epow[j];
	}
      }
    } else {
      for(G4int i = 0; i < ii; i++ ) {
	for(G4int j = 0; j < jj; j++) {
	  result = result + coefn[j][i]*apow[i]*epow[j];
	}
      }
    }
  } else {
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
  G4double eta = 0.0, rho = 0.0;
  const G4double c2 = 0.00516;
  const G4double c3 = 0.007165;
  const G4double uma = 938.0;
  G4double ml = 0.0;
  
  G4double ecm = E*fm2/(fm1+fm2);
  G4double fm = fm1*fm2*uma/(fm1+fm2);
  G4double r1 = 0.0;
  G4double r2 = 0.0;
  if(fm1 >= 2.0) {
    r1 = 1.2*std::pow(fm1,0.33333333);
  }
  r2 = 1.2*std::pow(fm2,0.33333333);
  eta = c2*z1*z2*std::sqrt(fm/ecm);
  rho = c3*(r1+r2)*std::sqrt(fm*ecm);

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

  static G4int i0 = 0, j0 = 0;
  
  static G4double prob = 0.0;
  static G4double dumm = 0.0, x = 0.0, cx = 0.0;
  static G4double t = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
  static G4double f = 0.0, g = 0.0;
  static G4double temp = 0.0, temp1 = 0.0, temp2 = 0.0;
  static G4double xk = 0.0, delp0 = 0.0, delp1 = 0.0, delx0 = 0.0, delx1 = 0.0;

  if (rho > y) {                                                 
    if (psi > dp4 && psi < 50.0) {                           
      prob = clmb2(rho,eta,&dumm);                                       
    } else {
      if(eta <= 1.0e-6) { // Safeguard against a floating point exception
	x = 0.0;
      } else {
	x = std::exp(std::log(eta)/6.0);
      }
      prob = std::sqrt(dp1 - y*x/(c0 + c1 * std::pow(x,3) + rho * x));
    } 
    (*ml) = 0;
  } else {
    x = rho/y;                                                          
    if (psi <= psi0[0]) {
      t = min(pi*y,prm1);                                              
      cx = t/(std::exp(t) - dp1);                                             
      t1 = std::cos(rho) * (dp1-0.75*std::pow(psi,2) + dp5 * x * std::pow(psi,2)) - dph * psi * rho * std::sin(rho);
      t2 = dp1 + dph * psi * (dp1 - x/6.0);                                   
      if (eta > dp1) {                                          
	t3 = std::log(psi)+dp2*gamma - dp1 + dp1/(12.e0*std::pow(eta,2))+dp1/(12.e1*std::pow(eta,4));
      } else {                                                          
	t3 = std::log(dp2*rho) + gamma - dp1/(dp1 + std::pow(eta,2)) + s3*std::pow(eta,2) + s4*std::pow(eta,4);   
      }
      g = t1 + psi*t2*t3;                                                
      f = cx*rho*t2;                                                    
      prob = cx/(std::pow(g,2)+std::pow(f,2));                                           
      (*ml) = 3;
    } else if (psi <= psi0[ln0-1]) {
      if (x <= x0[0]) {
	temp = std::log(psi/psi0[0]);
	j0 = 1 + int(temp/delp0);                                        
	j0 = min(max(j0,1),(ln0-1));                                     
	temp = temp - (j0-1)*delp0;                                      
	t = f0[j0] + (f0[j0+1] - f0[j0])*temp/delp0;                       
	xk = x*std::sqrt(psi);                                              
	prob = (dp1+3.33e-1*x+3.e-1*xk+1.e-1*std::pow(xk,2))*std::exp(t);            
	t = min(pi*y,prm1);                                            
	cx = t/(std::exp(t)-dp1);                                           
	prob = cx/std::pow(prob,2);                                             
	(*ml) = 1;
      } else {                                                          
	temp1 = std::log(x/x0[0]);
	i0 = min(max(1 + int(temp1/delx0),1),lt0-1);                                     
	temp1 = temp1 - (i0 - 1)*delx0;                                    
	temp2 = std::log(psi/psi0[0]);
	j0 = min(max(1+int(temp2/delp0),1),ln0-1);                                     
	temp2 = temp2-(j0-1)*delp0;                                    
	t1 = trans0[i0][j0] + (trans0[i0+1][j0] - trans0[i0][j0]) * temp1/delx0;
	t2 = trans0[i0][j0+1] + (trans0[i0+1][j0+1] - trans0[i0][j0+1]) * temp1/delx0;
	prob = std::exp(t1 + (t2 - t1)*temp2/delp0);                                   
	(*ml)=2;                                                        
      }
    } else if (psi <= psi1[ln1-1]) {
      if (x <= x1[0]) {
	temp = std::log(psi/psi1[0]);
	j0 = min(max(1+int(temp/delp1), 1), ln1-1);
	t = f1[j0]+(f1[j0+1]-f1[j0])*(temp - (j0 - 1)*delp1)/delp1;                       
	xk = x*std::sqrt(psi);                                              
	prob = (dp1+3.33e-1*x+3.0e-1*xk+1.e-1*std::pow(xk,2))*std::exp(t);            
	t = min(pi*y,prm1);                                            
	cx = t/(std::exp(t)-dp1);                                           
	prob = cx/std::pow(prob,2);                                             
	(*ml) = 1;
      } else {                                                         
	temp1 = std::log(x/x1[0]);
	i0 = min(max(1+int(temp1/delx1),1),lt1-1);                                    
	temp1 = temp1-(i0-1)*delx1;                                    
	temp2 = std::log(psi/psi1[0]);
	j0 = min(max(1+int(temp2/delp1),1),ln1-1);                                     
	temp2 = temp2 - (j0-1)*delp1;                                    
	t1 = trans1[i0][j0] + (trans1[i0+1][j0] - trans1[i0][j0])*temp1/delx1;
	t2 = trans1[i0][j0+1] + (trans1[i0+1][j0+1] - trans1[i0][j0+1])*temp1/delx1;                                                      
	prob = std::exp(t1 + (t2-t1)*temp2/delp1);                                    
	(*ml)=2;                                                        
      }
    } else {
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
  static G4double x = 0.0, temp = 0.0, prob = 0.0;
  static G4int i = 0;

  x = dp1/(dp1 + std::sqrt(dph*rho*eta));                                     
  if (x < x1) {                                                 
    temp = t0[2] * std::pow((x/x1),dpth);
  } else {
    i = int(std::floor(xi*x));                                                          
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

// C****************************************************************************

void G4Incl::projo_spec(G4int, G4int ips,
		   G4double, G4double, G4double) // max_a_proj is a global #define'd as a global constant
{
  //      subroutine projo_spec(ia1,ips,max_a_proj,
  //     sn_projspec,fmpinc,pinc,tlab)
  //      common/bl1/p1(300),p2(300),p3(300),eps(300),ind1(300),ind2(300),ta
  // spectators of the projectile treated here (nucleus for fermi-breakup):
  // ips is the number of spect, their label is in n_projspec(ips).
  // e,p are outside the target potential for these nucleons.
  // la transmission de max_a_proj=60 ne marche pas....explicite ici!
  //      integer a_projspec,z_projspec,n_projspec(60)
  //      real*4 m_projspec
  //      common/proj_spect/a_projspec,z_projspec,ex_projspec,t_projspec,
  //     s p1_projspec,p2_projspec,p3_projspec,m_projspec
  ///      common/quatvectprojo/eps_c(60),p3_c(60),
  //     s  p1_s(60),p2_s(60),p3_s(60),
  //     s  t_c(60)
  // (eps_c,p1_s,p2_s,p3_s,eps_c used to store the kinematics of nucleons
  // for composit projectiles before entering the potential)
  //      common/fermi/tf,pf,pf2     

  //  G4double ttrou[21];
  G4double fja = 0.0;
  G4double fjz = 0.0;
  //  G4double ex2;
  //  G4double rmmass;
  //  const G4double um = 931.20793;

// c verifs: ground state
// c      es=0.
// c      p1s=0.
// c      p2s=0.
// c      p3s=0.
// c      gsmass=0.
// c      do iloc=1,ia1
// c      es=es+eps_c(iloc)
// c      p1s=p1s+p1_s(iloc)
// c      p2s=p2s+p2_s(iloc)
// c      p3s=p3s+p3_s(iloc)
// c      gsmass=gsmass+sqrt(eps_c(iloc)**2
// c     s        -p1_s(iloc)**2-p2_s(iloc)**2-p3_s(iloc)**2)
// c      enddo
// c      sgs=sqrt(es**2-p1s**2-p2s**2-p3s**2)
// c      write(6,*) 'sqs, gsmass :',sgs,gsmass
// c      write(6,*) 'es,p1s,p2s,p3s:',es,p1s,p2s,p3s
// c end verifs
      
//  G4cout <<"proj_spect called with: " << G4endl;
//  G4cout <<"ips = " << ips << G4endl;
  ps->a_projspec=0;
  ps->z_projspec=0;
  ps->t_projspec=0.0;
  G4double e_spec=0.;
    G4double p1_spec=0.;
    G4double p2_spec=0.;
    G4double p3_spec=0.;
    G4double gs_mass=0.;

  if(ips > 1) {
    for(G4int i_spec=1; i_spec <= ips; i_spec++) {
      G4int i_c = ps->n_projspec[i_spec];
      e_spec = e_spec+bl1->eps[i_c];
      p1_spec = p1_spec+bl1->p1[i_c];
      p2_spec = p2_spec+bl1->p2[i_c];
      p3_spec = p3_spec+bl1->p3[i_c];
      G4double pi_spec2= std::pow(bl1->p1[i_c],2) + std::pow(bl1->p2[i_c],2) + std::pow(bl1->p3[i_c],2);
      // current ground state mass:
      gs_mass=gs_mass+std::sqrt(std::pow(bl1->eps[i_c],2)-pi_spec2);
      if(bl1->ind2[i_c] == 1) ps->z_projspec=ps->z_projspec+1;
    } // enddo
    G4double p_spec2=std::pow(p1_spec,2)+std::pow(p2_spec,2)+std::pow(p3_spec,2);
    G4double s_spec = std::sqrt(std::pow(e_spec,2)-p_spec2);

    // no projectile spectator if a>=4 and a=z or a=n (no dresner breakup)
    if(ips >= 4 && (ps->z_projspec == ips || ps->z_projspec == 0)) {
      for(G4int i_spec = 1; i_spec <= ips; i_spec++) { //do i_spec=1,ips
	G4int i_c = ps->n_projspec[i_spec];
	// single nucleon forced on shell! momentum kept.
	// single nucleon forced on shell! energy kept.
	G4double arg=(std::pow(bl1->eps[i_c],2)-std::pow(938.2796,2))/
	  (std::pow(bl1->p1[i_c],2)+std::pow(bl1->p2[i_c], 2)+std::pow(bl1->p3[i_c],2));
	G4double coef=0.0;
	if(arg < 0.) { // Single nucleon is forced on shell. We have
		       // to choose between energy and momentum
		       // conservation... so we choose energy:
	  coef = 1.0;
	  bl1->eps[i_c]=std::sqrt(std::pow(bl1->p1[i_c],2) + std::pow(bl1->p2[i_c], 2) + 
				  std::pow(bl1->p3[i_c], 2) + std::pow(938.2796,2));
	} else {
	  coef=std::sqrt(arg);
	  bl1->p1[i_c]=bl1->p1[i_c]*coef;
	  bl1->p2[i_c]=bl1->p2[i_c]*coef;
	  bl1->p3[i_c]=bl1->p3[i_c]*coef;
	}
      } // enddo
      return;
    } // endif
	
    // excitation energy from particle-hole excitation:
    //   seconde methode de calcul ex4:
    for(G4int i = 1; i <= bl3->ia1; i++) { // do i=1,ia1
      G4int nb_p=i;
      //      call ordered(t_c(i),nb_p,tpart)
      ordered(qvp->t_c[i], nb_p);
    } // enddo
    G4double ex4=0.;
    for(G4int i = 1; i <= ips; i++) { //do i=1,ips
      G4int j=ps->n_projspec[i];
      ex4=ex4+qvp->t_c[j]-ps->tab[i]; //ex4=ex4+qvp->t_c[j]-ps->tpart[i];
    } // enddo

    ps->ex_projspec =ex4;

    // coherent rest mass:
    ps->m_projspec = s_spec - ps->ex_projspec;

    ps->a_projspec=ips;

    if((ps->ex_projspec < 0.) || (ps->m_projspec < 0.)) {
      G4cout <<"ex_projspec is negative,a,z: " << ps->ex_projspec << " " << fja << " " << fjz << G4endl;
      G4cout <<"...or m is negative,a,z:" << ps->m_projspec << " " << fja << " " << fjz << G4endl;
    } // endif

    ps->t_projspec = e_spec-s_spec;
    ps->p1_projspec = p1_spec;                  
    ps->p2_projspec = p2_spec;
    ps->p3_projspec = p3_spec;
      
// c check module...      
// c      write(6,*) 'projectile spectator a,z:',a_projspec,z_projspec
// c      write(6,*) 'm,t,e*:',m_projspec,t_projspec,ex_projspec
// c      write(6,*) 'p1,p2,p3:',p1_projspec,p2_projspec,p3_projspec
// c      write(6,*) '   constituants:',ips
// c      do i_spec=1,ips
// c	        i_c = n_projspec(i_spec)
// c      write(6,*) i_c,eps(i_c),p1(i_c),p2(i_c),p3(i_c)       
// c      enddo
// c ...end check module 

    return;
  } else if(ips == 1) { // then
    ps->a_projspec = 0;	// number of nucleons in the spectator nucleus
    G4int i_c = ps->n_projspec[1];
    //c single nucleon forced on shell! momentum kept.
    //c      eps(i_c)=sqrt(p1(i_c)**2+p2(i_c)**2+p3(i_c)**2+938.2796**2)     
    //c single nucleon forced on shell! energy kept.
    G4double arg=(std::pow(bl1->eps[i_c],2)-std::pow(938.2796,2))/
      (std::pow(bl1->p1[i_c],2)+std::pow(bl1->p2[i_c],2)+std::pow(bl1->p3[i_c],2));
    G4double coef=0.0; // Initialize the coefficieant
    if(arg < 0.) {
      G4cout <<"off shell problems in projo_spec, forced coef = 1.0 to recover from the error" << G4endl;
      coef = 1.0;
    } else {
      coef = std::sqrt(arg);
    }

    bl1->p1[i_c]=bl1->p1[i_c]*coef;
    bl1->p2[i_c]=bl1->p2[i_c]*coef;
    bl1->p3[i_c]=bl1->p3[i_c]*coef;
  } //endif

  ps->a_projspec = 0; // number of nucleons in the spectator nucleus
  return;
}
      
void G4Incl::ordered(G4double t, G4int nb)
{
  //      SUBROUTINE ORDERED(T,nb,TAB)
  //C This subroutine will ordered in TAB() from min TAB(1) to max
  //C TAB(nb), nb values of T sent successively. 
  //C CALL ORDERED (T1,1), CALL ORDERED (T2,2),..CALL ORDERED (Tnb,nb) 
  //      DIMENSION TAB(1)

  //  G4cout <<"Ordered was called" << G4endl;
  if(nb == 1) {
    ps->tab[1]=t;
    return;
  } // endif

  G4int i=0;

 ordered100:
  i=i+1;
  if(t < ps->tab[i]) {
    //    for(G4int l = nb-1; l != i; l--) { //do l=nb-1,i,-1
    for(G4int l = nb-1; l >= i; l--) { //do l=nb-1,i,-1
      //      G4cout <<"Ordered: L = " << l << " nb = " << nb <<  " i = " << i << G4endl;
      ps->tab[l+1]=ps->tab[l];
    } // enddo
    ps->tab[i]=t;
    return;
  } else {
    if(i <= nb-1) goto ordered100;
    ps->tab[i]=t;
    return;
  } // endif
  return;
}

// Utilities

G4double G4Incl::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  } else {
    return b;
  }
}

G4int G4Incl::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  } else {
    return b;
  }
}

G4double G4Incl::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  } else {
    return b;
  }
}

G4int G4Incl::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  } else {
    return b;
  }
}

G4int G4Incl::nint(G4double number)
{
  G4double intpart;
  G4double fractpart;
  fractpart = std::modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return int(std::floor(number));
    } else {
      return int(std::ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return int(std::floor(number));
    } else {
      return int(std::ceil(number));
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
  return std::sqrt(d*d-a*a-b*b-c*c);
}

G4double G4Incl::pcm(G4double E, G4double A, G4double C)
{
  return (0.5*std::sqrt((std::pow(E,2)-std::pow((A+C),2))*(std::pow(E,2)-std::pow((A-C),2)))/E);
}

G4double G4Incl::sign(G4double a, G4double b)
{
  if(b >= 0) {
    return utilabs(a);
  }
  if(b < 0) {
    return (-1.0*utilabs(a));
  }

  if(verboseLevel > 2) {
    G4cout <<"Error: sign function failed. " << G4endl;
  }
  return a; // The algorithm is never supposed to reach this point.
}

G4double G4Incl::utilabs(G4double a)
{
  if(a > 0) {
    return a;
  }
  if(a < 0) {
    return (-1.0*a);
  }
  if(a == 0) {
    return a;
  }

  if(verboseLevel > 2) {
    G4cout <<"Error: utilabs function failed. " << G4endl;
  }
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
  return (std::sqrt(a*a+b*b+c*c+d*d));
}

G4int G4Incl::idnint(G4double a)
{
  G4int value = 0;

  G4int valueCeil = int(std::ceil(a));
  G4int valueFloor = int(std::floor(a));

  if(std::abs(value - valueCeil) <= std::abs(value - valueFloor)) {
    return valueCeil;
  }
  else {
    return valueFloor;
  }
}

// Initialization helper
void G4Incl::clearState() {
  G4int epSize = 300;
  for(G4int i = 0; i < epSize; ++i) {
    kind[i] = 0;
    ep[i] = 0.0;
    alpha[i] = 0.0;
    beta[i] = 0.0;
    gam[i] = 0.0;
  }
  inside_step = 0;
}

// C------------------------------------------------------------------------
// C Logging
// C

// C     Start a new simulation step log entry
//       subroutine print_log_start_step()
void G4Incl::print_log_start_step()
{
  if(inside_step != 1) {
    G4cout << "(simulation-step " << G4endl;
    inside_step = 1;
  } else {
    G4cout << ";; Error: Can not start a new step inside s step" << G4endl;
  }
}

// C     End a simulation step
//       subroutine print_log_end_step()
void G4Incl::print_log_end_step()
{
  if(inside_step == 1) {
    G4cout <<") ;;ENDSIMULATIONSTEP" << G4endl;
    inside_step = 0;
  } else {
    G4cout <<";; Error: Not inside a step" << G4endl;
  }
}

// C     The main logging routine
// C     Parameters are used to select what to print. Possible options are:
// C     list of avatars, selected avatar, list of particles and
// C     the index of the selected avatar
//       subroutine print_log_entry(iavatars, iselected, iparticles, imin)
void G4Incl::print_log_entry(G4int iavatars, G4int iselected, G4int iparticles, G4int imin)
{
//       COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TA
//       COMMON/BL2/CROIS(19900),K,IND(20000),JND(20000) 
//       COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22270

  if(iavatars == 1) {
    G4cout <<";; List of avatars" << G4endl;
    print_avatars();
  }

  if(iselected == 1) {
    G4cout <<";; Selected avatar" << G4endl;
    print_one_avatar(imin); // (avatar ...) ;; Selected avatar (index: imin)
  }
      
  if(iparticles == 1) {
    G4cout <<";; Particle avatar map" << G4endl;
    print_map();      // ((list particle-avatar-map particle-avatar-map ...))
  }
}

//       subroutine print_avatars(crois, k, ind, jnd, p1, p2, p3, eps, ind1
//      s,ind2, ta, r1, r2, x1, x2, x3, ia1, ia2, rab2)
void G4Incl::print_avatars()
{
// C      COMMON/BL1/P1(350),P2(350),P3(350),EPS(350),IND1(350),IND2(350),TA
// C     -,V(350)                                                           P-N22260
// C      COMMON/BL2/CROIS(19900),K,IND(20000),JND(20000) 
// C      COMMON/BL3/R1,R2,X1(350),X2(350),X3(350),IA1,IA2,RAB2             P-N22270
      
  G4cout <<"(list ;; List of avatars" << G4endl;
  for(int i = 1; i <= bl2->k; ++i) { //do i=1,k
    print_one_avatar(i);
  }
  G4cout << ")" << G4endl;
}

// C     Prints the information of one avatar (selected using the index
// C     (values between 1 and K) in BL2 common block.
//       subroutine print_one_avatar(index)
void G4Incl::print_one_avatar(G4int index)
{
//       COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TA
//       COMMON/BL2/CROIS(19900),K,IND(20000),JND(20000) 
//       COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22270

  if(index < 0) {
    G4cout <<";; Error! index < 0! index = " << index << G4endl;
    return;
  }

  G4int i_ind1 = bl1->ind1[bl2->ind[index]];
  //  G4int i_ind2 = bl1->ind2[bl2->ind[index]];
  G4int j_ind1 = -1;
  G4int j_ind2 = -1;

  if(bl2->jnd[index] != -1) {
    j_ind1 = bl1->ind1[bl2->jnd[index]];
    j_ind2 = bl1->ind2[bl2->jnd[index]];
  } else {
    j_ind1 = -1;
    j_ind2 = -1;
  }

  G4cout <<"(avatar " << bl2->crois[index] << G4endl;// Create avatar with time crois(index)

// C     Now we print the avatar type
  if(i_ind1 == 0 && j_ind2 == 0) { // NN collision
    G4cout <<"(quote nn-collision)" << G4endl;
  } else if(bl2->ind[index] != -1 && bl2->jnd[index] == -1) {
    G4cout <<"(quote reflection)" << G4endl;
  } else {
    G4cout <<"(quote unidentified-avatar)" << G4endl;
  }

// C     Print the list of particles in the avatar
  G4cout <<"(list ;; List of particles in the avatar" << G4endl;
  print_one_particle(bl2->ind[index]);
  if(bl2->jnd[index] != -1) {
    print_one_particle(bl2->jnd[index]);
  }
  G4cout <<") ;; End of the particle list" << G4endl;

  G4cout <<") ;; End of the avatar" << G4endl;
}

// C     Prints the information of one particle in common blocks BL1 and BL3.
// C     The particle is selected by giving the index of the particle entry
// C     as a parameter.
//       subroutine print_one_particle(index)
void G4Incl::print_one_particle(G4int index)
{
//       COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TA
//       COMMON/BL2/CROIS(19900),K,IND(20000),JND(20000) 
//       COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22270

// C     Particle knows its ID, type position, momentum and energy
  G4cout <<"(particle " << index << G4endl;
  if(bl1->ind1[index] == 0) { // Nucleon
    if(bl1->ind2[index] == 1) {
      G4cout << "(quote proton)" << G4endl;
    } else if(bl1->ind2[index] == -1) {
      G4cout << "(quote neutron)" << G4endl;
    } else {
      G4cout << "(quote unidentified-nucleon-ind1-" << bl1->ind1[index] <<"-ind2-" << bl1->ind2[index] << G4endl;
    }
  } else if(bl1->ind1[index] == 1) { // Delta
    if(bl1->ind2[index] >= -2 && bl1->ind2[index] <= 2) {
      G4cout << "(quote delta-"
	     << bl1->ind2[index]
	     << ")" << G4endl;
    } else {
      G4cout << "(quote unidentified-delta)" << G4endl;
    }
  } else {
    G4cout <<"(quote unidentified-particle)" << G4endl;
  }

  print_three_vector (bl3->x1[index], bl3->x2[index], bl3->x3[index]);
  print_three_vector (bl1->p1[index], bl1->p2[index], bl1->p3[index]);

  G4cout << bl1->eps[index] << G4endl;

  G4cout <<")" << G4endl; // Close the particle s-expression
}

//       subroutine print_three_vector(x, y, z)
void G4Incl::print_three_vector(G4double x, G4double y, G4double z)
{
  G4cout <<"(vector3 " << x << " " << y << " " << z << ")" << G4endl;
}

//       subroutine print_map()
void G4Incl::print_map()
{
//       COMMON/BL1/P1(300),P2(300),P3(300),EPS(300),IND1(300),IND2(300),TA
//       COMMON/BL2/CROIS(19900),K,IND(20000),JND(20000) 
//       COMMON/BL3/R1,R2,X1(300),X2(300),X3(300),IA1,IA2,RAB2             P-N22270

  G4cout <<"(list ;; Particles at the beginning of the avatar" << G4endl;
  G4int ia = bl3->ia1 + bl3->ia2;
  for(G4int i = 1; i <= ia; ++i) { // do i=1,IA
    print_one_particle(i);
  }
  G4cout <<")" << G4endl;
}
