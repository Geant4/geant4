#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
//#include <iostream.h>
//#include <fstream>
//#include "../include/G4Cascade.hh" 

//#include "G4NuclearFermiDensity.hh"
//#include "G4NuclearShellModelDensity.hh"
//#include "G4NucleiPropertiesTable.hh"
//#include "Randomize.hh"

#include "G4ios.hh"
//#include "g4rw/tvvector.h"

#include "G4BertiniCascade.hh"

G4BertiniCascade::G4BertiniCascade() {
ifirst = 0 ;
}

G4BertiniCascade::~G4BertiniCascade(){
  ;
}

G4VParticleChange* G4BertiniCascade::ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus) {
  

  G4VParticleChange *particles = new G4VParticleChange();
  G4cout << "particle update: " << particles << G4endl;

  return particles;
 
}

G4ReactionProductVector* G4BertiniCascade::Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{

  return new G4ReactionProductVector;

}

G4int G4BertiniCascade::operator == (G4BertiniCascade& right)
{

  return (this == &right);

}

G4int G4BertiniCascade::operator != (G4BertiniCascade& right)
{

  return (this != &right);

}



void G4BertiniCascade::test(){

  G4cout << "Testing G4 framework" << G4endl;

  // particles
  G4ParticleDefinition *aPiMinus = G4PionMinus::PionMinus();
  G4double m = aPiMinus->GetPDGMass()/MeV;
  G4double a(aPiMinus->GetPDGCharge()); //GetBaryonNumber, Get4Momentum
  G4cout << "pi- mass" << m << "charge " << a << G4endl;


  // dynamic particle 
  G4DynamicParticle *dp = new G4DynamicParticle();
  dp->SetDefinition( aPiMinus );
  dp->SetKineticEnergy( 3.0 - m );
  dp->SetMomentum( 3.0 );
  dp->SetMomentumDirection( 0.0, 0.0, 0.0);  
  G4cout << dp->GetDefinition() << ": " << "eKin: " << dp->GetKineticEnergy() << G4endl;
  G4cout << "momentum p: " << dp->GetMomentum() << " p direction: " << dp->GetMomentumDirection() << G4endl;

  // vectors
  G4ThreeVector v;
  v.setX(0.0); v.setY(0.0); v.setZ(0.0);
  G4cout << " v = " << v << G4endl;

#define G4Vector G4std::vector
  typedef G4Vector<G4DynamicParticle* > DPvector;

  typedef G4Vector<G4double> Dvector;
  Dvector w;
  w.reserve(6);  
  w.erase( w.begin(), w.end() );

  w.push_back( 1.0 );
  w.push_back( w[0] );
  G4cout << " vector w: " << w[0] << G4endl;
 

  // nucleus
  G4Nucleus *myNuc;
  myNuc->AddExcitationEnergy(500.0 * MeV);
  G4cout << "nucleus exitation e: " << myNuc->GetEnergyDeposit() / GeV << " GeV" << G4endl;


  /*

    theCascade->test();  
    *theCascade.Propagate(G4KineticTrackVector* theSecondaries, 
    G4V3DNucleus* theNucleus);
    G4ParticleDefinition* theProtonDefinition = G4Proton::ProtonDefinition();
    G4double              theProtonFormationTime = 0.0*ns;
    G4ThreeVector         theProtonInitialPosition 
    (0.0*cm, 0.0*cm, 0.0*cm);
    G4LorentzVector       theProton4Momentum 
    (0.0*MeV, 0.0*MeV, 1.0*MeV, 938.27*MeV);  
    G4KineticTrack        aProton(theProtonDefinition, 
    theProtonFormationTime,
    theProtonInitialPosition,
    theProton4Momentum);
                          
    G4Fancy3DNucleus      theCarbon;
    G4double              theA = 12.0;
    G4double              theZ =  6.0;
    theCarbon.Init(theA, theZ);
   
    G4cout << G4endl << " Nucleus mass number = " << theCarbon.GetMassNumber() << G4endl
    << " Nucleus mass = " << theCarbon.GetMass() << G4endl
    << " Nucleus charge = " << theCarbon.GetCharge() << G4endl;
             
    G4KineticTrackVector theProjectileParticles;
    theProjectileParticles.insert(new G4KineticTrack(aProton));
    //   theProjectileParticles.insert(new G4KineticTrack(aProton));

  */

}

void G4BertiniCascade::main() {
  i1 = 0;
  ncasca = 0;
  no5rca = 0;
  nevaph = 0;
  nevapl = 0;
  nerupl = 0;
  nerup  = 0;
  // inpcem(); //:::
  // input2(); //:::
  tmel1  = 0.0;
  jtimoh = 0;
  itim   = iclock();
  jtimo  = 0;
  xtimo  = 0.0;
  noevt  = 0;
  nomaxh = 0;
  nomax  = 0;
  ibert  = 0; // means hetc configuration tape is read at first call to bert
  // note that if imash. = 1 subroutine mashnk in sorcem.f may set icem = 0 if (some nuclei so that sub bert is called without being
  // the first nucleus collided with;hence need to set ibert = 0 if first appearance of these to read hetc tape. this is done in cascad at stm. 430.
  iqev   = 0; // zero some counters used in cascad
  ibrt   = 0;
  iqheh  = 0;
  ihec   = 0;
  iskl   = 0;
  ipcl   = 0;
  mmf    = 0;
  main1();
} 
void G4BertiniCascade::main1() {
  // input(); //:::
  // input2(); //:::
  nil    = 0;
  nomaxn = nomax;
  nspltt = 0;
  iqhist = 0;
  if (neutp > 0 && nbertp > 0) {
    // rewind neutp; //:::
  }
  for (G4int i = 0; i < 25000; ++i) gmstor[i] = 0.0;
  G4int ngeom = 69497;
  gmstor[1] = ngeom;
  // jomin(ngeom, in, io); //:::
  //  if (verboseLevel > 1) G4cout << " geom requires " << ngeom << "locations " << "i = " << i <<endl;
  if (ngeom > 25000) G4Exception("main1");   
  mainBodyInit();
} 

void G4BertiniCascade::mainBodyInit() {
  ncol = -1;
  for (G4int i = 0; i < 16; ++i) {
    pevts[i]  = 0.0;
    revts[i]  = 0.0;
    phevts[i] = 0.0;
    rhevts[i] = 0.0;
    pdevts[i] = 0.0;
    peevts[i] = 0.0;
    reevts[i] = 0.0;
    hpevt[i]  = 0.0;
    hrevt[i]  = 0.0;
    hpevth[i] = 0.0;
    hrevth[i] = 0.0;
    hpdevt[i] = 0.0;
    peiv[i]   = 0.0;
  }
  nkey = 1;
  ncol = -1;
  analz1(ncol); 
  sors(ncol);
  nocas = 1;
  neutno = 0;
  wtav = 0.0;
  cneuav = 0.0;
  edtotn = 0.0;
  nhcasc = 0;
  for (i = 0; i < 9; ++i) ncoutp[i] = 0;
  nobch = 1;
  negex = 0;
  lowaz = 0;
  nwtc = 0;
  ncas = 0;
  nsav = 0;
  nfirst = 1;
  nbelo = 0;
  nabov = 0;
  mainBody();
} 

void G4BertiniCascade::mainBody() {
  //
  //
  no = 1;
  // the new history starts here
  // also a new batch;  some zeroing if (batch is done after end
  // of batch message on the history tape and after call to analz1
  noo = no;
  ncc = maxcas * maxbch;
  if (nocas == 1) {
    maxcso = maxcas;
    maxbco = maxbch;
    ncco = ncc;
  }
  if (maxcas != maxcso || maxbch != maxbco || ncc != ncco) {
    if (verboseLevel > 1) {
      G4cout << "num = " << num << " nhist = " << nhist 
             << " noevt = " << noevt << "maxcas = " << maxcas 
             << " maxbch = " << maxbch << " ncc = " << ncc << " nocas = " 
             << nocas << " nobch = " << nobch << " nodat = " << nodat 
             << "main" <<endl;
    }
  }
  niil    = 20;
  nomaxn  = nomax;
  nil     = 10;
  nomaxn  = nomax;
  nil     = 8;
  name[1] = 1;
  namax   = 1;
  nomax   = 1;
  ncol    = 1;
  sors(ncol);
  nhis    = (nobch - 1) * maxcas + nocas;
  nameol  = 0;
  if (nfirst == 0) {
    nfirst = 0;
    nhkey = 0;
    for (G4int i = 0; i < 17; ++i) {
      if (denh[i] > 0.0) {
	nhkey = 1;
	break;
      }
    }
    tipstr = tip[no];
    ncmc  = 0;
    if (tipstr <= 1.0) {
      if (emax > ehin) {
	ncmosc = 1;
	q = 0.0;
	iq = 0;
      }
    } if (tipstr <= 4.0 && emax > ehipi) {
      ncmosc = 1;
      q  = 0.0;
      iq = 0;
    }
    oldwt = wt[no];
    wtcut = 0.0 * wt[no];
    if (energy[1] > emax) {
      if (verboseLevel > 1) G4cout << " the source energy e[1] = " << energy[1] << " computed in sors exceeds emax" << G4endl;
      G4Exception("main2");
    }
    if (ncol == -4) mainFinal();
    //gomsor(x, y, z, nmed, blz); //::: defined in geom.f 
    mat  = nmed[1];
    ncol = 1;
    analz1(1);
    noo = no;
    mainBody2();
  }
} 

void G4BertiniCascade::mainBody2() {
  ityp = static_cast<G4int>(tip[no]) + 1; 
  ++noevt;
  ic = ityp;
  if (ic <= 30) ic = ijevnt[ic];
  noo = no;
  kstor = 0;
  nofask = 0; // nofask is the number of particles not transported and put in datalo
  if (tip[no] == 5.0 || tip[no] == 6.0) kindo[no] = static_cast<G4int>(tip[no]) + 1; 
  if (tip[no] == 5.0 || tip[no] == 6.0) barr[no] = 0; // above is to account for muons
  if (nomaxn > nomaxh) nomaxh = nomaxn;
  // note: pi0 were put in comon/ but are knocked out at 102 if they come directly from sors.
  // note: other cases were given e(no)=-1 at 160 and come to 101 to end
  // note: particle with wt(no).le.wtlow are sent to 160 and then to 101 to be killed with e(no)=-1
  ic = kindo[no];
  if (ic <= 30) ic = ijevnt[ic];
  ic = kindo[no];
  if (ic <= 30) ic = ijevnt[ic];
  nhcas = (nobch - 1) * maxcas + nocas;
  if (nhcas != nhist) {
    ++nhcasc;
    nhistt = nhist;
  } if (energy[no] < 1.0) {
    if (main3()) mainBody2();
    mainFinal();
  } else if (energy[no] == 1.0) {
    if (ityp != 5) {
      if (main3()) mainBody2();
      mainFinal();
    }
    mat = nmed[no];
    bold = blz[no];
    if (npidk > 0) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    }
    ec[no] = energy[no];
    xc[no] = x[no];
    yc[no] = y[no];
    zc[no] = z[no];
    oldwt = wt[no];
    getrig();
    G4int test;
    test = main2();
    while (test == 0)
      if (test == 1) mainBody2();
    if (test == 2) {
      mmf = 5;
      m = 5;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      ncol = 5;
      nparto = nopart;
      if (ibertp < 0) mfpd2(5); 
      if (nseudo > 0) {
	ncol = 5;
	analz1(5); 
	nabov = 0;
	nbelo = 0;
      } else {
	nabov = 0;
	nbelo = 0;
	if (nsav <= 0) {
	  nsav = 1;
	  es = energy[no];
	  xs = x[no];
	  ys = y[no];
	  zs = z[no];
	}
      }
      if (wt[no] <= wtlow && wt[no] > 0.0) {
	energy[no] = -1.0;
	ic = kindo[no];
	if (ic <= 30) ic = ijevnt[ic];
	mainBody2();
      }
      if (wt[no] == 0.0) {
	energy[no] = -1.0;
	ic = kindo[no];
	if (ic <= 30) ic = ijevnt[ic];
	mainBody2();
      }
      mainBody4();
    }
    if (test == 3) mainBody8();
    if (test == 4) mainBody9();
  } else if (energy[no] > 1.0) {
    if (tip[no] == 3) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    }
    main5();
  }
  mainBody3();
}

void G4BertiniCascade::mainBody3a() {      
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (nspred > 0 && name[no] == 1 && ityp != 2 && mat != 6666) {
    mark = 1;
    dist = d;
    //sprd(ind, rng, rbtors, mark, dist);//:::
    // cascad(collisions), depending upon nbeta
    // NOTE: if ind = 4 mainBody3c() and then either to decay (m = 3) or to
    // getrig and cascad, depending upon nbeta
    switch(ind) {
    case 1:
      if (verboseLevel > 1) {
	G4cout << " ind = " << "ind << at stmt 97 in main(nmtc)" << G4endl;
	G4Exception("main2");
      }
    case 2:
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    case 3:
      main5();
      break;
    case 4:
      mainBody3c();
    case 5:
      iqhist = 1; 
      iltemp = nhstp;
      ibtemp = nobch + 1;
      --nquit;
      ncol = -4;
      analz1(-4); 
      if (ncmosc != 0) nkey = 4;
      if (iqhist != 1 && nquit > 0) mainBodyInit();
      // ::: if (neutp > 0 && nbertp > 0) endfile neutp
      // ::: if (nhstp > 0 && !bfhist) endfile nhstp
      main1();
    }
  }
  mainBody3b();
}

void G4BertiniCascade::mainBody3b() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  dist = d;
  if (name[no] != nameol) {
    mark = 1;
    nameol = name[no];
  }
  //gomprp(x[no], y[no], z[no], u[no], v[no], w[no], nmed[no] , blz[no], mark, dist, s, xc[no], yc[no], zc[no], ec[no], mat, nameol); //:::
  ind = mark + 3;
  switch(ind) {
  case 1:
    if (verboseLevel > 1) {
      G4cout << "ind  = " << ind << " at stmt 115 in main(nmtc)" << G4endl;
    }
    G4Exception("main2");
  case 2:
    ngam = 1;
    ncol = 4;
    nameol = 0; // escape system
    if (tip[no] != 1.0) {
      dpr = s;
      rr = rng - dpr;
      if (rr <= 0.001) {
	ec[no] = emin[ityp];
	nameol = 0;
      } else {
	rp = rr * rbtors;
	ecol(rp, mat, ityp, ec[no], energy[no]);//:::
	if (verboseLevel > 1) {
	  G4cout << "main after ecol: " << G4endl;
	  G4cout << "rr \t\t"       << rr      << G4endl;
	  G4cout << "rp \t\t"       << rp      << G4endl;
	  G4cout << "rbtors \t\t"   << rbtors  << G4endl;
	  G4cout << "nhist \t\t"    << nhist   << G4endl;
	  G4cout << "no \t\t"       << no      << G4endl;
	  G4cout << "name[no] \t\t" << name[no]<< G4endl;
	  G4cout << "ityp \t\t"     << ityp    << G4endl;
	  G4cout << "ec[no] \t\t"   << ec[no]  << G4endl;
	  G4cout << "ncol \t\t"     << ncol    << G4endl;
	}
      }
    }
    main6();
  case 3:
    ngam = 2;
    ncol = 7; // medium boundary crossing
    if (tip[no] != 1.0) {
      dpr = s;
      rr = rng - dpr;
      if (rr <= 0.001) {
	ec[no] = emin[ityp];
	nameol = 0;
      } else {
	rp = rr * rbtors;
	//ecol(rp, mat, ityp, ec[no], energy[no]);//:::
      }
    }
   main6();
  }
  // magnetic field stuff. straight line approximation
  mainBody3c();  // if ind = 4 mainBody3c() and then to getrig and cascad or to decay with m = 3 if nbeta > 1   
}

void G4BertiniCascade::mainBody3c() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (bodge != 0.0) modify(u[no], v[no], w[no], ec[no], d, tip[no]);
  if (nbeta > 1) {
    ngam = 3;
    ncol = 3;
    m = 3;
    nopart = 0;
    icno = kindo[no];
    if (kindo[no] <= 30) icno = ijevnt[kindo[no]];
    // note: call mfpd2(3) if nbeta > 1, or if ec[no] == emin[ityp]
    //if (ibertp < 0) mfpd2(3);//:::
    mmf = 3;
    main6();
  }
  if (mat == 6666) {
    //dklos();//:::
    lelem = 0;
    nabov = 0;
    nopart = -1;
    npsmd += static_cast<G4int>(wt[no]);
    peiv[ityp] += 1.0;
    mmf = 5;
    m = 5;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    ncol = 5;
    nparto = nopart;
    //    if (ibertp < 0) mfpd2(5); //:::
    if (nseudo > 0) {
      ncol = 5;
      analz1(5); 
      nabov = 0;
      nbelo = 0;
    } else {
      nabov = 0;
      nbelo = 0;
      if (nsav <= 0) {
	nsav = 1;
	es = energy[no];
	xs = x[no];
	ys = y[no];
	zs = z[no];
      }
    }
    if (wt[no] <= wtlow && wt[no] > 0.) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    } if (wt[no] == 0.0) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    }
    mainBody4();
  }
  mainBody5();
  mainBody3a();
}

void G4BertiniCascade::mainBody5() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  getrig();
  G4int test;

  test = main2(); // ::: verify these few lines
  while (test == 0)
    if (test == 1) mainBody2();
  if (test == 2) {
    mmf = 5;
    m = 5;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    ncol = 5;
    nparto = nopart;
    if (ibertp < 0) mfpd2(5);
    if (nseudo > 0) {
      ncol = 5;
      analz1(5); 
      nabov = 0;
      nbelo = 0;
    } else {
      nabov = 0;
      nbelo = 0;
      if (nsav <= 0) {
	nsav = 1;
	es = energy[no];
	xs = x[no];
	ys = y[no];
	zs = z[no];
      }
    }
    if (wt[no] <= wtlow && wt[no] > 0.0) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    }
    if (wt[no] == 0.0) {
      energy[no] = -1.0;
      ic = kindo[no];
      if (ic <= 30) ic = ijevnt[ic];
      mainBody2();
    }
    mainBody4();
  }
  if (test == 3) mainBody8();
  if (test == 4) mainBody9();
  main6();
  mainBody4();
}

void G4BertiniCascade::mainBody6() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (epart[kk][1] > emin[2]) {
    ++nevaph;
    // datahi(); //::: // datahi correctly rotates particles if icem = 1.if icem = 0 it calls gtiso
    if (verboseLevel > 1) {
      G4cout << 
	" nabov = "   << nabov  << " > 350 "  << 
	" (nopart = " << nopart << 
	" num = "     << num    << 
	" nhist = "   << nhist  << ")"<< G4endl;
    } 
    if (kk++ == npart[1]) return;
    mainBody6();
  }
  if (icem == 0) {
    gtiso(adc, bdc, gdc);
    ++nerupl; // nerupl is cumulative sum of npart[1] if  all collisions in a batch
  } else {
    adcc = ang1[kk][1][1];
    bdcc = ang1[kk][2][1];
    gdcc = ang1[kk][3][1];
    sumdc = sqrt(sqr(adcc)+sqr(bdcc)+sqr(gdcc));
    ++nevapl; // nevapl is cumulative sum of npart[1] if all collisions in a batch
    // if icem == 1, must not assume isotropy in lab system, but rotate
    ++n;
    alpha[n] = adcc;
    beta[n] = bdcc;
    gam[n] = gdcc;
    kind[n] = 1;
    ep[n] = epart[kk][1];
    if (verboseLevel > 1) {
      G4cout << "ep of pre-evaporation n = " << ep[n] << G4endl;
    }
    wtfas[n] = 1.0;
    // datalo();//:::
    // note that datalo rotates particles to put them into ub, vb, wb arrays.
    // calling datalo here is a departure from usual procedure and is done
    // only if icem = 1, it is possible to do same if icem = 0, but not necessary 
    // since analysis codes appear to assume o5r neutrons are isotropic
    // adc, bdc, gdc are if double-differential cross-sectiond, and if neutp tape.
    G4double t1 = costh * adcc + sinth * gdcc;
    adc = cosphi * t1 - sinphi * bdcc;
    bdc = sinphi * t1 + cosphi * bdcc;
    gdc = costh * gdcc -  sinth * adcc;
  }
  if (icem != 1) ++neutno;
  // if icem = 1, datalo was called and increases neutno if a neutron.
  // also wtav is incremented in datalo
  edtotn += epart[kk][1];
  if (verboseLevel > 1) {
    G4cout << "edtotn from pre-evaporation" << edtotn << G4endl;
  }
  if (icem != 1) wtav += wtn;
  if (neutp <= 0 || icem == 1) {
    if (kk++ == npart[1]) return;
    mainBody6();
  }
  eneut = epart[kk][1] * 1.0e6;
  if (verboseLevel > 1) {
    G4cout << 
      " nocas  = " << nocas  << 
      " eneut  = " << eneut  << 
      " adc    = " << adc    <<
      " bdc    = " << bdc    << 
      " gdc    = " << gdc    << 
      " xc[no] = " << xc[no] <<
      " yc[no] = " << yc[no] << 
      " zc[no] = " << zc[no] <<  
      " wtn    = " << wtn    << G4endl;
  }
  if (kk++ == npart[1]) return;
  mainBody6();
}

void G4BertiniCascade::mainBody7() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (epart[ll][2] > emin[1]) {
    //datahi();//:::
    if (ll++ == npart[2]) return;
    mainBody7();
  }
  if (icem == 1) {
    adcc = ang1[ll][0][0];
    bdcc = ang1[ll][1][0];
    gdcc = ang1[ll][2][0];
    sumdc = sqrt(sqr(adcc) + sqr(bdcc) + sqr(gdcc));
    G4double t1 = costh * adcc + sinth * gdcc;
    adc = cosphi * t1   - sinphi * bdcc;
    bdc = sinphi * t1   + cosphi * bdcc;
    gdc = costh  * gdcc -  sinth * adcc;
    ++n;
    ep[n] = epart[ll][2];
    if (verboseLevel > 1) G4cout << " ep of pre-evaporated p " << ep[n] << G4endl;
    alpha[n] = adcc;
    beta[n] = bdcc;
    gam[n] = gdcc;
    kind[n] = 0;
    wtfas[n] = 1.0;
    // datalo();//:::
  }
  if (ll++ == npart[2]) return;
  mainBody7();
}

void G4BertiniCascade::mainBody8() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  kstor = 1;
  if (verboseLevel > 1) G4cout <<  "nopart = " << nopart << " lelem = " << G4endl;
  for (G4int i = 0; i < nopart; ++i) {
    n = i;
    if (verboseLevel > 1) G4cout << "kind = " << kind[n] << " ep = " << ep[n] << G4endl;
    nnow = n;
    ic = kindi[n];
    if (ic <= 30) ic = ijevnt[ic];
    if (verboseLevel > 1) G4cout << " wtfas < wtlow" << G4endl;
    if (ep[n] <= emin[kind[n]+1]) {
      // datalo();//:::
      if (kind[n] == 1) {
	++no5rca;
	edtotn += ep[n];
	if (verboseLevel > 1) G4cout << " edtotn from cascad = " << edtotn << G4endl;
	if (kind[n] != 4) continue;
	ep[n] = 1.0;
	// note that this pi- is written into both the datahi and the
	// datalo arrays. it will not matter if one does not try to add
	// the two arrays (tipa and tipb) to check if energy conservation
      }
      //datahi();//:::
    } 
    if (nabov > MAXP) {
      if (verboseLevel > 1) {
	G4cout << 
	  " nabov  : " << nabov  << 
	  " nopart : " << nopart << 
	  " num    : " << num    << 
	  " nhist  : " << nhist  << G4endl;
      }
    }
  }
}

void G4BertiniCascade::mainBody9(){}; // evaporation and pre-equilibrium

G4int G4BertiniCascade::main2() {
  noo = no;
  if (wt[no] <= wtlow) {
    energy[no] = -1.0;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    return 1;
  }
  cascad();
  if (nopart > 0) {
    if (ipcl == 0 && iqheh == 0) ++ncasca;
  } else if (nopart < 0) {
    npsm += static_cast<G4int>(wt[no]);
    if (energy[no] != 1.0) return 2;
    return 0;
  }
  if (nsav > 0) {
    nsav = 0;
    energy[no] = es;
    x[no] = xs;
    y[no] = ys;
    z[no] = zs;
  } if (nopart > 0) return 3;
  return 4;
}

void G4BertiniCascade::main6() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (nsav > 0) {
    nsav = 0;
    energy[no] = es;
    x[no] = xs;
    y[no] = ys;
    z[no] = zs;
  }
  // ncol = 3 was set at call mfpd2(3) if decay and stopped particle
  analz1(ncol); 
  nabov  =  0;
  nbelo  =  0;
  if (ngam < 0) {
    energy[no] = -1.0;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    mainBody2();
  } else if (ngam == 0) {
    mainBody4();
  }
  if (tip[no] != 4.0) {
    energy[no] = -1.0;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    mainBody2();
  }
  if (npidk > 0) {
    energy[no] = -1.0;
    ic = kindo[no];
    if (ic <= 30) ic = ijevnt[ic];
    mainBody2();
  }
  ec[no] = 1.0;
  energy[no]  = 1.0;
  oldwt  =  wt[no];
  x [no]  = xc[no];
  y [no]  = yc[no];
  z [no]  = zc[no];
  // now, send pi- at rest to getrig and cascad(collision).
  coupi += wt[no];
  mainBody5();
}

void G4BertiniCascade::mainFinal() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  ncas += nocas;
  ktim = iclock();
  jtim = (ktim - itim) / 100;
  jtime = jtim - jtimo;
  xtim = (ktim - itim) / 100.0;
  xtime = xtim - xtimo;
  jtimo = jtim;
  xtimo = xtim;
  if (nbertp > 0) {
    nn = -1;
    if (neutp > 0) {
      if (verboseLevel > 1) {
	G4cout << "nn    : " << nn    << G4endl;
        G4cout << "edum  : " << edum  << G4endl;
	G4cout << "udum  : " << udum  << G4endl;
	G4cout << "dum   : " << vdum  << G4endl;
	G4cout << "xdum  : " << xdum  << G4endl;
	G4cout << "ydum  : " << ydum  << G4endl;
	G4cout << "zdum  : " << zdum  << G4endl;
	G4cout << "wtdum : " << wtdum << G4endl;
      }
    }
  }
  wttot = wtav;
  cneut = neutno;
  if (neutno > 0) {
    wtav = wttot / cneut;
    cneuav += cneut;
    coupav = ncoutp[0] / cneut;
  } else {
    wtav = 0.0;
    wttot = 0.0;
  }
  if (verboseLevel > 1) {
    G4cout << "end of batch                                                 : " << nobch     << G4endl;
    G4cout << "neutrons produced in this batch                              : " << neutno    << G4endl;
    G4cout << "cumulative number of cascades completed                      : " << ncas      << G4endl;
    G4cout << "cpu time (s)                                                 : " << xtime     << G4endl;
    G4cout << "number of qevent or qheh collision number                    : " << num       << G4endl;
    G4cout << "nhist                                                        : " << nhist     << G4endl;
    G4cout << "total weight of o5r neutrons                                 : " << wttot     << G4endl;
    G4cout << "average weight of o5r neutrons                               : " << wtav      << G4endl;
    G4cout << "number of real collisions (non-hydrogen)                     : " << ncasca    << G4endl;
    G4cout << "number of pre-equilibrium and evaporation neutrons in datahi : " << nevaph    << G4endl;
    G4cout << "number of pre-equilibrium and evaporation (erupcem) neutrons : " << ncoutp[0] << G4endl;
    G4cout << "number of pre-equilibrium and evaporation o5r neutrons       : " << nevapl    << G4endl;
    G4cout << "number of evaporation (non-cem) o5r neutrons                 : " << nerupl    << G4endl;
    G4cout << "number of o5r neutrons from sub cascad                       : " << no5rca    << G4endl;
  }
  // why is ncoutp[0] often less than nevapl? What is ncoutp[0]?
  nspltt = 0;
  ncol = -3;
  analz1(ncol);
  ncasca = 0;
  no5rca = 0;
  nevaph = 0;
  nevapl = 0;
  nerupl = 0;
  wtav   = 0.0;
  wttot = 0.0;
  for (G4int i = 0; i < 9; ++i) ncoutp[i] = 0;
  nomaxn = nomax;
  if (nobch == maxbch && nocas == maxcas) {
    cnhist = maxbch * maxcas;
    cnhis = nhist;
    if (verboseLevel > 1) {
      G4cout << "cnhist : " << cnhist << G4endl;
      G4cout <<  "cneuav : " << cneuav << G4endl;
      G4cout << "edtotn : " << edtotn << G4endl;
    }
    if (cnhist > 0.0) {
      cneuav /= cnhist;
      edtotn /= cnhist;
    }
    edo5r = 0.0;
    if (cneuav > 0.0) edo5r = edtotn/cneuav;
    if (verboseLevel > 1) {
      G4cout << "average number of o5r neutrons per history       : " << cneuav << G4endl;
      G4cout << "total number of histories                        : " << cnhist << G4endl;
      G4cout << "average energy of o5r neutrons per history (MeV) : " << edtotn << G4endl;
      G4cout << "average energy per o5r neutron                   : " << edo5r  << G4endl;
    }
  }
  neutno = 0;
  nocas = 1;
  if (++nobch <= maxbch) mainBody();
  --nquit;
  ncol = -4;
  analz1(-4);
  if (ncmosc != 0) nkey = 4;
  if (iqhist != 1 && nquit > 0) mainBodyInit();
  //if (neutp > 0 && nbertp > 0) endfile neutp; //:::
  //if (nhstp > 0 && !bfhist) endfile nhstp; //:::
  main1();
}

void G4BertiniCascade::hcol(G4int ib, 
		     G4int ityp,  
		     G4int hsig, 
		     G4int no1,  
		     G4double ec, 
		     G4int nopart, 
		     G4int *kind,  
		     G4double *ep,  
		     G4double *alf, 
		     G4double *bet,  
		     G4double *gam){ 
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int j; 
  G4int iqho;
  G4int ndiff;
  G4int ippo;
  G4int ii;
  G4int noadd;
  G4int i1i; 
  G4int i2i;
  G4int noparo;
  G4int nopnow;
  G4int nopdif;
  G4int noopar;
  G4int itinc;
  G4int icc[12];
  G4int itt[MAXP];
  G4int itp;
  G4int kno; 
  G4int igo;
  G4double esum;

  iqh = 1;
  iqk = 0;
  if (ifirst == 0) {
    for (G4int i = 0; i < 161; ++i) {
      // cdd = frinn[i] * 100.0; //:::use frinn owned by BertiniModel
      ii = static_cast<G4int>(cdd + 1.0);
    }
  }
  ++ifirst;
  icoun1 = 0;
  icoun2 = 0;
  icoumx = 50;
  // icoumx above is optional iteration count in subroutine 
  ndiff = 0;
  ke = 1;
  if (ib == 0) {
    ib = 1;
    sf = 1.0;
    i1i = 0;
    i2i = 599; // '600' in f77
    for (G4int j = 0; j < 4; ++j) {
      for (G4int i = i1i; i < i2i; ++i) {
	// read(nbertp)(crsc(i),i=i1i,i2i) //:::
      }
      i1i += 600;
      i2i += 600;
    }
    for (G4int i = 0; i < 29849; ++i) {
      //  f77: 'read(nbertp) (tapcrs(i),i=1,29849)' //:::
    }
    randi[0] = 16896;
    sqnm = sqr(massNucleon); //  
    rcpmv = 5.0613e10; // reciprocal cm per MeV  
    for (i = 0; i < 19; ++i) {
      poac[i] += poac[i];
      ppac[i] /= sf;
    }
    // poac(19),ppac(19)
  } 
  for (G4int i = 0; i < 60; ++i) ipec[i] = 0;
  for (i = 0; i < 2114; ++i) esps[i] = 0.0;
  // do 19 i = 4515,4849 
  // 4855 is correct for cray unicos version if G4int*4 
  // for random nos.rand,erand,randi, and double precision turned off 
  // for IBM with G4int*2 for rands and erand,4549 is right,not 4548; 
  // for epmnas cluster, with G4int*4 for rands,erand,and randi: 485 
  for (i = 4514; i < 4851; ++i) esps[i]  = 0.0;
  for (i = 0;    i < 12;   ++i) icc[i]   = 0;
  for (i = 0;    i < 4;    ++i) rands[i] = randi[i];
  massParticle[0] = massNucleon;

  // p  n  pi+ pi0 pi- mu+ mu- 
  if (ityp == 4 || ityp == 6 || ityp == 7 ) {
    // dolio(c3, c1, (char *)(*ityp), (ftnlen)sizeof(G4int)); generated by f2c?
    // 24 format(1h ,7hityp = ,i2,19h at stmt 23 in pcol) 
    G4cerr << "hcol1" << G4endl;
  }
  if  (ityp == 3 || ityp == 5 ) {
    massParticle[0] = massPionCharged;
  }
  inc = 1;
  clsm = 1.0;
  massParticle[1]   = massNucleon;
  energy[1]    = massNucleon;
  energy[0]    = ec * rcpmv + massParticle[0]; 
  pxyz[1] = 0.0;
  pxyz[5] = 0.0;
  pxyz[9] = 0.0;
  p1cli(); // calc's x,y,z mom coord's of inc. particle (pxyz(1-5-9)) and p1oe1 = pz1 
  rlke = (((energy[0] * energy[1] - pxyz[0] * pxyz[1] - pxyz[4] * pxyz[5] - pxyz[8] * pxyz[9]) / massNucleon) - massParticle[0]) / rcpmv;
  G4double r = G4UniformRand();
  G4double rinit = r;
  value1 = rlke;
  itp = ityp - ityp / 5;
  nopart = 2;
  switch (itp) {
  case 1:  goto L200;
  case 2:  goto L200;
  case 3:  goto L260;
  case 4:  goto L150;
  }
 L150:
  r -= xsec(4, itp, ec) / hsig; // try (pi-, p) exchange 
  if (r > 0.0) goto L260;       // a (pi-, p) exchange event has occurred 
  it = 5;
  kind[1]   = 3;
  kindi[0]  = 4;
  ibbarr[0] = 0; // pi0 
  ibbarr[1] = 1;
  kindi[1]  = 2;
  kind[2]   = 1; // neutron 
  pt[1]     = 4.0;
  pt[13]    = 2.0;
  ik = it;
  massParticle[2] = massPionZero;
  if (rlke <= 340.0) {
    // snt = mud(); //:::
    // crdet(51, pmdx, rlke);        // :::dif (pi-, p) exchange 
    goto L500;
  }
  i3 = 1;
  if (rlke < 1000.0) {
    i3 = 2;
    value1 = rlke - 340.0;
  }
  rou16(pmdx);                 // dif (pi-, p) exchange 2500-1000, 1000-340, 340- 0 
  goto L570;
 L200:
  rinit = r;
  ndiff = 0;
  iqho = iqh;
  ippo = ippp;
  r = rinit; 
  rr = r - xsec(2, itp, ec) / hsig; // try (N, p) double production 
  if (iqh == 0) goto L230;
  if (ec < eth[2][itp]) goto L260;
  r -= xsec(2, itp, ec) / hsig;
  if (r > 0.0) goto L260; // double production has occurred
  nopart = 4;
  ippp = 0;
 L230:
  nno = ityp; 
  e1nc = ec;
  if (iqh == 1) goto L270;
  ndiff = noparo - nopart;
  // nopart=3 
  if (nopart != 3) nopart = 3;
  // we have assumed single isobar production from isob;isob does 
  // not return nopart or the equivalent, unless it is dout(1). 
  //isob();//:::
  iqk = 0;
  if (icoun1 > icoumx || icoun2 > icoumx) nopart = -1;
  iqh = iqho;
  if (nopart == -1) return;
  nopnow = static_cast<G4int>(dout[0]);
  nopdif = nopnow - nopart;
  einit = MeV * (ginum[ityp - 1] + ginum[0]) + ec;
  noadd = 0;
 L240:
  esum = 0.0;
  iqk = 0; 
  kno = 2;
  i1 = nopart;
  G4double edhcol;
  G4int noo;
  for (noo = 0; noo < i1; ++noo) {
    kind[noo] = static_cast<G4int>(dout[kno    ]) - 1;
    ep[noo]   = static_cast<G4int>(dout[kno + 1]);
    alf[noo]  = static_cast<G4int>(dout[kno + 4]);
    bet[noo]  = static_cast<G4int>(dout[kno + 5]);
    gam[noo]  = static_cast<G4int>(dout[kno + 6]);
    kindi[noo] = kind[noo] + 1;
    ibbarr[noo] = 0;
    if (kindi[noo] < 2) ibbarr[noo] = 1;
    amm = ginum[kind[noo]];
    amm *= 1e3;
    if (ep[noo] > 0.0) esum += ep[noo] + amm;
    edhcol = einit - esum;
    if (ep[noo] <= 0.0) noopar = noo;
    kno += 7;
  }
  if (edhcol > 2.0  & noadd == 0) {
    nopart = 4;
    noadd = 1;
    goto L240;
  } else nopart = 3;
  if (ep[1] <= 0.0) nopart = -1;
  if (nopart == -1) return;
  goto L430;
 L260:
  //  try sngl productio  ::: verify
  if (ec < eth[0][itp]) goto L330;
  r -= xsec(1, itp, ec) / hsig;
  if (r > 0.0) goto L330; // (N, p), (pi+, p), or (pi-, p)  sngl production has occurred 
  nopart = 3;
  ippp = -1;
  goto L230;
 L270:
  itinc  = ityp - 1;
  eoprs  = ec;
  noparo = nopart;
  no     = ityp;
  entr   = 0.0;
  nofsk  = 0;
  coqheh += 1.0;
  ndiff  = noparo - nopart;
  einit  = 1000 * (ginum[0] + ginum[ityp - 1]) + ec;
  esum   = 0.0;
  // go to 4113 
  // above statement commented out 2-9-98 because tagsgi says that the 
  // next statement cannot be reached. 
  if (nopart > 0) {
    i1 = nopart;
    for (G4int jj = 0; jj < i1; ++jj) {
      j = jj;
      kind[j] = itt[j];
      ep[j] = eff[j]; 
      amm = MeV * ginum[kindi[j]];
      esum += amm + ep[j];
      alf[j] = alpp[j];
      bet[j] = bett[j];
      gam[j] = gamm[j];
    }
  }
  eloss = einit - esum - entr;
  eloss1 = eloss;
  if (nopart > 0) {
    esum = 0.0;
    for (G4int jj = 0; jj < nopart; ++jj) {
      j = jj;
      kind[j] = itt[j];
      ep[j]   = eff[j];
      amm     = ginum[kind[j]] * MeV;
      esum   += amm + ep[j];
      alf[j]  = alpp[j];
      bet[j]  = bett[j];
      gam[j]  = gamm[j];
    }
  }
  eloss = einit - esum - entr;
  iqh = iqho;
  iqk = 1;
  return;
 L330:
  //  an elastic event has occurred. 
  switch (ityp) {
  case 1:  goto L350;
  case 2:  goto L360;
  case 3:  goto L460;
  case 4:  goto L340;
  case 5:  goto L520;
  case 6:  goto L340;
  case 7:  goto L340;
  }
 L340:
  //::: f2c? dolio(c3, c1, (char *)(*ityp), (ftnlen)sizeof(G4int));
  // 602 format(1h ,6hityp =,i2,17h stmt 601 in hcol) 
  G4cerr << "hcol2" << G4endl;
  exit(0);
 L350:
  i3 = 7;
  it = 18;
  goto L370;
 L360:
  i3 = 4;
  it = 15;
 L370:
  //rou20(dcin, dcln, dchn, pdci, pdch); //:::
  // elas dif xsects-np 300-740, np 0-300, np 660-3500, pp 500-1000,pp 660-3 
  // sets pt(2)=1.-prot,2.-neut,pt(14)=1.,pm(3)= massNucleon 
  // calc's cst(scat cos) for ec gt (500 for prot,740 for neuts)
  igo = i3 - 4;   // i3 = 5 6 7 8 
  switch (igo) {
  case 1:  goto L450;
  case 2:  goto L380;
  case 3:  goto L390;
  case 4:  goto L400;
  }
 L380:
  snt = sqrt(1.0 - sqr(cst));
  goto L400;
 L390:
  // pol1(cst, snt); //:::
 L400:
  // ::: scatteringWithHydrogen(); // elastic
  if (it != 5) {
    kind[1] = ityp - 1;
    kindi[0] = kind[1] + 1;
    ibbarr[0] = 1;
    if (ityp > 2  & ityp < 8) ibbarr[0] = 0;  
    kind[1] = 0;
    kindi[0] = kind[0] + 1;
    ibbarr[1] = 1;
  }
  // p 
  ep[1]  = static_cast<G4double>(pt[2]);
  alf[1] = static_cast<G4double>(pt[7]);
  bet[1] = static_cast<G4double>(pt[8]);
  gam[1] = static_cast<G4double>(pt[9]);
  ep[2]  = static_cast<G4double>(pt[14]);
  alf[2] = static_cast<G4double>(pt[19]);
  bet[2] = static_cast<G4double>(pt[20]);
  gam[2] = static_cast<G4double>(pt[21]);
  // ea = ec * rcpmv; //::: fix clash against ea
 L430:
  for (i = 1; i < 4; ++i) randi[i] = rands[i];
  return;
 L450:
  i3 = 4;
  //rou15(); //:::
  // calc's cst(scat cos) for ec < 740 MeV 
  goto L380;
 L460:
  i3 = 1;
  it = 1;
  // (pi+, p) elastic 
  //rou11(); //:::
  // tapcrs(248) (pi+, p) <  340 MeV 
  switch (i3) {
  case 1:  goto L510;
  case 2:  goto L480;
  case 3:  goto L480;
  case 4:  goto L500;
  }
 L480:
  //::: dolio(c3, c1, (char *)i3, (ftnlen)sizeof(G4int));
  // 807 format(1h ,3hi3=,i3,20h at stmt 806 in hcol) 
 L490:
  G4cerr << "hcol3" << G4endl;
 L500:
  cst = crdt[1] - fabs(snt * (crdt[1] - crdt[0]));
  goto L380;
 L510:
  i3 = 1;
  //rou15(); //:::
  // (pi+, p) dif elastic for 1.0 - 2.5 GeV at 3401, and for 0.34 -1.0 GeV at 3231 
  i3 = i3;
  goto L380;
 L520:
  it = 3;
  // (pi-, p) elastic 
  i3 = 2;
  pt[1] = 5.0;
  pt[13] = 1.0;
  ik = it;
  massParticle[2] = massPionCharged;
  if (rlke > 340.0) goto L540;
  // snt = mud(); //:::
  // crdet(51, pmdd, rlke); //:::
  goto L500;
 L540:
  i3 = 1;
  if (rlke < 1000.0) {
    i3 = 2;
    value1 = rlke - 340.0;
  }
  rou16(pmdd);
 L570:
  i3 = i3;
  if (i3 < 0) goto L380;
  if (i3 > 0) goto L590;
  i3 = 3;
  //rou15(); //:::
  goto L380;
 L590:
  //::: dolio(c3, c1, (char *)i3, (ftnlen)sizeof(G4int));
  // 980 format(1h ,3hi3=,i3,20h at stmt 970 in hcol) 
  goto L490;
}

void G4BertiniCascade::cascad() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int l;
  G4double cou12(0.0);
  G4double coupic(0.0);
  G4double am1(940.0746);
  G4double am1sq(883740.253565);
  G4int ifirst(1);
  G4int iwr(0);
  G4int ibrcon(0);
  static G4bool initi(true);
  if (initi) {
    initi = false;
    echek  = 0.0;
    echekk = 0.0;
    echec  = 0.0;
    npseub = 0.0;
    npsub  = 0.0;
    nozerb = 0.0;
    cobert = 0.0;
    cobcem = 0.0;
    counsk = 0.0;
    coqevt = 0.0;
    npseus = 0.0;
    nozers = 0.0;
    nskcol = 0.0;
    npseuq = 0.0;
    nozerq = 0.0;
    pdevst = 0.0;
  }
  nok   = 0;
  iqev  = 0;
  ibrt  = 0;
  iqheh = 0;
  iheh  = 0;
  iskl  = 0;
  ipcl  = 0;
  icr   = 0;
  iphev = 0;
  eno   = energy[no];
  tipno = tip[no];
  ecno  = ec[no];
  icem  = 0;
  ncc   = maxbch * maxcas;
  nhist = (nobch - 1) * maxcas + nocas;
  if (nhist != nhist1 && nhcasc == 1) {
    if (verboseLevel > 1) {
      G4cout << "cascad " << G4endl;
      G4cout << "nhist  : " << nhist  << G4endl;
      G4cout << "nhist1 : " << nhist1 << G4endl;
      G4cout << "nobch  : " << nobch  << G4endl;
      G4cout << "nocas  : " << nocas  << G4endl;
      G4cout << "maxcas : " << maxcas << G4endl;
      G4cout << "need to zero nhist in block data and put nhist = nhist+1 after first return in sors " << G4endl;
      G4cout << "need common numprt / in sors. and in block data "<< G4endl;
      G4cout << "nhistt : " << nhistt << G4endl;
      G4cout << "nhcas  : " << nhcas  << G4endl;
      G4cout << "nhcasc : " << nhcasc << G4endl;
    }
  }
  nofask = 0;
  noo    = no;
  nooo   = no;
  nonow  = no;
  io     = iou;
  tipp   = tip[no];
  ectipp = ec[no];
  negzpr = 0;
  ehicut = ehipi;
  if (tip[no] <= 1.0) ehicut = ehin;
  ihie = 1;
  if (ec[no] <= ehicut) ihie = 0;
  itim = iclock();
  r = G4UniformRand(); 
  if (energy[no] == 1.0) { 
 // pi- capture

  }  else { // ::: verify
    dklos();
    // if (ce first collision with ist  material  in mat = 1
    // provided particle enters in mat == 1 (set in sors and governed
    // by x(1), y(1), z(1) and the input geometry)
    lelem = 1;
    if (ifircl >= 1 && name[no] == 1 && mat == 1) { // allow if (large pseudo-cross-section if (mag field-bodge
      r1 = r;
      flt = G4UniformRand();
      rat = bodge / delsig;
      r = flt;
      if (flt < bodge / delsig) {
	++pevts[mat];
	lelem = 0;
	nopart = -1;
	goto L5;
      }
      delsig -= bodge; // now back to normal
      if (ityp >= 6) {  // end elastic additions
	// pseudo collision //:::
	goto L5;
      }
      r -= hsigg[ityp][mat] / delsig;
      // note that denh[mat] = 0.0 means hsigg[ityp][mat] = 0.0 so cannot goto L3 unless r = 0.0; correct if this case
      if (denh[mat] > 0.0 && r <= 0) {
	if (ihie != 1) {
	  gthsig(2);
	  r = G4UniformRand();
	  if (r > (hsig * denh[mat] * 1.0e+24 / hsigg[ityp][mat])) {
	    // speudo collision 2
	    // :::
	    goto L5;
	  }
	  ipcl   = 1;
	  esum   = 0.0;
	  edpcol = 0.0;
	  entrr  = 0.0;
	  nofsk  = 0;
	  no1    = no;
	  weni   = wt[no];
	  ihit   = 1;
	  if ((energy[no] == 1.0 || ec[no] == 1.0) && tip[no] == 4) ++coupih;
	  // in either pcol or hcol, the target nucleon is a proton
	  if (ihcol == 1)
	    hcol(ibert, ityp, hsig, no1, ec[no], nopart, kind, ep, alpha, beta, gam);
	  else
	    anucco = 1.0;
	  znucco = 1.0;
	  ++copcol;
	  if (nopart > 0) {
	    for (G4int i = 0; i < nopart; ++i) {
	      if (ihcol == 0) ++kindi[i];
	      if (ic <= 30) ic = ijevnt[kindi[i]];
	      else ic = kindi[i];
	      ibb = iibar[ic];
	      ibbarr[i] = ibb;
	      wtfas[i] = wt[no1];
	    }
	    lelem = -1;
	    ++rhevts[mat];
	  } else  { // if nopart <= 0
	    ++coupsd;
	    // pseudo collision 2
	    // :::
	  }
	  goto L5;
	}
	inc    = ityp - 1; // get h cross-section at ec
	icr    = 1;
	esum   = 0.0;
	nopart = 0;
	nofask = 0;
	einit  = 0.0;
	r      = G4UniformRand();
	anumt  = heht * denh[mat] * 1.0e24;
	denn   = hsigg[ityp][mat];
	rat    = anumt / denn;
	dif    = rat - r;
	if (dif < r) ++coupsd;
	if (heht * denh[mat] * 1.0e24 / hsigg[ityp][mat] < r) {
	  ++phevts[mat];
	  if (ihie == 1) ++hpevth[mat];
	  iphev  = 1;
	  lelem  = 0;
	  nopart = -1;
	  goto L5;
	}
	++hrevth[mat]; // real event
	r = G4UniformRand();
	rath = hehel / heht;
	if (r < (hehel / heht))  { // elastic
	  edsc = 0.0;
	  icr = 2;
	  esum = 0.0;
	  ihit = 1;
	  //scatt(inc); //:::  scatt(inc, ec[no], kind, ep, alpha, beta, gam); 
	  nopart = 2;
	  wtfas[1] = 1.0;
	  wtfas[2] = 1.0;
	  // if ihie == 1 sub datahi multiplies wtfas by wt[no] to get weight
	  kindi[1]  = kind[1] + 1;
	  kindi[0]  = kind[0] + 1;
	  ibbarr[1] = iibar[kindi[1]];
	  ibbarr[0] = iibar[kindi[0]];
	  ++coelas;
	  lelem = -1;
	  ++rhevts[mat];
	  goto L5;
	}
	// nonelastic
	if (icon <= 1) ikin = -1;
	else if (icon == 2) ikin = kindo[no];
	if (ikin > 5 || ec[no] > eheh[static_cast<G4int>(tip[no]+0.01) + 1]) {
	  weni   = wt[no];
	  iqheh  = 1;
	  esum   = 0.0;
	  ihit   = 1;
	  edqheh = 0.0;
	  entrr  = 0;
	  nofsk  = 0;
	  nofask = 0;
	  apr    = 1;
	  zpr    = 1;
	  ++coqheh;
	  if (nopart < 0) ++coupsd;
	  else if (nopart > 0) {
	    for (G4int i = 0; i < nopart; ++i) esum +=  ginum[kind[i] + 1] * MeV+ ep[i];
	    // kindi[i] has been defined in qheh, however, ep[i] has been computed if kind[i]
	  }
	  einit = ec[no] + ginum[ityp] * MeV + ginum[ihit] * MeV;
	  edqheh = einit - esum;
	  // note that mass and kinetic energy of all leptons is omitted here!
	  // the same is true if (all particles sent to addnwk (included in nofask)
	  lelem = -1;
	  ++rhevts[mat];
	  goto L5;
	}
	// note: only got here ihie  == 1,so ec[no] > eheh(index+1) if
	//       eheh  <= ehicut as set in subroutine sors.
	//       routine heh is dummied out in this version of code
	iheh   = 1;
	ihit   = 1;
	nopart = 0;
	++couheh;
	if (nopart > 0) {
	  esum = 0.0;
	  for (G4int i = 0; i < nopart; ++i) {
	    kindi[i] = kind[i] + 1;
	    (kindi[i] > 2) ? ibb = 0 : ibb = 1;
	    ibbarr[i] = ibb;
	    esum += ep[i] + ginum[kindi[i]] * MeV; //::: move 1000 to MeV
	    einit = ec[no] + ginum[ityp] * MeV + ginum[ihit] * MeV;
	  }
	}
	lelem = -1;
	++rhevts[mat];
      L5:
	apr  = 0.0;
	zpr  = 0.0;
	ex   = 0.0;
	erec = 0.0;
	for (G4int i = 0; i < 6; ++i) npart[i] = 0;
	uu = 0.0;
	entr = 0.0;
	if (nofask > 0) {
	  for (G4int j = 0; j < nofask; ++j) {
	    entr += ginum[static_cast<G4int>(kindk[j] + 1.0)] * MeV + epk[j];
	    ic <= 30 ? ic = ijevnt[ic] : ic = kindk[j] + 1;
	  }
	}
	eloss = 0.0;
	if ((nopart > 0 || nofask > 0) && icr != 1) eloss = einit - esum - entr;
	return;
      }
      nnn = nel[mat]; // here, the choice of isotope(lelem) if given material is made
      for (G4int lem = 0; lem < nnn; ++lem) {
	lelem = lem;
	r -= sigg[lem][mat] / delsig;
	if (r <= 0) goto L12;
      }
      rcut = sigg[lelem][mat] / delsig;
      rdif = r - rcut;
      // revisions to routine cascad to implement neutron-nucleus elastic scat
      // add following cards after ftn stmt no. 11 isn 0038 on 68.122/14.19.4
      if (ityp != 2 || ec[no] > elas) { // end elastic additions
	++pdevts[mat];
	++pdevst;
	if (tip[no] > 1.0) {
	  if (ec[no] > ehipi) ++hpdevt[mat];
	} else {if (ec[no] > ehin) ++hpdevt[mat];
	}
	lelem  = 0;
	nopart = -1;
	goto L5;
      }
      nl = noel[mat];
      if (nl <= 0) { // end elastic additions
	++pdevts[mat];
	++pdevst;
	if (tip[no] > 1.0) {
	  if (ec[no] > ehipi) ++hpdevt[mat];
	} else {
	  if (ec[no] > ehin) ++hpdevt[mat];
	}
	lelem = 0;
	nopart = -1;
	goto L5;
      }
      for (G4int l = 0; l < nl; ++l) {
	lelem = l;
	r -= sgels[l] / delsig;
	if (r <= 0.0) {
	  idd = id[l][mat];
	  // interpolateElasticNeutronData(mat, idd, ec[no]); // gets fone //:::
	  ex = 0.0;
	  kind[0] = ityp - 1;
	  uu = 0.0;
	  apr = a[lelem][mat];
	  zpr = zz[lelem][mat];
	  for (G4int i = 0; i < 6; ++i) npart[i] = 0;
	  csthcm = 2.0 * G4UniformRand() - 1.0;
	  r1 = G4UniformRand();
	  ab3f1 = fabs(3.0 * fone);
	  if (ab3f1 > r1) {
	    if (ab3f1 <= 1.0 || ab3f1 <= 2.0 * G4UniformRand() + 1.0) {
	      cos1 = 2.0 * G4UniformRand() - 1.0;
	      if (cos1 >= csthcm) csthcm = cos1;
	      if (fone < 0.0) csthcm *= -1;
	    } else {
	      if (fone >= 0.0) {
		++peevts[mat];
		nopart = -1;
		goto L5;
	      }
	      csthcm = -1.0;
	    }
	  }
	  goto L30;
	}
      }
      idd = id[nl][mat];
      // interpolateElasticNeutronData(mat, idd, ec[no]); //:::
      ex = 0.0;
      kind[0] = ityp - 1;
      uu = 0.0;
      apr = a[lelem][mat];
      zpr = zz[lelem][mat];
      G4int i;
      for (i = 0; i < 6; ++i) npart[i] = 0;
      csthcm = 2.0 * G4UniformRand() - 1.0;
      r1 = G4UniformRand();
      ab3f1 = fabs(3.0 * fone);
      if (ab3f1 > r1) {
	if (ab3f1 <= 1.0 || ab3f1 <= 2.0 * G4UniformRand()+ 1.0) { //:::
	  cos1 = 2.0 * G4UniformRand() - 1.0;
	  if (cos1 >= csthcm) csthcm = cos1;
	  if (fone < 0.0) csthcm *= -1;
	} else {
	  if (fone >= 0.0) {
	    ++peevts[mat];
	    nopart = -1;
	    goto L5;
	  }
	  csthcm = -1.0;
	}
      }
    L30:
      ++reevts[mat];
      am2 = a[lelem][mat] * amu_c2;
      e1 = ec[no] + am1;
      e2 = am2;
      wlab = e1 + e2; // am1, e1 and am2, e2 are rest mass, total energy if (part. 1 (inc. neut.)
      // and part. 2(struck nucleus) respectively, wlab is total energy in lab i
      p1 = sqrt(e1* e1 - am1sq); // total momentum of particle 1 in lab,  p2 = 0.0
      bta = p1 / wtlab;
      zta = sqrt(1.0 - sqr(bta));
      gma = 1.0 / zta;
      wcm = zta * wlab; // particle 4 is secondary n particle 3 is recoiling N
      e4cm = am2 * gma; // this expression if e4cm is obtained from equation 29b, rev. mod. phys., vo no. 3,
      // july 1962, page 433,  kinematics of high energy particles, k. g.
      e3cm  = wcm - e4cm;
      p3cm  = sqrt(e3cm * e3cm - am1sq);
      pz3cm = p3cm * csthcm;
      pz3   = gma * (pz3cm + bta * e3cm);
      e3    = gma * (e3cm + bta * pz3cm);
      p3    = sqrt(sqr(e3) - am1sq);
      ep[1] = e3 - am1;
      csth  = pz3 / p3;
      e4    = wlab - e3;
      if (erec < 0.0) erec = 0.0;
      nopart = 1;
      snth = 1.0 - sqr(csth);
      if (snth < 0.0) {
	csth = 1.0;
	snth = 0.0;
      } else snth = sqrt(snth);
      azirn(snph, csph);
      alpha[0] = snth * csph;
      beta[0]  = snth * snph;
      gam[0]   = csth;
      kindi[0] = 2;
      ibbarr[0]= 1;
      return;
    }
  }
 L12:
  f[0] = a[lelem][mat];
  f[1] = zz[lelem][mat];
  f[2] = ec[no];
  f[3] = ctofe;
  f[4] = 1.0;
  f[5] = andit;
  f[6] = tip[no];
  f[7] = ctofen;
  j = 1;
  // NOTE the variable andit is not initialized in sub input(common/inpu/)
  // but is never actually used to determine the angular
  // distributions in the Sternheimer-Lingenfelter isobar model
  denlm =   den[lelem][mat];
  sigglm = sigg[lelem][mat];
  lele  = lelem;
  matt  = mat;
  anucc = f[1];
  znucc = f[2];
  icem = 0; // dummy routine mashnk(icem, anucc, znucc, ecno, tipno, ehipi, ehin, ehicut);
  if (tip[no] > 1.0) {
    if (ec[no] > ehipi) {
      ehicut = ehipi;
      ihie = 1;
      if (icem == 1) {
	ihie = 0;
	goto L430;
      }
      goto L440;
    }
    ihie = 0;
    goto L430;
  }
  if (ec[no] <= ehin) {
    ihie = 0;
    goto L430;
  }
  ehicut = ehin;
  ihie = 1;
  if (icem == 1) {
    ihie = 0;
    goto L430;
  }
  goto L440;
 L430:
  ehinn = ehin;
  if (f[7] > 1) ehinn = ehipi;
  eskal = eskale[static_cast<G4int>(f[7]) + 1]; //:::
  if (eskal < ehinn && ec[no] > eskal) ihie = 1;
  if (ec[no] == 1.0 || energy[no] == 1.0) ihie = 0;
  if (kindo[no] > 7 && (ifircl >= 1 && name[no] == 1 && mat == 1)) ihie = 1;
  if (icem == 1) ihie = 0;
  if (ihie == 1 || (ihecc == 1 && icem == 0)) goto L4440;
  if ((energy[no] == 1.0 || ec[no] == 1.0) && tip[no] == 4) ++counpi;
  // i.e. will call qevent which calls eventq which will call nucrin
  // but, if ihecc = 1 and f(3) < eskal, qevent will call hecc!
 L429:
  ibrt = 1;
  esum = 0.0;
  ++iwr;
  anuma = energyAZ(anucc, znucc) / MeV; 
  anumas = anuma + conver * anucc;
  no1 = no;
  if (imash == 1 && icem == 0) {
    ++ibrcon;
    if (ibrcon == 1) ibert = 0;
  } 
  if (icem == 1) {
    // bertcem(ib, f[1], nopart, kind, ep, alpha, beta, gam);//::: 
    if (mat == 1) cobcem += wt[no];
    cobc += wt[no];
  } else {
    // bert(ibert, f[1], nopart, kind, ep, alpha, beta, gam);//::: in g4bertini
    if (mat == 1) cobert += wt[no];
    cobb += wt[no];
  }
  if (nopart == -1) {
    if (mat == 1) npseub += wt[no];
    npsub += wt[no];
    if ((energy[no] != 1.0 || ec[no] != 1)) coupic += wt[no];
  }
  if (nopart == 0 && mat == 1) nozerb += wt[no];
  // if ced collision options cycle back to call bert (or bertcem) so as to avoid including
  // pseudo or non-productive collisions on tape
  switch (ifircl) {
  case 0:
    break;
  case 1:
    if ((nopart <= 0) && (name[no] == 1) && (mat == 1)) goto L429;
    break;
  default:
    if (nopart < 0) goto L429;
    break;
  }
  ibert = 1;
  if (nopart <= 0) {
    if (nopart < 0) {
      if (tip[no] != 4.0 || energy[no] != 1.0) {
	++pevts[mat];
	lelem = 0;
	nopart = -1;
	goto L5;
      }
      ++cou12;
      if (icem == 0) {
	cobert -= wt[no];
	if (mat == 1) npseub -= wt[no];
	npsub -= wt[no];
      } else {
	cobcem -= wt[no];
	if (mat == 1) npseub -= wt[no];
	npsub -= wt[no];
      }
      goto L12;
    }
    ++revts[mat];
    switch (ityp) {
    case 1:
    case 2:
      apr =  a[lelem][mat] + 1.0;
      zpr = zz[lelem][mat] + 1.0 - tip[no];
      ex = ec[no] + 7.0;
      recoil();
      // NOTE: nbogus is printed out as nexite in subroutine input
      if (nbogus > 0) ex -= erec;
      break;
    case 3:
    case 5:
      apr =  a[lelem][mat];
      zpr = zz[lelem][mat] + 3.0 - tip[no];
      ex = ec[no] + 139.9;
      recoil();
      if (nbogus > 0) ex -= erec;
      break;
    default:
      if (verboseLevel > 1) G4cout << " cascad: ityp " << ityp << G4endl;
      G4cerr << "cascad2" << G4endl;
    }
    return;
  }
  nhistm = maxcas * maxbch;
  G4int i;
  for (i = 0; i < nopart; ++i) {
    wtfas[i] = wt[no1];
    kindi[i] = kind[i] + 1;
    kindi[i] > 2 ? ibb = 0: ibb = 1;
    ibbarr[i] = ibb;
    amm = ginum[kindi[i]] * MeV;
    esum += amm + ep[i];
    kindi[i] <= 30 ? ic = ijevnt[ic] : ic = kindi[i];
  }
  ami = ginum[ityp] * MeV;
  einit = ec[no] + ami;
  edbrt = einit - esum;
  if (icem != 1 && iwr <= 20) {
    if (edbrt > 2.0) {
      if (nopart > 0) {
	esum2 = 0.0;
	for (G4int i = 1; i < nopart; ++i) {
	  amm = ginum[kindi[i]] * MeV;
	  esum2 = esum2 + amm + ep[i];
	  edbrt = einit - esum2;
	}
      }
    }
  }
  if (nopart < 0) {
    if (tip[no] != 4.0 || energy[no] != 1.0) {
      ++pevts[mat];
      lelem  = 0;
      nopart = -1;
      goto L5;
    }
    ++cou12;
    if (icem == 0) {
      cobert -= wt[no];
      if (mat == 1) npseub -= wt[no];
      npsub -= wt[no];
    } else {
      cobcem -= wt[no];
      if (mat == 1) npseub -= wt[no];
      npsub -= wt[no];
    }
    goto L12;
  } else if (nopart == 0) {
    ++revts[mat];
    switch (ityp) {
    case 1:
    case 2:
      apr =  a[lelem][mat] + 1.0;
      zpr = zz[lelem][mat] + 1.0 - tip[no];
      ex = ec[no] + 7.0;
      recoil();
      if (nbogus > 0) ex -= erec;
      break;
    case 3:
    case 5:
      apr =  a[lelem][mat];
      zpr = zz[lelem][mat] + 3.0 - tip[no];
      ex =  ec[no] + 139.9;
      recoil();
      if (nbogus > 0) ex -= erec;
      break;
    default:
      if (verboseLevel > 1) G4cout << " cascad: ityp " << ityp << G4endl;
      G4cerr << "cascad2" << G4endl;
    } 
    return;
  }
  pi0 = 0.0;
  ++revts[mat];
  sume = 0.0;
  prono = 0.0;
  pipos = 0.0;
  pineg = 0.0;
  G4int n;
  for (n = 0; n < nopart; ++n) {
    switch (kind[n]) {
    case 0:
      ++prono;
      break;
    case 1:
      break;
    case 2:
      ++pipos;
      break;
    case 3:
      ++pi0;
      break;
    case 4:
      ++pineg;
      break;
    default:
      if (verboseLevel > 1) G4cout << " cascad: kind[n]+1 " << kind[n] + 1 << G4endl;
      G4cerr << "cascad3" << G4endl;
    }
    sume += ep[n];
  }
  chgpis = pipos + pineg;
  fpt = nopart - chgpis - pi0;
  if (tip[no] <= 1.0) {
    apr =  a[lelem][mat] + 1.0 - fpt;
    zpr = zz[lelem][mat] + 1.0 - tip[no] - prono - pipos + pineg;
    if (zpr < 0.0) {
      ++negzpr;
      if (negzpr > 5) {
	G4cerr << "cascad4" << G4endl;
    }
      --revts[mat];
      goto L430;
    }
    ex = ec[no] + (1.0 - fpt) * 7.0 - sume - chgpis * 139.9 - pi0 * 135.1;
  } else {
    apr =  a[lelem][mat] - fpt;
    zpr = zz[lelem][mat] + 3.0 - tip[no] - prono - pipos + pineg;
    ia = static_cast<G4int>(apr + 0.5);
    iz = static_cast<G4int>(zpr + 0.5); // ::: verify
    nn = ia -iz;
    if (nn <= 0 || zpr < 0.0) {
      ++negzpr;
      if (negzpr > 5) {
	G4cerr << "cascad4" << G4endl;
    }
      --revts[mat];
      goto L430;
    }
    ex = ec[no] + (1.0 - chgpis) * 139.9 - sume - fpt * 7.0 - pi0 * 135.1;
  }
  recoil();
  if (nbogus > 0) ex -= erec;
  return;
  // high energy interaction
 L440:
  if (tip[no] <= 1) ehicut = ehin;
  else ehicut = ehipi;
  iqev = 0;
  weni = wt[no];
  if (icon == 2) ikin = kindo[no];
  if (icon <= 1) ikin = -1;
  if (ikin > 5) goto L4440;
  if (f[6] <= 4) eskal = eskale[static_cast<G4int>(f[6])+1];
  if (f[2] >= eskal || ihecc == 1) goto L4440;
  jtim = iclock();
 L439:
  iskl   = 1;
  esum   = 0.0;
  iqev   = 0;
  apr    = 0.0;
  zpr    = 0.0;
  edqev  = 0.0;
  edqevv = 0.0;
  edskl  = 0.0;
  edskll = 0.0;
  //  skale() // not used in this hetc version
  nskcol += wt[no];
  coss = nskcol;
  if (mat == 1) counsk += wt[no];
  if (nopart == -1 && mat == 1) npseus += wt[no];
  if (nopart == -1) npsus += wt[no];
  if (nopart == 0  && mat == 1) nozers += wt[no]; // if ce collision in medium 1
  switch (ifircl) {
  case 0:
    break;
  case 1:
    if (nopart <= 0 && name[no] == 1 && mat == 1) goto L439;
    break;
  default:
    if (nopart == -1) goto L439;
    break;
  }
  ktim = iclock();
  if (nopart > 0) {
    esumw = 0.0;
    esum = 0.0;
    for(G4int i = 0; i < nopart; ++i) {
      ik = kind[i] + 1;
      ic = ijevnt[ik];
      kindi[i] = ik;
      if (kindi[i] > 2) ibb = 0;
      else ibb = 1;
      ibbarr[i] = ibb;
      esum += (ep[i] + ginum[kindi[i]] * MeV) * wtfas[i];
      esumw += ep[i] + ginum[kindi[i]] * MeV;
      nhistm = maxcas * maxbch;
      if (wtfas[i] < 0.9) ++nskwt;
    }
    einit  = ec[no] + ginum[ityp] * MeV;
    edskl  = einit - esum;
    edskll = einit - esumw;
  }
  goto L4441;
 L4440:
  j      = 21;
  iqev   = 1;
  iskl   = 0;
  apr    = 0.0;
  zpr    = 0.0;
  edskl  = 0.0;
  edskll = 0.0;
  weni   = wt[no];
  ktim   = iclock();
  zpr    = zress;
  apr    = aress;
  nn     = static_cast<G4int>((apr - zpr) + 0.5);
  if (nopart > 0) {
    einit = ec[no] + ginum[ityp] * MeV;
    esum = 0.0;
    esumw = 0.0;
    for (G4int i = 0; i < nopart; ++i) {
      kindi[i] <= 30 ? ic = ijevnt[kindi[i]] : ic = kindi[i];
      // kindi is defined in qevent to include kaons, etc
      // != kind[i] + 1 if icon == 2 || icon == 0
      esumw += ginum[kindi[i]] * MeV + ep[i];
      esum += ginum[kind[i]+1] * MeV + ep[i];
      // note that if icon = 1 or icon = 2 particles have been converted to original 
      // seven types, conserving energy, so esum is the quantit to use wtfas[i] contains weni
      sumsi = sqr(alpha[i]) + sqr(beta[i]) + sqr(gam[i]); //::: replace pow -> square etc.
      icp = static_cast<G4int>(kind[i] + 1.01);
      if (icp <= 30) icp = ijevnt[icp];
    }
    edqev = einit - esum;
    edqevv = einit - esumw;
  }
 L4441:
  ic = ijevnt[ityp];
  if (nopart < 0)  { // pseudo event
    ++pevts[mat];
    ++hpevt[mat];
    lelem = 0;
    nopart = -1;
    goto L5;
  }
  // real event
  ++revts[mat];
  ++hrevt[mat];
  if ((iqev == 1 && nopart <= 0 && nofask <= 0) || (iqev == 0 && nopart <= 0)) {
    switch (ityp) {
    case 1:
    case 2:
      apr = a[lelem][mat] + 1.0;
      zpr = zz[lelem][mat] + 1.0 - tip[no];
      ex = ec[no] + 7.0;
      recoil();
      if (nbogus > 0) ex -= erec;
      break;
    case 3:
    case 5:
      apr = a[lelem][mat];
      zpr = zz[lelem][mat] + 3.0 - tip[no];
      ex = ec[no] + 139.9;
      recoil();
      if (nbogus > 0) ex -= erec;
      break;
    default:
      if (verboseLevel > 1) {
	G4cout << " cascad: ityp " << ityp << G4endl;
      }
      G4cerr << "cascad2" << G4endl;
    
    return;
  }
  // productive event
  for (i = 0; i < 5; ++it) prd[i] = 0.0;
  if (nopart > 0) {
    for (G4int n = 0; n < nopart; ++n) {
      it = kind[n] + 1;
      if (it <= 5) prd[it] += wtfas[n];
    }
  }
  // compute residual z
  if (iqev != 1) {
    tip[no] <= 1.0 ? zdum = 1.0 : zdum = 3.0;
    zzz = zz[lelem][mat] + zdum - tip[no] - prd[1] - prd[3] + prd[5];
    izz = static_cast<G4int>(zzz);
    yzz = static_cast<G4double>(izz);
    xzz = zzz - yzz;
    xzz < 0.5 ? zpr = static_cast<G4double>(izz) : zpr = static_cast<G4double>(izz + 1);
    // computes residual a
    if (zpr < 0.0) {
      ++negzpr;
      if (negzpr > 5) {
	G4cerr << "cascad5" << G4endl;
}
      --revts[mat];
      --hrevt[mat];
      goto L440;
    }
    aaa = rmfas / (ginum[1] * 500.0 + ginum[2] * 500.0);
    iaa = static_cast<G4int>(aaa);
    yaa = static_cast<G4double>(iaa);
    xaa = aaa - yaa;
    xaa >= 0.5 ? apr = static_cast<G4double>(iaa) : apr = static_cast<G4double>(iaa + 1);
  }
  ex = exfas;
  erec = refas;
  mtim = iclock();
  ittim = static_cast<G4int>((static_cast<G4double>(mtim - itim) / 100.0));
  // attempt to check energy conservation
  // einit + mass of target = esum + ex + erec + mass of compound nucleus
  // amasw = energyAZ(anucc, znucc) + conver * MeV * anucc in MeV	
  // amscom = energyAZ(apr, zpr) + conver * MeV * apr
  // conver = 1.0
  // convert to MeV
  convrr = conver * MeV;
  nn = static_cast<G4int>(apr - zpr);
  // nn <= 0.0 ? amscom = rmfas : amscom = energyAZ(apr, zpr) + convrr * apr; //::: fix
  if (ibrt == 1) rmfas = amscom;
  amasdf = rmfas - amscom;
  wapsen = energyAZ(anucc, znucc);
  amasw = wapsen + convrr * anucc;
  amasd = amasw - rmfas - ex - erec;
  eloss = elab + amasw - (esum + amscom + ex + erec);
  elosss = einit + amasd - esum;
  // above leaves out mass and ke. of particles included in nofask
  entr = 0.0;
  if (nofask > 0) {
    for (G4int j = 0; j < nofask; ++j) entr += epk[j] + ginum[kindk[j] + 1] * MeV;
  }
  eloss -= entr;
  elosss -= entr;
  if (iskl == 1) elos = edskl + amasd - entr;
  if (iqev == 1) elos = edqev + amasd - entr;
  if (ibrt == 1) elos = edbrt + amasd - entr;
  adiff = anucc - apr;
  zdiff = znucc - zpr;
  return;
  }
}

//---------------------------------------------------------------------------------------
// methods ready for testing

void G4BertiniCascade::update() {
  // decription: transfers particles produced above cut-off to the transport bank
  //             the first (highest energy) produced particle replacing the current
  //             transport particle
  // parameters:
  // uses:
  // changes:
  static G4bool first = true;
  if (first) {
    nomaxp = 0;
    first = false;
  }
  if (no > MAXPART) G4Exception("update: bank size exeeded");
  noo = no;
  io = iou;
  if (nabov <= 0) { // if no particles were were produced above cut-off
    energy[no] = -1.0; // the current particle is changed to a dead particle
    return;
  }
  wtaa = wta[1];
  tip[no] = tipa[1];
  name[no] = namea[1];
  energy[no] = ea[1];
  u [no] = ua[1];
  v [no] = va[1];
  w [no] = wa[1];
  wt[no] = wta[1];
  x[no] = xc[no];
  y[no] = yc[no];
  z[no] = zc[no];
  i151[no] = 0;
  kindo[no] = kinda[1];
  barr[no] = ibarra[1];
  ic = kindo[no];
  if (ic <= 30) ic = ijevnt[ic];
  ic = kindo[no];
  if (ic <= 30) ic = ijevnt[ic];
  if (nabov - 1 <= 0) return;
  mxcess = 0;
  for (G4int ll = 1; ll < nabov; ++ll) {
    if (wta[ll] <= wtlow) continue;
    m = nomax + ll - 1;
    if (m > MAXPART) {
      ++mxcess;
      mobch = nobch + 1;
      if (verboseLevel > 1) {
	G4cout 
	  << m << " particles produced above cutoff in cascade " << nocas 
	  << " of batch " << mobch
	  << G4endl;
      }
      continue;
    }
    name[m] = namea[ll];
    tip[m] = tipa[ll];
    energy[m] = ea[ll];
    u[m] = ua[ll];
    v[m] = va[ll];
    w[m] = wa[ll];
    wt[m] = wta[ll];
    x[m] = xc[no];
    y[m] = yc[no];
    z[m] = zc[no];
    blz[m] = blz[no];
    nmed[m] = nmed[no];
    icp = static_cast<G4int>(tip[m] + 1.01);
    if (icp <= 30) icp = ijevnt[icp];
    i151[m] = 0;
    barr[m] = ibarra[ll];
    kindo[m] = kinda[ll];
    ic = kindo[m];
    if (kindo[m] <= 30) ic = ijevnt[kindo[m]];
  }
  // here, particle is scratched
  nomax += nabov - 1;
  ktim = iclock();
  jtim = (ktim - itim) / 100;
  if (nomax > nomaxp + 100) nomaxp += 100;
  nomaxn = nomax;
  nomax -= mxcess;
  return;
}

void G4BertiniCascade::range() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double barik[3];
  G4double pkm[8];
  G4double sum1;
  G4double sum2;
  G4double sum3;
  G4double f4;
  G4double ek[8] = {1., 2., 5., 10., 50., 100., 200., 400.};
  // pk has been flipped
  G4double pk[3][8] = { {1.37, 1.24, 1.12, 1.07, 1.02, 1.0, 1.0, 1.0},
			{1.47, 1.33, 1.27, 1.19, 1.06, 1.04, 1.00, 1.0},
			{1.60, 1.45, 1.31, 1.22, 1.19, 1.16, 1.11, 1.0} };
  G4int jtyp[3] = {1, 3, 6};
  elow = 1.0;
  G4double elowd=1.0;
  G4double emx = emax * 10.0;
  G4int npts = 900;
  G4int npts5 = npts * 5;
  einc = log(emx / elowd) / npts5;
  ++npts5;
  barik[1] = log(zfoi(4.0) * 1.0e-06); //::: MeV<> eV ?
  barik[2] = log(zfoi(13.0) * 1.0e-06);
  barik[3] = log(zfoi(82.0) * 1.0e-06);
  for (G4int n = 0; n < mxmat; ++n) {
    sum1 = denh[n];
    sum2 = denh[n] * log(zfoi(1.0) * 1.0e-06);
    sum3 = denh[n];
    nnn = nel[n];
    G4double t;
    for (G4int i = 0; i < nnn; ++i) {
      t = den[i][n] * zz[i][n];
      sum3 =+ den[i][n];
      sum1 =+ t;
      sum2 =+ t * log(eion[i][n]);
    }
    G4double constm = 0.26057 * sum1;
    G4double bari = sum2 / sum1;
    if (bari <= barik[1]) {
      ik = 1;
      for (G4int i = 0; i < 8; ++i) pkm[i] = pk[ik][i];
    }
    else if (bari <= barik[2]) {
      ik = 2;
      fact = (bari - barik[ik - 1]) / (barik[ik] - barik[ik - 1]);
      for (G4int i = 0; i < 8; ++i)
	pkm[i] = pk[ik - 1][i] + fact * (pk[ik][i] - pk[ik - 1][i]);
    }
    else if (bari >= barik[3]) {
      ik = 3;
      for (G4int i = 0; i < 8; ++i) pkm[i] = pk[ik][i];
    }
    G4double rng1 = 0.0;
    G4double sig1 = 0.0;
    G4int j  = 0;
    G4int k = 0;
    G4int ie = 2;
    G4double pke = 0.0;
    G4double f1(0.0);
    G4double f2;
    G4double f3(0.0); //::: default value added
    for (i = 0; i < npts5; ++i) {  
      e2 = elowd * exp((i - 1) * einc); 
      f2 = dxde(e2, sum1, sum2);
      G4double estar = e2 / 938.232;
      G4double temp1 = (estar + 1.0) * (estar + 1.0);
      G4double bsq = (temp1 - 1.0) / temp1;
      G4double g2 = 1.0 / sqrt(1.0 - bsq);
      G4bool pkeTest = true;
      if (pke != 1.0) {
	//::: fixdo {   
	if (e2 < ek[ie]) {
	  pke = pkm[ie - 1] + (e2 - ek[ie - 1]) / (ek[ie] - ek[ie - 1]) * 
	    (pkm[ie] - pkm[ie - 1]);
	  pkeTest = false;
	  break;
	}
	++ie;
	//:: fix   } 
	while (ie <= 8);
	if (pkeTest) pke = 1.0;
      }
      f4 = (1.0 - 0.5 * bsq) * pke / 
	((1.0 - bsq) * (1.0 + 0.0010893 * g2)) * pow(f2, 3);
      if (i > 0) {
	rng1 += (e2 * f2 + e1 * f1) * einc / 2.0;
	sig1 += (e2 * f4 + e1 * f3) * einc / 2.0;
	++k;
	if (k == 5) {
	  ++j;
	  rnge[j][n] = rng1;
	  sigr2[j][n] = sig1 * constm;
	  k = 0;
	}
      }
      f1 = f2;
      e1 = e2;
      f3 = f4;
    }
  }
  einc *= 5;
  for (G4int i = 0; i < npts; ++i) erng[i] = elowd * exp(i * einc);
  G4int j1 = 0;
  if (emin[1] == elow) j1 = 1;
  for (n = 0; n < mxmat; ++n) {
    for (G4int j = 0; j < 7; ++j) {
      rngelo[j][n] = 0.0;
      sigelo[j][n] = 0.0;
    }
    G4int itp;
    for (j = j1; j < 3; ++j) { 
      itp = jtyp[j];
      rainge(emin[itp], n, itp);
      rngelo[itp][n] = rbar;
      sigelo[itp][n] = sigr;
    }
    rngelo[5][n] = rngelo[3][n];
    rngelo[7][n] = rngelo[6][n];
    sigelo[5][n] = sigelo[3][n];
    sigelo[7][n] = sigelo[6][n];
  }
  return;
}

void G4BertiniCascade::gene(G4double *z) {
  // dimension z(101);
  G4double cd = com * 1.0e2;
  G4int i = static_cast<G4int>(cd + 1.0);
  G4double az = z[i];
  G4double bz;
  G4double cz;
  G4double xz;
  G4double sca;
  G4double disc; 
  if (i == 0) {
    cz = z[i + 1];
    sca = (cz - az);
    bz = az + sca * 0.7071067812;
    xz = cd + cd;
  } else {
    bz = z[i + 1];
    if (i + 1 > 100) {
      G4Exception("gene1-22");
      return;
    }
    else if (i + 1 == 100) {
      cz = bz + 0.50 * (bz - az);
    } else {
      cz = z[i + 2];
    }
    xz = cd - (i - 1);
  }
  G4double sba = bz - az;
  G4double sqa = az * az;
  G4double sqac = sqa - cz * cz;
  G4double sqba = bz * bz - sqa;
  G4double rb = sqac + sqba + sqba;
  G4double rc = az * cz * sca - sba * (2.0 * az * bz + xz * (bz - cz) * sca);
  G4double ra = sca - sba - sba;
  if (ra == 0.0) {
    com = az + xz * sba;
    return;
  }
  disc = sqr(rb) - 4.0 * ra * rc;
  if (disc < 0.0) {
    G4Exception("gene1-22");
    return;
  }
  disc = sqrt(disc);
  G4double plus = (disc - rb)/(2* ra);
  G4double aminus = (-disc -rb)/(2* ra);
  if (i == 1) {
    bz = cz;
    xz -= cd;
    sba = cz-az;
  }
  if (plus > bz || plus < az) {
    if (aminus > bz || aminus < az) {
      G4Exception("gene1-22");
      return;
    }
    com = aminus;
    return;
  }
  if (aminus > bz) {
    com = plus;
    return;
  }
  if (aminus >= az) {
    ra = xz * sba + az;
    rb = fabs(ra - aminus);
    rc = fabs(ra - plus);
    if (rb > rc){
      com = plus;
    } else {
      com = aminus;
    }
    return;
  }
  com = plus;
  return;
}

void G4BertiniCascade::relas(G4double minEnergy) { //::: verify hole routine
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int npoits = 500;
  G4int idt[4];
  G4double hol[10];
  // ::: if (nelstp != in) rewind nelstp;
  locsig[1] = 1;
  locf1[1] = 1;
  G4int loc1;
  G4int loc2;
  G4int isave;
  G4int isavef;
  G4int k;
  for (k = 0; k < noelas; ++k) {
    //:::read(nelstp,10,end=999)(idt(i),i=1,4),fmas,(hol(i),i=1,10);
    //:::10 format(4i5,f10.4,10a3);//:::
    G4double fmas;
    if (verboseLevel > 1) {
      G4cout 
	<< " and endfile has been encountered in subroutine relas while reading nelstp " << nelstp
	<< " the (n, N) elasic data "
	<< G4endl;
    }
    G4Exception("relas1");
    if (k != 1) locsig[k] = locsig[k - 1] + isave;
    loc1 = locsig[k];
    loc2 = loc1 + idt[4] - 1;
    isave = idt[4];
    if (verboseLevel > 1) {
      G4cout << " idt:" ;
      for (G4int i = 0; i < 4; ++i) G4cout << " " << idt[i];
      G4cout << " fmas " << fmas << G4endl;
      G4cout << " hol:";
      for (i = 0; i < 10; ++i) G4cout << " " << hol[i];
      G4cout << G4endl;
    }
    if (locsig[k] - 1 + idt[4] > npoits) {
      if (verboseLevel > 1) {
	G4cout 
	  << " total room for elastic (n, N) data exceeded "
	  << " dimension on commons ef1, esge, f1s, totsge " << npoits
	  << " it must be increased and data statement for npoits should be changed"
	  << G4endl;
      }
      G4Exception("relas1");
    }
    // read(nelstp,12,end=999)(es(i),sige(i),i=loc1,loc2);
    // 12 format(2e15.5);
    // read(nelstp,14);
    // 14 format(1h);
    for (G4int i = loc1 - 1; i < loc2; ++i) es1[i] *= 1.000001e-06; //::: indexes?
    //ck(k, idt[4], &es1, elas, minEnergy); //::: fix
    if (nledit != 0) {
      if (verboseLevel > 0) {
	G4cout << "es1[i], sige[i] " << G4endl;
	for (G4int i = loc1 - 1; i < loc2; ++i) G4cout << energy[i] << sige[i]  << G4endl; 
      }
    }
  }
  // read(nelstp,10,end=999)(idt(i),i=1,4),fmas,(hol(i),i=1,10); //:::
  if (k != 1) locf1[k]= locf1[k - 1] + isavef;
  loc1   = locf1[k];
  loc2   = loc1 + idt[4] - 1;
  isavef = idt[4];
  locsig[noelas + 1] = locsig[noelas] + isave;
  locf1[noelas + 1] =   locf1[noelas] + isavef;
  // ::: if (nelstp != in) rewind nelstp;
  return;
}
void G4BertiniCascade::dklos() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // dklos calculates the fraction dkwt and multiplies wt[no] by dkwt
  G4int mt;
  G4int ikinpr;
  static G4bool first = true;
  G4double velc(2.997925);
  G4int k;
  if (first) {
    first = false;
    ikinpr = 0;
    velc *= 1.0e+10; // velc is velocity of light
    for (G4int k = 1; k < 30; ++k) {
      ic = ijevnt[k];
      tauinv[ic] = 1.0;
      if (tau[ic] > 0.0) tauinv[ic] = 1.0 / (velc * tau[ic]);
    }
    for (k = 30; k < 110; ++k) {
      ic = k;
      tauinv[k] = 1.0;
      if (tau[ic] > 0.0) tauinv[k] = 1.0 / (velc * tau[ic]);
    }
    ikind = 0;
    ineg = 0;
  }
  i310 = 0;
  i300 = 0;
  if (mat == 6666)
    mt = mxmat + 1;
  else
    mt = mat;
  igo = 0;
  if (kindo[no] != ikinpr) {
    if (kindo[no] > 7 || (kindo[no] == 4 && tip[no] != 3)) {
      ++ikind;
      ic = kindo[no];
      if (ic < 31) ic = ijevnt[ic];
      ikinpr = kindo[no];
    }
  }
  ineg = 0;
  igo = 2;
  dkwt = 1.0;
  k = static_cast<G4int>(kindo[no] + 0.001);
  ic = k;
  if (ic <= 30) ic = ijevnt[k];
  G4int ikt = static_cast<G4int>(tip[no] + 1.001);
  ityp = ikt;
  amm = aam[ic] * MeV;
  sigdk = tauinv[ic] / sqrt((ec[no] / amm + 1.0) * (ec[no] / amm + 1.0) - 1.0);
  if (ic == 1 || ic == 8 || ic == 9 || ic == 2) {
    // real nucleons and antinucleons. no decay
    dkwt = 1.0;
    delsig = sigmx[ityp][mt] + totels + bodge; // add magnet bodge to delsig
    return;
  }
  if (tauinv[ic] == 1.0) {
    // tau[ic] = 0.0; tauinv[ic] was set to 1.0. or,
    // sigmas and lambdas are sent here,
    // with small tau[ic]
    dkwt = 1.0e-08;
    ++i300;
    // since wt[no] = wt[no] * dkwt, and dklos is called in cascad, the factor
    // dkwt must be considered the probability a collision (either real or pseudo)
    // will occur in cascad:pcol, or hcobert, skale, or qevent
    // note this practically forces a decay of the no particle in mfpd2
    // (with m = 2 or m = 5) rather than a collision in cascad
    // so in mfpd2, weight=oldwt * (1.0 - dkwt) when = oldwt = original wt[no]
    //
    // NOTE: that in revised decay, kaons with zero lifetime are decayed immediately
    // also, many resonances are decayed immediately
    // also, at the moment, k-long is decayed immediately, when decay called
    // however, k-longs can be formed in nucrin, mfpd2, and hco
    // and are not yet decayed immediately but converted to neutrons
    weig   = oldwt - dkwt * wt[no];
    wto    = wt[no];
    wt[no] = wt[no] * dkwt;
    weight = oldwt - wt[no];
    if (wt[no] <= wtlow) {
      wtno = wt[no];
      ++nwtc;
    }
    return;
  }
  delsig = sigmx[ityp][mt] - sigdk + bodge;
  dkwt   = delsig / (sigmx[ityp][mt] + bodge);
  dkwto  = dkwt;
  if (tauinv[ic] < 1.0e-28) {
    dkwt = 0.9999;
    ++i310;
    // this value of dkwt is arbitrary. it hardly changes the weight of the particle
    // in mfpd2, weight = oldwt - wt[no] is close to 0.0
    // for the particles resulting from decay
    // here, tau[ic] is very large; tauinv[ic] is very small
    weig   = -dkwt * wt[no] + oldwt;
    wto    =         wt[no];
    wt[no] =  dkwt * wt[no];
    weight =        -wt[no] + oldwt;
    if (wt[no] <= wtlow) {
      wtno = wt[no];
      ++nwtc;
    }
    return;
  }
  if (dkwt <= 0.0) {
    ++ineg;
    if (kindo[no] > 15 && kindo[no] <= 20) {
      dkwt = 1.0e-08;
      ++i300;
      weig   = - dkwt * wt[no] + oldwt;
      wto    =          wt[no];
      wt[no] =   dkwt * wt[no];
      weight =        - wt[no] + oldwt;
      if (wt[no] <= wtlow) {
	wtno = wt[no];
	++nwtc;
      }
      return;
    }
    if (ineg == 1) {
      switch(ityp) {
      case 1:
      case 2:
	dkwt =1.0;
	delsig = sigmx[ityp][mt] + totels + bodge; // add magnet bodge to delsig
	return;
      case 3:
      case 5:
	sigdk = 0.001259 / sqrt((ec[no] / 139.9 + 1.0) * (ec[no] / 139.9 + 1.0) - 1.0);
	sigdki = sigdk;
	break;
      case 4:
	G4Exception("dklos1");
      case 6:
      case 7:
	sigdk = 1.587e-05 / sqrt((ec[no] / 107.0+ 1.0) * (ec[no] / 107.0 + 1.0) - 1.0);
	sigdki = sigdk;
      }
      // add magnet bodge to delsig and allow for it in dkwt
      delsig = sigmx[ityp][mt] - sigdk + bodge;
      dkwt = delsig / (sigmx[ityp][mt] + bodge); // thus, 1.0 - dkwt = sigdk / (sigmx[ityp][mt] + bodge this is the probablity factor for a decay
      if (dkwt <= 0.0) {
	++ineg;
	if (kindo[no] > 15 && kindo[no] <= 20) {
	  dkwt = 1.0e-08;
	  ++i300;
	  weig   = - dkwt * wt[no] + oldwt;
	  wto    =          wt[no];
	  wt[no] =   dkwt * wt[no];
	  weight =        - wt[no] + oldwt;
	  if (wt[no] <= wtlow) {
	    wtno = wt[no];
	    ++nwtc;
	  }
	  return;
	}
      }
      wtnew  = dkwt * wt[no];
      wto    =        wt[no];
      wt[no] = dkwt * wt[no];
      weight =      - wt[no] + oldwt;
      if (wt[no] <= wtlow) {
	wtno = wt[no];
	++nwtc;
      }
      return;
    }
  }
  wtnew  = dkwt * wt[no];
  wto    =        wt[no];
  wt[no] = dkwt * wt[no];
  weight =      - wt[no] + oldwt;
  if (wt[no] <= wtlow) {
    wtno = wt[no];
    ++nwtc;
  }
  return;
}

void G4BertiniCascade::go(){
  // description : start testing
  // parameters  : -
  // uses        :
  // changes     :
  mainBodyInit();
}

void G4BertiniCascade::shxd() {
  // description : initialize data matrices locx and etc
  // parameters  : -
  // uses        :
  // changes     :
  if (verboseLevel > 1) G4cout << "shxd: initializing data vector index ..." << G4endl;
  de = 20.0;
  locx[0][0] = 995;
  locx[1][0] = 1153;
  locx[2][0] = 3793;
  locx[3][0] = 0;
  locx[0][1] = 1283;
  locx[1][1] = 1441;
  locx[2][1] = 3969;
  locx[3][1] = 0;
  locx[0][2] = 2009;
  locx[1][2] = 0;
  locx[2][2] = 3667;
  locx[3][2] = 0;
  locx[0][3] = 2243;
  locx[1][3] = 0;
  locx[2][3] = 3541;
  locx[3][3] = 3415;
  for (G4int it = 0; it < 4; ++it) {
    for (G4int id = 0; id < 4; ++id) eth[id][it] = 0.0;
  }
  for (it = 1; it < 2; ++it) {
    eth[0][it] = 360.0;
    eth[1][it] = 920.0;
  }
  for (it = 2; it < 4; ++it) eth[0][it] = 180.0;

  if (verboseLevel > 2) {
    G4cout << "shxd: eth initialized to :" << G4endl;
    for (G4int it = 0; it < 4; ++it) {
      for (G4int id = 0; id < 4; ++id) cout << "\t" <<eth[id][it];
      cout << G4endl;
    }
  }
}

G4double G4BertiniCascade::xlamb(const G4double x, 
			  const G4double y, 
			  const G4double z) {
  // description : initialize data matrices locx and etc
  // parameters  : -
  // uses        :
  // changes     :
  return sqrt(fabs(sqr(x) - 2.0 * x * (y + z) + sqr(y - z)));
} 


G4bool G4BertiniCascade::main3() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  if (no < nomax) {
    ++no;
    noo = no;
    ic = kindo[no];
    if (ic <= 30)ic = ijevnt[kindo[no]];
    return true;
  }
  if (nocas < maxcas) {
    ++nocas;
    mainBody();
  }
  return false;
}


void G4BertiniCascade::main4() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // increase narriv, call mfpd2 to decay the no particle
  // call update to put nabov particles into COMON.F
  // call analz1 to write both nabov and nbelo arrays on tape
  // nabov = 0, so do not need to call update
  // however, dklos was called if this collision, so why not call mfpd2(2)
  // to decay the no particle with proper probability?
  // actually, neutrons and protons are not decayed, nor are antinucleons
  ++narriv;
  if (ibertp < 0) mfpd2(2);
  analz1(2);
  if (ifircl <= 1) update();
  nabov = 0;
  nbelo = 0;
  if (nomax > MAXPART) {
    if (nocas >= maxcas) mainFinal();
    ++nocas;
    mainBody();
  }
  mainBody2();
}

void G4BertiniCascade::main5() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  mat = nmed[no];
  if (tip[no] == 1.0) {
    nbeta = 1;
    ec[no] = energy[no];
  } else {
    rainge(energy[no], mat, ityp);  
    rbtors = rbar / rng;
  }
  getflt();
  dd = d;
  oldwt = wt[no];
  if (ityp == 2) mainBody3b();
  rr = rng - d;
  if (rr > 0.0) {
    nbeta = 1;
    rp = rr * rbtors;
    ecol(rp, mat, ityp, ec[no], energy[no]);//:::
    if (ec[no] >= elop || ityp != 1) mainBody3a();
  }
  mainBody3();
}

void G4BertiniCascade::mainBody3() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  nbeta = 2;
  d = rng;
  ec[no]  = emin[ityp];
  // particles at rest do goto gomprp
  mainBody3a();
}

void G4BertiniCascade::mainBody4() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  x[no]  = xc[no];
  y[no]  = yc[no];
  z[no]  = zc[no];
  energy[no]  = ec[no];
  if (ityp == 2) main5();
  if (ec[no] == emin[ityp]) {
    nbeta = 2;
    ngam = 3;
    ncol = 3;
    m = 3;
    nopart = 0;
    icno = kindo[no];
    if (kindo[no] <= 30) icno = ijevnt[kindo[no]]; 
    if (ibertp < 0) mfpd2(3); // note: call mfpd2(3) if nbeta > 1, or if ec[no] == emin[ityp]
    mmf = 3;
    main6();
  }
  if (ncol == 7) {
    if (mat == 6666 || medeq[mat] != nmed[no]) main5();
  }
  rng -= d;
  getflt();
  dd = d;
  oldwt = wt[no];
  if (ityp == 2) mainBody3b();
  rr = rng - d;
  if (rr <= 0.0) mainBody3();
  nbeta = 1;
  rp = rr * rbtors;
  ecol(rp, mat, ityp, ec[no], energy[no]);//:::
  if (ec[no] < elop && ityp == 1) mainBody3();
}


void G4BertiniCascade::p1cli() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // momentum x and y coordinates, particle 1 = 0.0 
  // z coordinate = sqrt(sqr(total energy)-sqr(mass)) 
  // for particle one.  p1oe1 = current(moment/total)
  // USES: energy[0], massParticle[0], pxyz[8]
  // MODIFIES: pxyz[0,4,8], p1oe1 
  pxyz[0] = 0.0;
  pxyz[4] = 0.0;
  pxyz[8] = sqrt(sqr(energy[0]) - sqr(massParticle[0]));
  p1oe1   = pxyz[8] / energy[0];
}

G4double G4BertiniCascade::xsec(G4int id, 
			 G4int itp,  
			 G4double eCurrent) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // interpolate cross section data for the energy eCurrent
  G4int i, lc;
  if (ifirst == 0) ifirst = 1;
  G4double  et = eth[id][itp];
  lc = locx[id][itp];
  i = static_cast<G4int>((eCurrent - et) / de + 1.0);
  G4double e = static_cast<G4double>(i - 1) * de + et;
  G4double temp = static_cast<G4double>(hadronCrossSection[i + lc - 1]);
  G4double xsec = temp + (eCurrent - e) / de * (static_cast<G4double>(hadronCrossSection[i + 1 + lc - 1]) - temp);
  return xsec;
}

void G4BertiniCascade::rou16(G4double *t) {
  // decription: 
  // parameters:
  // uses: i3
  // changes: i2, cst, i3
  // snt = mud(); /:::
  //  crdet(51, t[1], rlke);
  i2 = 0;
  // cst = crdt[1] - fabs(snt * (crdt[1] - crdt[0])); //::: uncomment crdet
  switch (i3) {
  case 1: i3 = -1;
  case 2: i3 = 0;;
  case 3: i3 = 1;;
  }
}

void G4BertiniCascade::modify(G4double &u,
		       G4double &v,
		       G4double &w,
		       const G4double energy, 	
		       const G4double d, 	
		       const G4double particle) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int particleType = static_cast<G4int>(particle + 1.01);
  if (particleType == eNeutron) return;
  G4double p = energy * sqrt(1.0 + 2.0 * particleMass[eNeutron] / energy);
  G4double uvw = sqrt(sqr(u) + sqr(v) + sqr(w));
  if (fabs(uvw - 1.0) >= 0.000001) {
    u /= uvw;
    v /= uvw;
    w /= uvw;
  }
  fact = 0.2997 * particleCharge[particleType] * d / p;
  u += (v * bfield[3] - w * bfield[2]) * fact;
  v += (w * bfield[1] - u * bfield[3]) * fact;
  w += (u * bfield[2] - v * bfield[1]) * fact;
  uvw = sqrt(sqr(u) + sqr(v) + sqr(w));
  u /= uvw;
  v /= uvw;
  w /= uvw;
}

G4double G4BertiniCascade::enrg(const G4double a, 
			 const G4double z) {
  // decription: calculate eneergy ::: for atom
  // parameters: a atomic weight, z number of protons
  // uses: cam2, cam3
  // changes:
  // note: this function is taken from hetc.kfa
  // compare with deltas function in hetc. ornl function energy
  // a = n + z, atomic weight = neutrons + protons 
  G4double am2zoa = sqr((a - 2 * z) / a);
  G4double a13    = pow(a, oneThird);
  G4double am13   = 1.0 / a13;
  G4double ev = -17.0354 * (1.0 - 1.84619  * am2zoa) * a;
  G4double es =  25.8357 * (1.0 - 1.712185 * am2zoa) * (1.0 - 0.62025 * sqr(am13)) * (sqr(a13) - 0.62025);
  G4double ec = 0.799 * z * (z - 1.0) * am13 * (((1.5772 * am13 + 1.2273) * am13 - 1.5849) * sqr(am13) + 1.0);
  G4double eex = -0.4323 * am13 * pow(z, fourThirds) * (((0.49597 * am13 - 0.14518) * am13 - 0.57811) * am13 + 1.0);
  return 8.367 * a - 0.783 * z + ev + es + ec + eex + cam2[static_cast<G4int>(z)] + cam3[static_cast<G4int>(a - z)];
}

G4double G4BertiniCascade::energyAZ(const G4double a, 
			     const G4double z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // a = n + z, atomic weight = neutrons +protons 
  // this version of function energy is from hetc.kfa
  if (static_cast<G4int>(a - z) <= 0) G4Exception("energy-1");
  G4int ia    = static_cast<G4int>(a);
  G4double zz = waps[ia][0];
  G4int iz    = static_cast<G4int>(z - zz + 2.0);
  G4double energy;
  if (zz <= z && z <= (zz + 9) && waps[ia][iz] != 0.0) energy = waps[ia][iz];
  else energy = enrg(a, z);
  return energy;
}

G4double G4BertiniCascade::zfoi(const G4double z) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  const G4double a1(9.76);
  const G4double a2(58.8);
  const G4double xp(-0.19);
  G4double zr;
  if (z > 12.0) zr = a1 * z + a2 * pow(z, xp); 
  else zr = fli[static_cast<G4int>(z)]; 
  return zr;
}

void G4BertiniCascade::getrig() {
  // decription : calculate trigonometric function of polar and azimuthal angles
  // parameters : -
  // changes    : sinth, costh, cosphi, sinphi
  // uses       : direction cosines u, v, w
  G4double rt = sqr(u[no]) + sqr(v[no]);
  if (rt <= 0.0) {
    sinth  = 0.0;
    costh  = 1.0;
    cosphi = 1.0;
    sinphi = 0.0;
  } else {
    rt = sqrt(rt);
    sinth = rt;
    cosphi = u[no] / rt;     
    costh  = w[no];
    sinphi = v[no] / rt;
  }
}

void G4BertiniCascade::gtiso(G4double &u, 
		      G4double &v, 
		      G4double &w) {
  // decription: produce the direction cosines of an isotropic unit vector
  // parameters:
  // uses:
  // changes: u, v, w
  G4double x; 
  G4double y;
  G4double z(0.0);
  G4double xsq;
  G4double ysq;
  G4double zsq;
  G4double d(1.0);
  while (sqr(d) > z) {
    z = G4UniformRand();
    x = 0.687368 * frnd();
    y = 0.687368 * frnd();
    xsq = sqr(x);
    ysq = sqr(y);
    zsq = sqr(z);
    d = xsq + ysq + zsq;
  }
  u = 2.0 * x * z / d;
  v = 2.0 * y * z / d;
  w = (zsq - xsq - ysq) / d;
  return;
}

G4double G4BertiniCascade::gaurn() {
  // decription:   // gaussian distributed random variable  // ::: move to utils 
  // param`<eters:
  // uses`:
  // changes:

  // uses: exprnf()
  G4double test(1.0);
  G4double z(0.0);
  G4double ans;
  while (test > z) {
    ans = exprnf();
    z = exprnf();
    test = sqr(static_cast<G4int>(ans - 1.0)) / 2.0;
  }
  G4double r1 = 2.0 * G4UniformRand() - 1.0;
  if (r1 < 0.0) ans = -ans;
  return ans;
}

G4double G4BertiniCascade::exprnf() {
  // decription: exponentially distributed random number  // ::: move to utils 
  // parameters:
  // uses: 
  // changes:
  G4double randExp;
  G4double r;
  G4double ex;
  G4double whole(0.0);

  for(;;){
    randExp   = G4UniformRand();
    ex = randExp;

    for(;;) {
      r       = G4UniformRand();

      if (r  >= randExp) return ex + whole; 

      randExp = G4UniformRand();
      if (randExp > r) break;
    } 

    whole  += 1.0;
  }
}
	
G4double G4BertiniCascade::dxde(G4double ee, G4double sum1, G4double sum2) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // includes density effect. returned units are cm / MeV
  // note: k1 = log(2.0 * electronRestMass) (MeV)
  // k2 = 4 * pi * pow(4.803e-10 (esu)), 4) * 1.0e24 / ((9.109e-28 (gm/electron)) *
  // sqr(3.0e10 (cm/sec))) * (1.602e-06 (erg/mev))), final units-(MeV-cm^2) / electron
  G4double k1(0.0217615);
  G4double k2(0.50985);
  G4double alf = ee / 938.232; //::: use G4 mass
  G4double a2  = (alf + 1.0) * (alf + 1.0);
  G4double a3  = a2 - 1.0;
  G4double a4  = a3 / a2;
  G4double a5  = log(a3);
  G4double a6  = a5 - a4 + k1; 
  G4double del = adel(ee, sum2, sum1);
  G4double a7  = k2 * ((a6 - del / 2.0) * sum1 - sum2) / a4;
  return 1.0 / a7;
}

G4double G4BertiniCascade::adel(G4double e, 
			 G4double sum2, 
			 G4double sum1) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // computes asymptotic density effect correction
  const G4double xm    = 938.2;
  const G4double alpha = 1.378e-09;
  G4double del = log(alpha * sum1 * (e / xm) * (e / xm + 2.0)) - 2.0 * sum2 / sum1 - 1.0;
  if (del < 0.0) del = 0.0;
  return del;
}

void G4BertiniCascade::azirn(G4double &s,   // ::: move to utils
		      G4double &c) { 
  // decription: produce the sine and cosine of a random azimuthal angle   
  // parameters:
  // uses:
  // changes:
  G4double r1;
  G4double r2;
  G4double r1sq;  
  G4double r2sq;
  G4double rsq(2.0);  
  while (rsq > 1.0) {
    r1   = frnd(); // [-1,1]
    r1sq = sqr(r1);
    r2   = G4UniformRand();
    r2sq = sqr(r2);
    rsq  = r1sq + r2sq;
  }
  s = 2.0 * r1 * r2 / rsq;
  c = (r2sq - r1sq) / rsq;
}

void G4BertiniCascade::getflt() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int mt;
  totels = 0.0;
  mat = nmed[no];
  eps = exprnf();
  if (mat == 6666) {
    mt = mxmat + 1;
  } else {
    mt = mat;
    if (ityp == 2 && energy[no] <= elas && noel[mt] > 0) {
      if (ncol == 5) { // if previous event was not pseudo collision ...
	totels = savels;
      } else {
	j = 0;
	// interpolateElasticNeutronData(mt, j, energy[no]); //:::
	savels = totels;
      }
    }
  }
  G4double sigtrn = sigmx[ityp][mt] + totels;
  bodge = 0.0;
  // magnetic field. if field present, and particle charged add whopping amount to total sigma to ensure small stepsize
  // p.s.  this is the only place hfield gets called  
  if (ityp != 2) hfield(x[no], y[no], z[no], u[no], v[no], w[no], tip[no], energy[no]);
  sigtrn += bodge; // end magnetic change
  d = eps / sigtrn; // force collision for source particle
  if (ifircl >= 1 && name[no] == 1 && mat == 1) d = 0.001;
}

void G4BertiniCascade::hfield(G4double x, 
		       G4double y, 
		       G4double z, 
		       G4double u, 
		       G4double v,
		       G4double w, 
		       G4double tip, 
		       G4double e) {
  // decription: compute field vector (kilogauss) and set bodge if field nonzero 
  // parameters:
  // uses:
  // changes:
  G4double rmass[20];
  for (G4int i = 0; i < 20; ++i) rmass[i] = ginum[i] * MeV;
  bodge = 0.0;
  bfield[0] = 0.0;
  bfield[1] = 0.0;
  bfield[2] = 0.0;
  if (sqrt(sqr(x) + sqr(y)) > 3.047) return;
  if (z >= -590.0 && z <=  590.0) return;
  bodge = 1.0;
  bfield[0] = 39.4;
  bfield[1] = 0.0;
  bfield[2] = 0.0;
  it = static_cast<G4int>(tip + 1.01);
  G4double rm = rmass[it];
  // e is apparently the kinetic energy of particle in getflt calling hfield
  G4double etot = e + rm;
  G4double ee = e;
  if (kindo[no] >= 8) {
    it = kindo[no];
    rm = rmass[it];
    e = etot - rm;
  }
  G4double p = e * sqrt(1.0 + 2.0 * rm / e);
  bodge = 0.01 * p / bfield[1];
  bodge = 1.0 / bodge;
  if (bodge > 1.0) bodge = 1.0;
  e = ee;
}

void G4BertiniCascade::ck(G4int k, 
		   G4int npts, 
		   G4double *energy, 
		   G4double elas, 
		   G4double minEnergy) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  for (G4int i = 1; i < npts; ++i) {
    if (energy[i] < energy[i-1]) break;
    if (verboseLevel > 1) {
      G4cout << "elastic (N, N) energies out of order in ck: " << G4endl;
      G4cout << " i     : " << i      << G4endl;
      G4cout << " npts  : " << npts   << G4endl;
      G4cout << " k     : " << k      << G4endl;
      G4cout << " energy[i]  : " << energy[i]   << G4endl;
      G4cout << " energy[i-1]: " << energy[i-1] << G4endl;
    }
    G4Exception("ck1");
  }
  if (energy[0] >= elas) {
    if (energy[npts - 1] <= minEnergy) return;  //::: fix: following Exection is not done
    if (verboseLevel > 1) {
      G4cout << " neutron cutoff enrgy is below elastic (n, N) cross-section data in ck" << G4endl;
      G4cout << " npts        : " << npts        << G4endl;
      G4cout << " k           : " << k           << G4endl;
      G4cout << " energy[npts - 1] : " << energy[npts - 1] << G4endl; 
      G4cout << " minEnergy   : " << minEnergy   << G4endl;
    }
    G4Exception("ck1");
  }
  if (verboseLevel > 1) {
    G4cout << " elas out of elastic (n, N) cross-section range in ck: " << G4endl;
    G4cout << " npts : " << npts << G4endl;
    G4cout << " k    : " << k    << G4endl;
    G4cout << " energy[1] : " << energy[1] << G4endl; 
    G4cout << " elas : " << elas << G4endl;
  }
  G4Exception("ck1");
}

void G4BertiniCascade::rainge(G4double eng, 
		       G4int mat, 
		       G4int ityp) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double rng;
  static G4bool first = true;
  if (first) {
    first = false;
    for (G4int i = 0; i < 7; ++i) {
      if (i == 1 || i == 3) continue;
      rms[i]    = 1.0 / ginum[i];
      rmsinv[i] = ginum[i];
    }
  }
  if (ityp > 7 || ityp < 1) G4Exception("error in rainge");
  elowd = elow;
  if (mat == 6666) rng = 1.0e19;
  G4double e = eng;
  if (ityp > 1) e *= rms[ityp];
  G4int i = static_cast<G4int>(log(e / elowd) / einc);
  if (i >= 1) {
    fact = log(e / erng[i]) / log(erng[i + 1] / erng[i]);
    xbar = log( rnge[i][mat]) + fact * log( rnge[i + 1][mat] /  rnge[i][mat]);
    sigr = log(sigr2[i][mat]) + fact * log(sigr2[i + 1][mat] / sigr2[i][mat]);
    xbar = exp(xbar);
    sigr = exp(sigr);
  } else {
    xbar = log(e / elowd) / log(erng[0] / elowd) *  rnge[0][mat];
    sigr = log(e / elowd) / log(erng[0] / elowd) * sigr2[0][mat];
  } if (ityp > 1) {
    xbar = xbar * rmsinv[ityp];
    sigr = sigr * rmsinv[ityp];
  } 
  xbar = xbar - rngelo[ityp][mat];
  sigr = sigr - sigelo[ityp][mat];
  if (sigr <= 0.0) {
    rbar = xbar;
    rng = xbar;
  }
  rng = -1.0;
  while (rng < 0.0) {
    r1 = gaurn();
    if (G4UniformRand() <= 0.5) r1 = -r1;
    rng = r1 * sqrt(sigr) + xbar;
  }
  rbar = xbar;
}

void G4BertiniCascade::gthsig(G4int isgnl) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  // hsigmx  maximal cross-section :::
  G4int i;
  G4double ei;
  if (isgnl != 1) {
    it = ityp - static_cast<G4int>(ityp / 5);
    is = isa[it];
    sgm(2, it, is, de, ec[no], i, ei, hsig);
    return;
  }

  ifstream dataFile("chetc.dat");
  if (verboseLevel > 1) G4cout << "reading bertini data ..." << G4endl;
  if (dataFile.bad()){
    G4cerr << "chetc.dat file not found" << G4endl;
  }
  for (i = 0; i < 29850; i++) dataFile >> crossSection[i];
  if (verboseLevel > 1) G4cout << "bertini data OK" << G4endl;

  shxd();

  for (i = 0; i < 4; ++i) {
    isa[i] = locx[2][i];
    if (verboseLevel > 1) G4cout << "index isa[" << i << "]: " << isa[i] << G4endl;
  }
  // compute total (p, p) and (n, p) cross-sections (single + double production + elastic)
  // and store (p, p) beginning at bcs(3794) and (n, p) beginning bcs(3970)
  G4int isp = 3395; // 995 + 2400 index for  (p, p) single production 
  G4int idp = 3553; // 1153 + 2400) index for (p, p) double production 

  for (G4int it = 0; it < 2; ++it) {
    is = isa[it];
    cout << "is: " << is << G4endl;
    for (i = 0; i < 158; ++i) {
      crossSection[is + 18 + i] += crossSection[isp + i];
      cout << "to" << is + 18 + i << " added" << crossSection[isp + i] << " index " << isp + i << G4endl;
    }
    isp += 288;
    for (i = 0; i < 130; ++i) crossSection[is + 46 + i] += crossSection[idp + i];
    idp += 288;
  }
  
  if (verboseLevel > 1) {
    G4cout << "total (p, p) cross-section:" << G4endl;
    for (G4int i = 3792; i < 3970; ++i) cout <<  i << " " << crossSection[i] << G4endl;
  }
  
  if (verboseLevel > 1) {
    G4cout << "total (n, p) cross-section:" << G4endl;
    for (G4int i = 3968; i < 4150; ++i) cout <<  i << " " << crossSection[i] << G4endl;
  }
  // compute total pi+ (single production + elastic) and pi- (single production + exchange + elastic)
  // and store (pi+, p) beginning at bcs(3668) and (pi-, p) beginning at bcs(3542)
  isp = 4410; // 2009 + 2400 +1 this index points to ppscl vector
  for (it = 2; it < 4; ++it) {
    is = isa[it];
    for (G4int i = 0; i < 117; ++i) crossSection[is + 9 + i] += crossSection[isp + i];
    isp += 234;
  }
  G4int  iex = 3415 + 2400;
  is = isa[3];
  for (i = 0; i < 126; ++i) crossSection[is + i] += crossSection[iex + i];


  if (verboseLevel > 1) {
    G4cout << "total (pi+, p) cross-section:" << G4endl;
    for (G4int i = is-2; i  < (is - 2 +125)  ; ++i) cout <<  i << " " << crossSection[i] << G4endl;
  }
  
  for (it = 0; it < 2; ++it) {
    is = isa[it];
    for (G4int i = 0; i < 176; ++i) npsg[it][i] = crossSection[is + i];
  }
  for (it = 0; it < 2; ++it) {
    is = isa[it + 2];
    for (G4int i = 0; i < 126; ++i) pipsg[it][i] = crossSection[is + i];
  }
  if (emax < ehin) { // select maximum total cross-sections for x on prot. in energy range emin[x] to where x = p, n, pi+-
    emx[0] = emax;
    emx[1] = emax;
    emx[2] = 2500.0;
    emx[3] = 2500.0;
  } else {
    emx[0] = ehin;
    emx[1] = ehin;
    emx[2] = ehipi;
    emx[3] = ehipi;
  }
  for (G4int itp = 0; itp < 4; ++itp) {
    is = isa[itp];
    it = itp + itp / 4;
    sgm(1, itp, is, de, emin[it], il, el, sl); // *l = lower
    sgm(1, itp, is, de, emx[itp], ih, eh, sh); // *h = higher
    hsigmx[it] = max(sl, sh);
    if (il >= ih) break;
    ++il;
    for (G4int i = il - 1; i < ih; ++i) {
      if (crossSection[i+is] <= hsigmx[it]) break;
      hsigmx[it] = crossSection[i + is];
    }
  }
  hsigmx[3] = 0.0;

}

void G4BertiniCascade::readh(G4double *geosig) { 
  // description: set geometric cross-section
  //              read first 2400 numbers form chetc.dat and initialize geosig-vector
  // parameters:
  // uses:
  // changes:
  if (verboseLevel > 1) G4cout << "readh: reading bertini data ..." << G4endl;
  ifstream dataFile("chetc.dat");
  G4int k = 1;
  for (G4int j = 0; j < 4; ++j) {
    if (dataFile.bad()) {
      G4cerr << "chetc.dat file not found" << G4endl;
    }
    G4int i;
    for (i = 0; i < 600; i++) dataFile >> crsc[i];
  
    for (i = 3; i < 593; i = i + 10) {  // data organized into blocks of 10 numbers
      geosig[k] = crsc[i];
      geosig[k] = sqr(pi * geosig[k]) / millibarn;
      if (verboseLevel > 2) G4cout << "geosig: " << k << " " << geosig[k]<< G4endl;
      ++k;
    }
  }
  if (verboseLevel > 1) G4cout << "readh: geometric cross-section OK" << G4endl;
  return;
}

void G4BertiniCascade::sgm(
		    const G4int nsgnl, 	
		    G4int it, 
		    const G4int is, 
		    const G4double de, 	
		    const G4double em, 
		    G4int &i, 		    
		    G4double &e, 
		    G4double &s) {
  // decription: 
  // parameters:
  //      nsgnl :::
  //      it interaction type ?
  //      is used in cross-section vector to point correct partice data
  //      em energy used for cross-section calculation
  //      de difference of discrete energy points i cross-section tables 
  // output:
  //      e  fixed energy position used in cross-section tabulation
  //      s  cross-section at energy position em
  // uses: nsgml, it, is, de, em
  // changes: i, e, s
  G4int iemx = 9;
  i = static_cast<G4int>(em / de + 1.0);
  e = (i - 1) * de; 
  if (em >= 20.0 || it > 2) {
    if (nsgnl == 1) {
      s = crossSection[i + is] + (em - e) / de * (crossSection[i + 1 + is] - crossSection[i + is]); // interpolation
      return;
    }
    if (it > 2) {
      if (i > 125) {
	if (verboseLevel > 1) {
	  G4cout << " i = " << i << " in sgm" << G4endl;
	}
	G4cerr << "sgm1" << G4endl;
      }
      it = it - 2;
      s = pipsg[it][i] + (em - e) / de * (pipsg[it][i + 1] - pipsg[it][i]);
      // low energy ( < 20 MeV) (p, p) of (n, p) cross-sections for elastic scattering with h
      return;
    }
    if (i < 176) {
      s = npsg[it][i] + (em - e) / de * (npsg[it][i + 1] - npsg[it][i]);
      return;
    }    
    G4cerr << "sgm1" << G4endl;
  }
  if (em < enrgy[1]) G4cerr << "sgm2" << G4endl;
  for (G4int ie = 1; ie < iemx; ++ie) {
    if (em <= enrgy[ie]) {
      G4double s = log(ppnp[ie - 1][it]) + (em - enrgy[ie - 1]) /  (enrgy[ie] - enrgy[ie - 1]) * (log(ppnp[ie][it]) - log(ppnp[ie - 1][it])); // ::: indexes ?
      s = exp(s) * millibarn;
      return;
    }
  }
  G4cerr << "sgm3" << G4endl;
}

void G4BertiniCascade::recoil() {
  // decription : sums product particle momentum and sets recoil kinetic energy erec
  // parameters :
  // uses       :
  // changes    :
  if (nopart < 0 || ec[no] <= 0.0) {
    erec = 0.0;
    return;
  }
  if (apr <= 4.0) {
    erec = 0.0;
    return;
  }
  G4double pi;
  G4double px = 0.0;
  G4double py = 0.0;
  G4double pz = 0.0;
  G4double tm; // particle mass 
  if (nopart != 0) {
    for (G4int i = 0; i < nopart; ++i) {
      G4int kd = kind[i] + 1;
      switch(kd) {
      case 1:
      case 2:
	tm = 940.1;
	break;
      case 3:
      case 5:
	tm = 139.9;
	break;
      case 4:
	tm = 135.1;
      }
      if (ep[i] <= 0.0) {
	pi = 0.0;
	continue;
      } else {
	pi = ep[i] * sqrt(1.0 + 2.0 * tm / ep[i]);
	px = pi * alpha[i] + px;
	py = pi *  beta[i] + py;
      }
      pz = pi * gam[i] + pz;
    } 
  }
  G4int kt = static_cast<G4int>(tip[no] + 1.0);
  switch(kt) {
  case 1:
  case 2:
    tm = 940.1;
  case 3:
  case 5:
    tm = 139.9;
  case 4:
    tm = 135.1;
  }
  pz = ec[no] * sqrt(1.0 + 2.0 * tm / ec[no]) - pz;
  G4double aa = 931.162 * apr;
  erec = sqrt(sqr(aa) + sqr(px) + sqr(py) + sqr(pz)) - aa;
  return;
}

G4double G4BertiniCascade::ecol(G4double rr, 
			 G4int    mat, 
			 G4int    ityp, 
                         G4double ec,
			 G4double e) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  static G4double rmsinv[7], rms[7];
  static G4bool first = true;
  if (first) {
    for (G4int i = 0; i < 7; ++i) {
      if (i == 1 || i == 3) continue;
      rms[i]    = 1.0 / ginum[i];
      rmsinv[i] = ginum[i];
    }
  }
  if (ityp > 6 || ityp < 0) {
    G4cerr << "ecol" << G4endl;
  }
  G4double elowd = elow;
  if (mat == 6666) {
    return e;
  }
  G4double rngp = rr + rngelo[ityp][mat];
  if (ityp > 0) rngp *= rms[ityp];
  G4double ebar;
  for (G4int i = 0; i < 900; ++i) {
    if (rngp <= rnge[i][mat]) {
      if (i == 0) {
	ebar = log(elowd) + log(erng[1] / elowd) * rngp / rnge[1][mat];
      } else {
	ebar = log(erng[i - 1]) + log(erng[i] / erng[i - 1]) *
	  log(rngp / rnge[i - 1][mat]) / log(rnge[i][mat] / rnge[i - 1][mat]);
      }
      ebar = exp(ebar);
      if (ityp > 0) ebar *= rmsinv[ityp];
      return ebar;
    }
  }
  ebar = log(elowd) + log(erng[1] / elowd) * rngp / rnge[1][mat];
  ebar = exp(ebar);
  if (ityp > 0) ebar *= rmsinv[ityp];
  return ebar;
}

void G4BertiniCascade::analz1(int i){// :::functions to be made
;
}

void G4BertiniCascade::mfpd2(G4int i) {
  // description: mfpd2 performs pion decay to muona and neutrino 
  // by choosing random center-of-mass polar angle and random aximuthal angle
  // parameter: i = 2 reaction, i = 3 carged particle death, i = 5 pseudo collision                    
  ;
}

void G4BertiniCascade::sors(G4int i) {
  ;
}

G4int G4BertiniCascade::iclock() {
  // this is made to replace the ibm/370-lib function with same for cray
  // it should return iclock = cpu time in hundredths of seconds;
  // that is, time, t, in seconds: t = iclock() / 100.00.
  // call second(tsec)
  // iclock = 100.0 * tsec
  // for the ibm "cluster machine" in epm division:ornl:
  G4double tsec;
  // ::: mclock( tsec );
  return int(tsec);
}


void G4BertiniCascade::scatt(G4int inc) {
  // decription: calculates elastic scattering with hydrogen above scaling model treshhold
  // parameters:
  // uses:
  // changes:
  // rands is location of random number sequence
  // subroutines called G4UniformRand and azio
  // both are used in mecc
  // scatt(inc, ec[no], kind, ep, alpha, beta, gam); 
  G4double eke1 = ec[no]; 
  // itfas = static_cast<G4double>(kind); 
  //efas = ep; //::: 
  //alpfas = alpha; //::: 
  //betfas = beta; //::: 
  //gamfas = gam; //::: 
  G4double m2 = MeV * ginum[1];
  // NOTE: ginum values are in GeV but calculations here are done in MeV
  // it has been assumed the target particle is a p.
  G4double t2 = sqr(m2);
  ityp = inc + 1;
  ityp = kindo[noo];
  G4double m1;
  G4double t1;
  G4double p1;
  G4double b;
  if (icon == 1 || ityp <= 7 || ityp > 26) {
    switch(in) {
    case 1:
    case 2:
      m1 = ginum[1] * MeV;  //... MeV<>GeV
      break;
    case 3:
      b = 7.040e-06;
      m1 = ginum[5] * MeV;
      break;
    case 4:
      G4Exception("scatt1");
      break;
    case 5:
      b = 7.575e-06;
      m1 = ginum[5] * MeV;
    }
    e1 = eke1 + m1;
    t1 = sqr(m1);
    p1 = sqrt(sqr(e1) - t1);
    if (in <= 2) b = 7.26e-06 + 3.13e-11 * p1;
  } else {
    if (icon != 1) {
      G4int ino = in;
      switch(ino) {
      case 1:
      case 2:
	m1 = ginum[ityp] * MeV;
      case 3:
	b = 7.040e-06;
	m1 = ginum[ityp] * MeV;
      case 4:
	G4Exception("scatt2");
	break;
      case 5:
	b = 7.575e-06;
	m1 = ginum[ityp] * MeV;
      }
      e1 = eke1 + m1;
      t1 = sqr(m1);
      p1 = sqrt(sqr(e1) - t1);
      if (ino <= 2) b = 7.26e-06 + 3.13e-11 * p1;
    }
  }
  // note that the values of b have not been changed from the icon values
  // this should be rectified later, perhaps
  G4double etlab = e1 + m2;
  G4double beta = p1 / etlab;
  G4double gamma = 1.0 / sqrt(1.0 - sqr(beta));
  G4double etcm = etlab / gamma;
  G4double t7 = b * (pow(etcm, 4) - 2.0 * sqr(etcm) * (t1 + t2) + sqr(t1 - t2)) / (2.0 * sqr(etcm));
  G4double tb = t7;
  if (t7 > 50.0) tb = 50.0;
  cst = 1.0 + (log(1.0 - G4UniformRand() * (1.0 - exp( -tb)))) / t7;
  snt = sqrt(1.0 - sqr(cst));
  G4double t9;
  G4double t10;
  // azio(t9, t10); //::: 
  G4double ecm1 = (sqr(etcm) + t1 - t2) / (2.0 * etcm);
  G4double pcm1 = sqrt(sqr(ecm1) - t1);
  itfas[1] = inc;
  itfas[2] = 0.0;
  G4double t11 = gamma * (ecm1 + pcm1 * cst * beta) - m1;
  G4double t12 = etlab - t11 -m1 - m2;
  G4double plab1 = sqrt(2.0 * m1 * t11 + sqr(t11));
  G4double plab2 = sqrt(2.0 * m2 * t12 + sqr(t12));
  efas[1] = t11;
  efas[2] = t12;
  alpfas[1] = pcm1 * snt * t9;
  alpfas[2] = -alpfas[1];
  alpfas[1] = alpfas[1] / plab1;
  alpfas[2] = alpfas[2] / plab2;
  betfas[1] = pcm1 * snt * t10;
  betfas[2] = -betfas[1];
  betfas[1] = betfas[1] / plab1;
  betfas[2] = betfas[2] / plab2;
  gamfas[1] = gamma * (pcm1 * cst + beta * ecm1);
  gamfas[2] =  -gamfas[1] + p1;
  gamfas[1] = gamfas[1] / plab1;
  gamfas[2] =  gamfas[2] / plab2;
  return;
}



G4double G4BertiniCascade::frnd() {
  G4double result = 2.0 * G4UniformRand();
  if ( result <= 1.0 )return result;
  return 1.0 - result;
}


G4double G4BertiniCascade::ppnp[9][2] = 
{ 
  {1.20,		0.890},
  {0.620,               0.392},
  {0.310,		0.250},	
  {0.208,		0.178},	
  {0.155,		4.25},	
  {2.28,		1.62},	
  {1.14,		0.940},	
  {0.765,		0.645},	
  {0.555,		0.480 },
}; 

G4double G4BertiniCascade::enrgy[9] =
{ 
  1.0,  	3.0,  	5.0,  
  8.0,  	10.0, 	12.5,
  15.0, 	17.5,	20.0 
};

G4double G4BertiniCascade::ginum[30] = 
{
  938.2796e-03,	939.5731e-03,		139.5669e-03,
  134.9626e-03,	139.5669e-03,		105.65932e-03,
  105.65932e-03,493.688e-03,		497.67e-03,
  493.688e-03,	497.67e-03,             497.66e-03,
  497.66e-03, 	938.2796e-03,		939.5731e-03,
  1115.6e-03,   1115.6e-03,             1189.37e-03,
  1192.47e-03,	1197.35e-03,		0.5110034e-03,
  0.5110034e-03,	0.0,	        0.0,
  0.0,			0.0,		0.0,
  0.0,			0.0,	       0.0
};
	
G4double G4BertiniCascade::fli[] = 
{
  18.7, 	42.0, 	39.0, 
  60.0,	        68.0,	78.0,
  99.5,	        98.5,	117.0,
  140.0,	150.0,	157.0,
  163.0
};

G4double G4BertiniCascade::particleMass[7] = 
{
  938.3,        939.6,	 139.89,
  135.143,	139.89,	 105.7,
  105.7
};

G4double G4BertiniCascade::particleCharge[7] = 
{ 
  1.0,	   0.0,	   1.0,    0.0,	-1.0,	1.0,	-1.0
};
 
G4double G4BertiniCascade::massNucleon = 4.758e13;  
G4double G4BertiniCascade::massPionCharged = 7.08e12;
G4double G4BertiniCascade::massPionZero = 6.84e12;   
G4double G4BertiniCascade::oneThird = 1.0 / 3.0;
G4double G4BertiniCascade::twoThirds = 1.0 / 3.0;
G4double G4BertiniCascade::fourThirds = 4.0 / 3.0;







