#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4Nucleus.hh"
#include "G4BertiniInelasticCollision.hh"

G4BertiniInelasticCollision::G4BertiniInelasticCollision() {
  ;
}

G4BertiniInelasticCollision::~G4BertiniInelasticCollision(){
  ;
}


void G4BertiniInelasticCollision::bert(G4int ibert, 
		     G4double *finput, 
		     G4int nopart, 
		     G4int *kind, 
		     G4double *eray, 
		     G4double *aray, 
		     G4double *bray, 
		     G4double *gray) {
  G4double c3040 = 0.0;
  // description: a.c. 3526 (3410-44) cascade calculation
  G4double zero;
  G4double xinc;
  G4int    isix;
  G4double plvc1;
  G4int    ifive;
  G4int    k;
  G4int    m;
  G4int    klmn;
  G4int    i11;
  G4int    i31;
  G4int    i18;
  G4int    ne;
  G4int    i19;
  G4int    kk;
  G4int    i1i;
  G4int    i2i;
  G4int    icc[12];
  G4double abz;  
  --gray;
  amasno = 1.0;
  amasno = finput[1];
  zee    = finput[2];
  einc   = finput[3];
  ctofe  = finput[4];
  casesn = finput[5];
  andit  = finput[6];
  ctofen = finput[8];
  prtin  = finput[7];
  ke = 0;
  if (ibert != 0) goto L16;
  //::: al1.aerr = 0;
  //::: al1.aunit = nbertp;
  //::: frew(&al1);
  nrt = 0; // nrt number of files 
  sf = 1.0; // sf is a scale factor subject to change
  i1i = 0;
  i2i = 600;
  G4int j;
  for (j = 0; j < 4; ++j) {
    i1 = i2i;
    for (G4int i = i1i; i < i1; ++i) {
      // f77: 'read'
    }
    i1i += 600;
    i2i += 600; 
  }
  G4int i;
  for (i = 0; i < 29850; ++i) {
    // f77: 'read'
  }
  ln = 2;
  randi[0] = 16896;
  sqnm = sqr(massNucleon);  
  rcpmv = 5.0613e10; // reciprocal cm per mv  
  ifive = 5;
  isix  = 6;
  zero  = 0.0;
  ne = 0;
  for (i = 0; i < 19; ++i) {
    poac[i] += poac[i]; 
    ppac[i] /= sf;
  }
  // p0ac(19), ppac(19) 
  // andit 
  // isobar angular distribution [0, 50] % isotropic is 50% case 
  // percent forward-backward-1, all isotropic-2, all forward-backward is  0 % 
  // all of word in crdet to be considered.  not 0, only particle inpt 
  // escaping particle storage  esps 
 L16:
  for (i = 0; i < 60; ++i) ipec[i] = 0;
  i18 = 0;
  i19 = 0;
  for (i = 0; i < 2114; ++i) esps[i] = 0.0;
  // do15 i = 4515, 4849 add 3 for rands(4) & 3 for erand(4) -from i * 2  to i * 4 
  for (i = 4514; i < 4850; ++i) {
    // above is correct for code on epmnas cluster; 4848 is correct for i 
    // with int*2 for rands and erand, randi, 4854 is correct for cra 
    // with int*4 for rands, erand, and randi and double precision turned off  
    esps[i] = 0.0;
  }
  for (i = 0; i < 10; ++i) count[i] = 0.0;
  for (i = 0; i < 12; ++i) icc[i] = 0;
  space[12] = einc;
  no = static_cast<G4int>(amasno);
  nmas = (no - 1) * 10 + 1;
  for (i = 0; i < 10; ++i) {
    out[i] = crsc[nmas]; 
    ++nmas;
  }
  //  bertiniCascade->rout1(); //::: 
  if (prtin > 4.0) goto L135;
  no = static_cast<G4int>(prtin + 1.0); 
  value1 = einc + space[3];
  switch (no) {
  case 1:  goto L1000;
  case 2:  goto L1000;
  case 3:  goto L121;
  case 4:  goto L3015;
  case 5:  goto L121;
  }
 L3015:
  ne = 3; // incident particle is pi0, therefore in error 
  goto L135;
 L121:
  //  bertiniCascade->rout2(ppac); //:::
  if (i1 != 0) goto L200;
  ne = 4; // value1 > 2500 
  goto L135;
 L200:
  if (no < 4) {
    goto L3005;
  } else if (no == 4) {
    goto L3025;
  } else {
    goto L202;
  }
 L3025: 
  ne = 5;
  goto L135;
 L202:
  iv = 0;
 L203:
  xyi(1, 14, 30);
  ip = 1;
 L201:
  //  bertiniCascade->rout3(); //:::
  if (begru != 0.0) goto L220;
  else goto L5090;
  // if begru = 0, last run completed--bg6a 
 L220:
  kk = i1;
  xinc = coordinate[0]; // xinc is x-coordinate of incoming particle 
  // bertiniCascade->rout4(); //:::
  if (i1 >= 0) goto L224;
  else goto L1117;
 L224:
  i1 = kk;
 L225:
  if (in != 0) goto L4480;
  if (ex <= d[1]) goto L230;
  else goto L560;
 L230:
  curr[1] = out[12];
  wkrpn[2] = curr[1]; // kinetic energy with respect to protons region 3 
  wkrpn[5] = out[15]; // kinetic energy with respect to neutrons region 3 
  ifca = 0;
 L231:
  bg6ca(3, 3);
 L234:
  ifcc = 0; 
  abran(6);
  knot = not;

  if (not == 4) goto L315;
  else goto L235;
 L235:
  bg6c(isw[10]);
  value1 = rlke;
  if (in != 0) goto L4830;
  if (not < 4) {
    goto L23501;
  } else if (not == 4) {
    goto L315;
  } else {
    goto L23502;
  }
 L23501:
  any = space[not + 12];
  goto L236;
 L23502:
  any = s[not - 5];
 L23503:
  if (not < 5) {
    goto L236;
  } else if (not == 5) {
    goto L540;
  } else {
    goto L550;
  }
  // it = 1 - 6  pipps(20051), bg129(21011), pimpd(21004), pipnd(20644) 
 L236:
  // bertiniCascade->rout5(ppec, pmec, pmxc); //:::
 L254:
  if (clsm < 2.0) {
    goto L640;
  } else if (clsm == 2.0) {
    goto L580;
  }
  // (pi-, p) exchange scattering cross-section 
  if (value1 <= value2) goto L260;
  else goto L290;
 L260:
  if (isw[0] != 0) goto L275;
  ifc = ifcc + 1;
  if (in != 0) goto L4500;
 L269:
  c[2] = 0.0;
 L270:
  c[0] = curr[3];
  c[1] = curr[4];
  c[2] = c[2] + curr[5] + ex + d[0];
 L274:
  switch (it) {
  case 1:  goto L670;
  case 2:  goto L905;
  case 3:  goto L910;
  case 4:  goto L930;
  case 5:  goto L955;
  case 6:  goto L975;
  case 7:  goto L980;
  case 8:  goto L990;
  case 9:  goto L995;
  case 10:  goto L980;
  case 11:  goto L2000;
  case 12:  goto L4240;
  case 13:  goto L4245;
  case 14:  goto L995;
  case 15:  goto L4250;
  case 16:  goto L4280;
  case 17:  goto L4285;
  case 18:  goto L4325;
  case 19:  goto L4330;
  case 20:  goto L4360;
  case 21:  goto L4365;
  case 22:  goto L4410;
  case 23:  goto L5050;
  case 24:  goto L5070;
  case 25:  goto L5075;
  case 26:  goto L5080;
  case 27:  goto L980;
  case 28:  goto L5085;
  }
 L275:
  if (isw[1] != 0) goto L285;
  ifc = ifcc + 2; 
  if (in != 0) goto L4515;
 L284:
  c[2] = d[1] + d[2];
  goto L270;
 L285:
  ifc = ifcc + 3;
  if (in != 0) goto L480;
 L289:
  c[2] = d[1] + d[2] + d[3] + d[4];
  goto L270;
  // ifc(1-3), bg6c(1502), bg6f(3243), bg6k(4055) collision 
 L290:
  signex();
 L291:
  if (isw[0] == 0) goto L225;
  if (in != 0) goto L4535;
  if (ex <= d[5]) goto L230;
  if (isw[1] != 0) goto L310;
 L305:
  ++ipec[6]; // number of escaped particles on region 2 
  goto L201;
 L310:
  ++ipec[10]; // number of escaped particles on region 1 
  goto L201;
 L345:
  i3 = 1;
  goto L316;
 L315:
  i3 = -1;
 L316:
  // bertiniCascade->rout6(); //:::
  if (i3 < 0) {
    goto L3110;
  } else if (i3 == 0) {
    goto L350;
  } else {
    goto L365;
  }
 L3110:
  ne = 22; // impossible for 6--7, curr(1) is < 3 or > 5 
  goto L135;
 L350:
  // bertiniCascade->rout6a(); //:::
  if (clsm < 2.0) {
    goto L655;
  } else if (clsm == 2.0) {
    goto L660;
  } else {
    goto L4870;
  }
 L356:
  ifca = 1; 
  if (isw[0] != 0) goto L360;
  else goto L269;
 L360:
  if (isw[1] != 0) goto L289;
  else  goto L284; 
  // non-deuteron absorption 
 L365:
  // bertiniCascade->rout7(); //:::
  if (i3 >= 0) goto L425;
  else goto L3110;
 L425:
  //  bertiniCascade->rout7a(); //:::
  i3 = i3; //:::
  switch (i3) {
  case 1:  goto L480;
  case 2:  goto L3030;
  case 3:  goto L535;
  case 4:  goto L510;
  case 5:  goto L269;
  case 6:  goto L289;
  case 7:  goto L485;
  case 8:  goto L490;
  case 9:  goto L284;
  }
 L3030:
  ne = 6; // pwd or nwd < 7.0 
  goto L135;
 L480:
  value1 = ex + d[3] + d[4];
 L484:
  if (curr[9] - 2.0 >= 0.0) goto L490;
 L485:
  c[0] = value1 * curr[6] + curr[3];
  c[1] = value1 * curr[7] + curr[4];
  c[2] = value1 * curr[8] + curr[5];
  goto L274;
 L490:
  value1 += d[2];
  goto L500;
 L500:
  if (curr[9] <= 2.0) goto L485;
  value1 += d[1];
  goto L485;
 L510:
  if (inc != 0) goto L515;
  else goto L525;
 L515:
  c[2] = d[1];
  if (isw[2] != 0) goto L520; 
  else goto L270;
 L520:
  c[2] = c[2] + d[2] + d[3];
  goto L270;
 L525:
  if (isw[2] != 0) goto L530; 
  else goto L500;
 L530:
  value1 = ex + d[3];
  goto L484;
 L535:
  value1 = ex;
  if (inc != 0) goto L284;
  else goto L484;
  
 L540:
  if (rlke > 2500.0) {
    ne = 26; // rlke > 2500 
    goto L135;
  }
  if (rlke > 180.0) goto L545;
 L541:
  signex();
  if (clsm < 2.0) {
    goto L629;
  } else if (clsm == 2.0) {
    goto L601;
  } else {
    goto L291;
  }
 L545:
  value1 += -180.0;
  crjab(1, ppscl[1]); // (pi+, p) single production cross-section low energy 
  goto L254;
 L550:
  if (rlke <= 2500.0) goto L554;
  ne = 27; // rlke > 2500 
  goto L135;
 L554:
  if (rlke <= 180.0) goto L541;
  value1 += -180.0;
  crjab(1, pmscl[1]); // (pi-, p) single production cross-section low energy 
  goto L254;
 L560:
  if (d[2] != 0.0) goto L570;
  ++ipec[1]; // number of particles incident on region 3 escaping 
  goto L201;
 L570:
  isw[0] = 1;
  spac32(31);
 L574:
  if (in != 0) goto L4825;
  if (ex <= d[2]) goto L575;
  else goto L615;
 L5750:
  if (in != 0) goto L576;
 L575:
  curr[1] = out[13];
  wkrpn[1] = curr[1];
  wkrpn[4] = out[16]; // kinetic energy for protons and neutrons in region 2 
 L576:
  bg6ca(2, 2);
  goto L234;
 L585:
  iv = -1;
  goto L581;
 L580:
  iv = 0;
 L581:
  //  bertiniCascade->rout8(); //:::
  i3 = i3;
  switch (i3) {
  case 1:  goto L601;
  case 2:  goto L4640;
  case 3:  goto L270;
  case 4:  goto L530;
  }
 L601:
  if (isw[2] == 0) goto L574;
 L605:
  if (ex <= d[4]) goto L5750;
  if (in != 0) goto L4840;
  spac32(32);
 L611:
  if (ex <= d[5]) goto L230;
  else goto L310;
 L615:
  if (d[3] != 0.0) goto L625;
  spac32(32);
 L621:
  if (ex <= d[5]) goto L230;
  else goto L305;
 L625:
  isw[1] = 1;
  isw[2] = 1;
  spac32(30);
 L629:
  if (in != 0) goto L4855;
 L6290:
  //bertiniCascade->rou10();//:::
  if (i3 >= 0) goto L636; 
  else goto L605;
 L636:
  bg6ca(1, 1);
  goto L234;
 L640:
  if (value1 > value2) goto L650; 
 L645:
  ifc = ifcc + 9;
  if (in != 0) goto L4620;
  else goto L284;
 L650:
  signex();
  if (in != 0) goto L4855;
  else goto L6290;
 L655:
  if (in == 0) {
    ifca = 6;
    goto L284;
  }
  ifca = abs(i6 - 2) * 9 + (i6 - 1) * 13 * (3 - i6);
  goto L4620;
 L660:
  if (in != 0) goto L4885;
  ifca = 7;
  goto L515;
 L670:
  i3 = 1;
  goto L673;
 L6700:
  i3 = 4;
  goto L673;
 L671:
  i3 = 2;
  goto L673;
 L672:
  i3 = 3;
 L673:
  // ::: compile hcol.f bertiniCascade->rou11();
  // ppdcl(378) 
  i3 = i3; // :::
  switch (i3) {
  case 1:  goto L865;
  case 2:  goto L925;
  case 3:  goto L970;
  case 4:  goto L677;
  }
 L677:
  cst = crdt[1] - fabs(snt * (crdt[1] - crdt[0]));
 L680:
  snt = sqrt(1.0 - sqr(cst));
 L681:
  // bertiniCascade->rou12();//:::
  if (i3 == 0) goto L685;
  if (i3 > 0) goto L710;
  ne = 7; // com < -5.0e-06 
  goto L135;
 L685:
  if (efrn < value1) goto L720;
  fcn += 1.0; // fcn is the number of forbidden collisions for neutrons
 L695:
  iv = -1;
  goto L696;
 L705:
  iv = 0;
 L696:
  i1 = 0;
  //  bertiniCascade->rou13(); //:::
  if (i3 < 0) goto L425;
  if (i3 > 0) goto L350;  
  ifc = ifc;
  switch (ifc) {
  case 1:  goto L225;
  case 2:  goto L621;
  case 3:  goto L611;
  case 4:  goto L1119;
  case 5:  goto L1236;
  case 6:  goto L1251;
  case 7:  goto L574;
  case 8:  goto L605;
  case 9:  goto L629;
  case 10:  goto L1219;
  case 11:  goto L1249;
  case 12:  goto L1244;
  case 13:  goto L4480;
  case 14:  goto L4535;
  case 15:  goto L4535;
  case 16:  goto L4611;
  case 17:  goto L4544;
  case 18:  goto L4631;
  case 19:  goto L4825;
  case 20:  goto L605;
  case 21:  goto L4480;
  case 22:  goto L4535;
  case 23:  goto L4535;
  case 24:  goto L4855;
  case 25:  goto L4480;
  case 26:  goto L4535;
  case 27:  goto L4535;
  case 28:  goto L4825;
  case 29:  goto L4955;
  case 30:  goto L4855;
  }
 L710:
  if (efrp < value1) goto L720;
  fcp += 1.0; // fcp is the number of forbidden collisions for protons
  goto L695;
 L850:
  i3 = 0;
  goto L721;
 L861:
  i3 = -1;
  goto L721;
 L720:
  i3 = 1;
 L721:
  i11 = i1;
  i31 = i3;
  plvc1 = plvc[0];  // plvc is particle with velocity less than criterion
  //  bertiniCascade->rou14();//:::
  i3 = i3;
  switch (i3) {
  case 1:  goto L710;
  case 2:  goto L685;
  case 3:  goto L3040;
  case 4:  goto L1340;
  case 5:  goto L4415;
  }
 L3040:
  ne = 8;
  // pinst(med-)  punp(pgvc or plvc value)  stpr(pgvc > 40  or plvc > 80) 
  c3040 += 1;
  if (c3040 > 100.0) goto L135;
  return;
 L865:
  i3 = 1;
  goto L866;
 L882:
  i3 = 3;
  goto L866;
 L883:
  i3 = 4;
 L866:
  //::: bertiniCascade->rou15();
  // hppdci(45), ppdci(170) 
  switch (i3) {
  case 1:  goto L915;
  case 2:  goto L960;
  case 3:  goto L680;
  case 4:  goto L3045;
  case 5:  goto L920;
  case 6:  goto L965;
  }
 L3045:
  ne = 9;
  // snn(rlke >= 1000) dcintp(rlke >= cross-section values) 
  goto L135;
 L905:
  pt[1] = 5.0;
  ik = it; // bg129  (pi-, n) 
  pt[13] = 2.0;
  goto L671;
 L910:
  pt[1] = 5.0; // pipnx  dir. scatterin 
  pt[13] = 1.0;
  ik = it;
  goto L671;
 L915:
  i3 = 1;
  goto L926;
 L920:
  i3 = 2;
  goto L926;
 L925:
  i3 = 3;
 L926:
  //::: bertiniCascade->rou16(pmdd);
  // hpmddi(45), pmddi(170), pmddl(378) 
 L927:
  if (i3 != 0) goto L680;
  else goto L882;
 L930:
  pt[1] = 3.0;
  pt[13] = 2.0;
  ik = 3;
  goto L671;
 L955:
  pt[13] = 2.0;
 L9550:
  ik = it;
 L956:
  pt[1] = 4.0;
  massParticle[2] = massPionZero; // pi0 mass/cm 
  goto L672;
 L960:
  if (ik == 23) goto L5055;
  i3 = 1;
  goto L972;
 L966:
  i3 = 2;
  goto L972;
 L971:
  i3 = 3;
 L972:
  //::: bertiniCascade->rou16(pmdx);
  // hpmdxi(45), pmdxi(170), pmdxl(378) 
  goto L927;
 L965:
  if (ik != 23) goto L966;
  else goto L5060;
 L970:
  if (ik != 23) goto L971;
  else goto L5065; // (pi-, p) cross-section
 L975:
  pt[13] = 1.0;
  goto L9550;
 L980:
  pt[1] = 1.0; // pi- (pp) abs, it = 10, pi+ (nn)  abs 
 L985:
  pt[13] = 2.0;
 L986:
  pol1(cst, snt);
  goto L681;
 L990:
  pt[1] = 2.0; // pi- (nn)  abs 
  goto L985;
 L995:
  pt[1] = 1.0; // pi- (pp)  abs, also pi+ abs 
  pt[13] = 1.0;
  goto L986;
 L2000:
  isw[8] = 0;
  isw[9] = 0; 
  i3 = 0;
  goto L2002;
 L2001:
  i3 = -1;
 L2002:
  //::: bertiniCascade->rou17(fripn, pnmi, fmxsp, pcfsl, pnfsl);
  // fripn(117), pnmi(101), fmxsp(117), pcfsl(234), pnfsl(234) 
  if (i3 > 0) goto L2055;
  if (i3 == 0) goto L2084; 
  ne = 10; // coll(com < -5.0e-06) ecpl(error in curr , strkp, pt(26), pt(2), pt(14) or pt(37)) isw10 = 0 
  goto L135;
 L2055:
  k = 3;
 L2056:
  if (pt[k - 1] != 1.0) goto L2085;
  if (pt[k] <= 0.0) goto L2070;
  if (pt[k] > efrp) goto L2075;
 L2070:
  fcp += 1.0; // number of forbidden collisions involving protons 
 L2071:
  massParticle[3] = massNucleon;
  goto L705;
 L2075:
  m = static_cast<G4int>(pt[k - 1]);
  if (pt[k] > eco[m]) goto L2081;
  pt[k] = 0.0;
  pnbc[m] += 1.0;
 L2081:
  if (col[14] < 1.0) { 
    goto L2084;
  } else if (col[14] == 1.0) {
    goto L4035;
  }
  if (col[14] < 3.0) { 
    goto L4071;
  } else if (col[14] == 3.0) {
    goto L4070;
  } 
  if (col[14] < 5.0) {
    goto L4080;
  } else if (col[14] == 5.0) goto L4145;
  ne = 11; // col(15) > 5 
  goto L135;
 L2084:
  collm(0);
  if (pt[37] != 0.0) goto L4015;
  else goto L4010;
 L2085:
  if (pt[k - 1] != 2.0) goto L4225;
  if (pt[k] <= 0.0) goto L4000;
  if (pt[k] > efrn) goto L2075;
 L4000:
  fcn += 1.0;
  goto L2071;
 L4010:
  i3 = 1;
  goto L4036;
 L4015:
  i3 = 2;
  goto L4036;
 L4070:
  i3 = 4;
  goto L4036;
 L4071:
  i3 = 5;
  goto L4036;
 L4080:
  i3 = 6;
  goto L4036;
 L4035:
  i3 = 3;
 L4036:
  //  bertiniCascade->rou18(); //:::
  i3 = i3;
  k = iv;
  switch (i3) {
  case 1:  goto L2056;
  case 2:  goto L3060;
  case 3:  goto L2081;
  case 4:  goto L4235;
  case 5:  goto L3061;
  }
 L3060:
  ne = 12; // pt[k-1] <= 2 
  goto L135;
 L3061:
  ++i18;
  goto L2071;
 L4145:
  // bertiniCascade->rou19(); //:::
  if (i3 == 0) goto L4235;
  if (i3 > 0) goto L3066;
  ne = 13; // pt[k-1] < 3, > 6 - k < 27 
  goto L135;
 L3066:
  ++i19;
  goto L2071;
 L4225:
  if (col[14] < 1.0) goto L2084;
  ne = 14; // col(15) >= 1.0 
  goto L135;
 L4240:
  i3 = 2;
  goto L4236;
 L4245:
  i3 = 3;
  goto L4236;
 L4250:
  i3 = 4;
  goto L4236;
 L4280:
  i3 = 5;
  goto L4236;
 L4285:
  i3 = 6;
  goto L4236;
 L4325:
  i3 = 7;
  goto L4236;
 L4241:
  i3 = 8;
  goto L4236;
 L4235:
  i3 = 1;
 L4236:
  //:::bertiniCascade->rou20(dcin, dcln, dchn, pdci, pdch);
  i3 = i3;
  switch (i3) {
  case 1:  goto L3075;
  case 2:  goto L861;
  case 3:  goto L2001;
  case 4:  goto L985;
  case 5:  goto L883;
  case 6:  goto L680;
  case 7:  goto L986;
  case 8:  goto L681;
  }
 L3075:
  ne = 15; // pinst -med-  rlke > 3500 
  goto L135;
 L4360:
  i3 = 2;
  goto L4341;
 L4365:
  i3 = 3;
  goto L4341;
 L4410:
  i3 = 4;
  goto L4341;
 L4330:
  i3 = 1;
 L4341:
  // bertiniCascade->rou21(frinn, dmin, fmxsn, fmxdn, fsln); //:::
  goto L2002;
 L4830:
  iv = 2;
  goto L4481;
 L4500:
  iv = 3;
  goto L4481;
 L4515:
  iv = 4;
  goto L4481;
 L4535:
  iv = 5;
  goto L4481;
 L4825:
  iv = 6;
  goto L4481;
 L4640:
  iv = 7;
  goto L4481;
 L4840:
  iv = 8;
  goto L4481;
 L4855:
  iv = 9;
  goto L4481;
 L4620:
  iv = 10;
  goto L4481;
 L4885:
  iv = 11;
  goto L4481;
 L4611:
  iv = 12;
  goto L4481;
 L4544:
  iv = 13;
  goto L4481;
 L4631:
  iv = 14;
  goto L4481;
 L4955:
  iv = 15;
  goto L4481;
 L4415:
  iv = 16;
  goto L4481;
 L4696:
  iv = 17;
  goto L4481;
 L4670:
  iv = 18;
  goto L4481;
 L4479:
  iv = 19;
  goto L4481;
 L4650:
  iv = 20;
  goto L4481;
 L4610:
  iv = 21;
  goto L4481;
 L4870:
  iv = 22;
  goto L4481;
 L5005:
  iv = 23;
  goto L4481;
 L4480:
  iv = 1;
 L4481:
  // bertiniCascade->rou22(ppac, poac, pnec, pmxc, pnnec);
  iv = iv;
  if (i1 < 0) goto L1117;
  switch (iv) {
  case 1:  goto L5045;
  case 2:  goto L3080;
  case 3:  goto L605;
  case 4:  goto L850;
  case 5:  goto L1141;
  case 6:  goto L485;
  case 7:  goto L500;
  case 8:  goto L484;
  case 9:  goto L1290;
  case 10:  goto L1161;
  case 11:  goto L1157;
  case 12:  goto L1270;
  case 13:  goto L576;
  case 14:  goto L23503;
  case 15:  goto L636;
  case 16:  goto L356;
  case 17:  goto L480;
  case 18:  goto L530;
  case 19:  goto L11410;
  case 20:  goto L231;
  case 21:  goto L5035;
  case 22:  goto L1170;
  case 23:  goto L3105;
  }
 L3080:
  ne = 16;
  // com > 3500, ccpes, esps(1) >= 30.0, com > 2500 
  goto L135;
 L3105:
  ne = 21; // iv > 22 
  goto L135;
 L5010:
  abz = 1.0;
  value1 = G4UniformRand();
  if (value1 >= ppnda) goto L365; // probability (pi-, D) abs 
  it = 27; // bg117(20040) pi0 abs 
  absec = -hvp[medium - 1];
  goto L345;
 L5020:
  if (rlke <= 2500.0) goto L5021;
  ne = 25; // rlke > 2500 
  goto L135;
 L5021:
  if (rlke <= 180.0) goto L5030;
  if (not > 6) goto L5040;
  value1 = rlke - 180.0;
  crjab(1, pnscl[1]); // (pi-, p) single production cross-section low energy 
  goto L1170;
 L5030:
  if (clsm < 2.0) {
    goto L1335;
  } else if (clsm == 2.0) {
    goto L1315;
  } else {
    goto L1180;
  }
 L5035:
  if (not >= 6) goto L5020;
  else goto L5005;
 L5040:
  value1 = rlke - 180.0;
  crjab(1, pnnsl[1]); // (pi-, n) single production cross-section low enenergy 
  goto L1170;
 L5045:
  isw[10] = 0;
  goto L4696;
 L5050:
  pt[13] = 1.0;
  goto L9550;
 L5060:
  i3 = 2;
  goto L5056;
 L5065:
  i3 = 3;
  goto L5056;
 L5055:
  i3 = 1;
 L5056:
  //::: bertiniCascade->rou16(pndd);
  // hpnddi(45)
  // pnddi(170) (pi-, p) direct cross-section intermediate energy
  // pnddl(378) (pi-, p) direct differential cross-section low energy  
  goto L927;
 L5070:
  pt[1]  = 3.0;
  pt[13] = 2.0;
  ik = it;
  goto L671;
 L5075:
  pt[13] = 2.0;
  ik = 23;
  goto L956;
 L5080:
  pt[1]  = 5.0;
  goto L6700;
 L5085:
  isw[8] = 2;
  goto L4241;
 L5090:
  itote = ipec[1] + ipec[6] + ipec[10];
  i1 = itote - 1;
  if (i1 < 0) {
    goto L7004;
  } else if (i1 == 0) {
    goto L7002;
  }
  G4Exception("bert1"); //::: fix bug ?
 L7002:
  nopart = -1;
 L7003:
  return;
 L7004:
  nopart = static_cast<G4int>(esps[0]);
  if (nopart <= 350) goto L7007; //::: upper limit for number of particles
  nopart = 350;
 L7007:
  if (nopart <= 0) goto L7003;
  i1 = nopart;
  G4int ndex;
  for (ndex = 0; ndex < i1; ++ndex) {
    klmn = 8 * (ndex - 1) + 1;
    kind[ndex] = static_cast<G4int>(esps[klmn] - 1.0);
    eray[ndex] = esps[klmn + 1];
    aray[ndex] = esps[klmn + 2];
    bray[ndex] = esps[klmn + 3];
    gray[ndex] = esps[klmn + 4];
  }
  goto L7003;
 L135:
  G4cout << " ne = " << ne << G4endl;
  G4Exception("bert2");
 L1000:
  value2 = einc + space[11];
  if (value1 > 160.0) goto L1015;
  space[32] = 1.4e-24; // no production possible 
  fmax[1]   = 1.4e-24;
  space[33] = 4.6e-25;
  fmax[0]   = 4.6e-25;
  for (i = 8; i < 12; ++i) s[i] = 0.0;
  // einc + 50.0 < 160.0 
  goto L1090;
 L1015:
  ans = bovera(value2, massNucleon); // nucleon mass = constant ans = p1 / e1 
  if (value1 > 560.0) goto L1055;
  s[10] = 0.0;
  s[11] = 0.0; // s[11] single production possible, s[12] double production 
  if (value1 < 400.0) goto L1030;
  s[8] = ans * 2.26e-26;
  s[9] = ans * 1.4e-26;
  space[43] = 5.6e-26;
  space[44] = ans * 2.7e-26;
  goto L1054;
 L1030:
  if (value1 < 300.0) goto L1040;
  s[8] = ans * 2e-26;
  s[9] = ans * 1.4e-26;
  space[43] = 1.06e-25;
  space[44] = ans * 3.6e-26;
  goto L1054;
 L1040:
  if (value1 < 200.0) goto L1050;
  s[8] = ans * 1.14e-26;
  s[9] = ans * 1.12e-26;
  space[43] = 3.13e-25;
  space[44] = 1.03e-25;
  goto L1054;
 L1050:
  s[8] = ans * 1.95e-27;
  s[9] = ans * 1.7e-27;
  space[43] = 5.2e-25;
  space[44] = 1.76e-25;
 L1054:
  space[32] = space[43];
  space[33] = space[44];
  goto L1090;
 L1055:
  if (value1 > 3600.0){
    ne = 17; // value1 > 3600 
    goto L135;
  }
  s[8] = ans * 2.26e-26; // double production possible 
  s[9] = ans * 1.4e-26;
  if (value1  >= 800.0) goto L1070;
  s[10] = ans * 1.9e-27;
  s[11] = ans * 9e-27;
  space[45] = ans * 3.84e-26;
  space[46] = ans * 2.72e-26;
  goto L1089;
 L1070:
  if (value1 >= 1680.0 ) goto L1085;
  s[10] = ans * 1.08e-26;
  s[11] = ans * 1.74e-26;
  space[45] = ans * 3.3e-26;
  space[46] = ans * 2.72e-26;
  goto L1089;
 L1085:
  space[45] = ans * 2.5e-26;
  space[46] = ans * 2.65e-26;
  s[9] = ans * 1.36e-26;
  s[10] = ans * 1.8e-26;
  s[11] = ans * 2.36e-26;
 L1089:
  space[32] = space[45];
  space[33] = space[46];
 L1090:
  switch (no) {
  case 1:  goto L1091;
  case 2:  goto L2100;
  }
 L1091:
  iv = 1;
 L1092:
  xyi(9, 33, 41);
  ip = 2;
 L1095:
  if (no >= 2) goto L1100;
  else goto L1105;
 L1100:
  isw[3] = 0;
  goto L1110;
 L1105:
  isw[3] = 1;
 L1110:
  undis();
  if (begru == 0.0) goto L5090;
  abz = 0.0;
  xinc = coordinate[0]; // xinc is x-coordinate inc. particle 
  inc = 1; // 0 if particle cascade 
  curr[0] =  no;
  curr[2] = massNucleon;  
  geo();
  if (i1 >= 0) goto L1118;
 L1117:
  ne = 20; // error in geometry
  goto L135;
 L1118:
  partin();
  spac32(43);
 L1119:
  if (ex <= d[1]) goto L1120;
  else goto L1205;
 L1120:
  wkrpn[2] = out[12];
  wkrpn[5] = out[15];
  curr[1] = wkrpn[5];
  if (isw[3] == 0) goto L1135;
  curr[1] = wkrpn[2]; // kinetic energy with respect to neutrons (protons) in region 3 
 L1135:
  bg6ca(3, 0);
  ifca = 3;
 L1140:
  ifcc = 3;
 L1141:
  ka = 6;
 L11410:
  abran(ka);
  knot = not + 6;
  if (in <= 0) goto L1143;
  knot += 6;
 L1143:
  if (knot == 17) goto L5010;
  bg6c(isw[3]);
  if (rlke > 0.0) goto L1146;
  ne = 24; // rlke <= 0 
  goto L135;
 L1146:
  value1 = rlke;
  if (in  != 0) goto L4670;
  if (not >= 5) goto L1290;
  if (not >  2) goto L1270;
  if (not == 2) goto L1160;
  any = space[32];
 L1157:
  crjab(1, ecn[1]); // (n, p) elastic cross-section scattering 
  goto L1170;
 L1160:
  any = space[33];
 L1161:
  crjab(1, pec[1]); // (p, p) elastic scattering crossection 
 L1170:
  if (clsm < 2.0) {
    goto L1330;
  } else if (clsm == 2.0) {
    goto L1310;
  }
  if (value1 <= value2) goto L260;
 L1180:
  signex();
  if (isw[0] != 0) goto L1184;
  if (in == 0) goto L1119;
  if (curr[0] <= 2.0) goto L4480;
  else goto L4479;
 L1184:
  if (in !=     0) goto L4535;
  if (ex <=  d[5]) goto L1120;
  if (isw[1] != 0) goto L1200;
 L1195:
  ++ipec[6]; // number of escaped particles on region 2 
  goto L1095;
 L1200:
  ++ipec[10]; // number of particles escaped on region 1 
  goto L1095; 
 L1205:
  if (d[2] != 0.0) goto L1215;
  ++ipec[1];  // number of particles incident on region 3 escaping 
  goto L1095;
 L1215:
  isw[0] = 1;
  spac32(42);
 L1219:
  if (ex > d[2]) goto L1230;
 L1220:
  wkrpn[1] = out[13];
  wkrpn[4] = out[16];
  curr [1] = wkrpn[4];
  if (isw[3] == 0) goto L1225;
  curr[1] = wkrpn[1]; // kinetic energy with respect to neutrons (protons) region 2 
 L1225:
  bg6ca(2, 0);
  goto L1140;
 L1230:
  if (d[3] != 0.0) goto L1240;
  spac32(43);
 L1236:
  if (ex <= d[5]) goto L1120;
  else goto L1195;
 L1240:
  isw[1] = 1;
  isw[2] = 1;
  spac32(41);
 L1244:
  if (ex <= d[3]) goto L1255;
  spac32(42);
 L1249:
  if (ex <= d[4]) goto L1220;
  spac32(43);
 L1251:
  if (ex <= d[5]) goto L1120;
  else goto L1200;
 L1255:
  wkrpn[0] = out[14];
  wkrpn[3] = out[17];
  curr[1] = wkrpn[3];
  if (isw[3] == 0) goto L1265;
  curr[1] = wkrpn[0]; // kinetic energy with respect to neutrons (protons) region 1 
 L1265:
  bg6ca(1, 0);
  goto L1140;
 L1270:
  if (rlke <= 3500.0) goto L1274;
 L3115:
  ne = 23; // rlke > 3500 
  goto L135;
 L1274:
  if (rlke <= 360.0) goto L1325;
  value1 = rlke - 360.0;
  if (in != 0) goto L1278;
  any = s[knot - 1];
 L1278:
  if (not < 4) {
    goto L1280;
  } else if (not == 4) {
    goto L1285;
  } 
  ne =18;
  // not = 5 
  goto L135;
 L1280:
  crjab(1, pspcl[1]); // (p, p) single production cross-section low energy 
  goto L1170;
 L1285:
  crjab(1, spcln[1]); // (n, p) single production cross-section low energy 
  goto L1170;
 L1290:
  if (rlke > 3500.0) goto L3115;
  if (rlke <= 920.0) goto L1325;
  value1 = rlke - 920.0;
  if (not < 6) {
    goto L1300;
  } else if (not == 6) {
    goto L1303;
  }
  ne = 19; // not > 6 
  goto L135;
 L1300:
  if (in != 0) goto L1302;
  any = s[10];
 L1302:
  crjab(1, pdpcl[1]); // (p, p) double production cross-section low energy 
  goto L1170;
 L1303:
  if (in != 0) goto L1305;
  any = s[11];
 L1305:
  crjab(1, dpcln[1]); // (n, p) double production cross-section low energy 
  goto L1170;
 L1310:
  if (value1 <= value2) goto L585;
 L1315:
  signex();
  if (in != 0) goto L4650;
  if (isw[2] != 0) goto L1249;
  else goto L1219;
 L1325:
  if (clsm < 2.0) {
    goto L1335;
  } else if (clsm == 2.0) {
    goto L1315;
  } else {
    goto L1180;
  }
 L1330:
  if (value1 <= value2) goto L645;
 L1335:
  signex();
  if (in != 0) goto L4610;
  else goto L1244;
 L1340:
  if (esps[0] != 0.0) goto L1345;
  nwds = 1;
  goto L1349;
 L1345:
  nwds = static_cast<G4int>(esps[0] * 8. + 1.5); // total number of words (escaping particles) 
  if (count[5] >= 0.0) goto L1349; // minus, record not representative, skip 
  i1 = nwds;
  for (i = 0; i < i1; ++i) esps[i] = 0.0;
  for (i = 0; i < 5; ++i) count[i] = 0.0;
 L1349:
  ++nor; // record number  
  in = 0;
  switch (ip) {
  case 1:  goto L201;
  case 2:  goto L1095;
  }
 L2100:
  iv = -1;
  goto L1092;
 L3005:
  iv = 2;
  goto L203;
  // 1370 return 
}

void G4BertiniInelasticCollision::stor(G4double p0, 
		     G4int i, 
		     G4double e, 
		     G4double a, 
		     G4double b, 
		     G4double g, 
		     G4double w, 
		     G4int j, 
		     G4double erem, 
		     G4double ps, 
		     G4double wmass, 
		     G4int itype) {
  // description: 
  G4int istart(0);
  G4double ampi;
  G4int incp;
  G4double tamp;
  G4double plet;
  G4double delp0;
  G4double ampi2;
  G4int k;
  G4double delps[20];
  G4double perem;
  G4double tampi;
  G4int icoun[10000]; // was [20][50][10] 
  G4int n2;
  G4int n5;
  G4double cm;
  G4double r;
  G4double cn2;
  G4double cp0;
  G4double cn5;
  G4double hp0[20];
  G4int jp0;
  G4double pp0;
  G4double pslowr[1000]; // was [20][50] ;
  G4int khi;
  G4double amp;
  G4int kps;
  G4double amp2;

  if (istart < 1) goto L10; // calculate jp0, index of incident momentum bin 
 L1:
  plet = log(1000.0 / p0);
  jp0 = static_cast<G4int>(plet / delp0);
  if (p0 <= 1000.0) {
    if (jp0 < n2) ++jp0;
    incp = i;
    if (itype > 2) incp = i + 5;
    switch (j) {
    case 1:  goto L100;
    case 2:  goto L200;
    case 3:  goto L300;
    case 4:  goto L400;
    }
  }
  if (j > 1) return;
  w = 0.0;
  return;
  // initialization 
 L10:
  n2 = 20; // if n2 and n5 are to be increased, be sure to change dimensions 
  n5 = 50;
  cn2 = n2;
  cn5 = n5;
  delp0 = log(400.0) / cn2;
  cp0 = exp(-delp0);
  pp0 = 1000.0 / cp0;
  i1 = n2;
  for (k = 0; k < i1; ++k) {
    pp0 *= cp0; // pp0 is highest p0 of bin 
    delps[k] = pp0 / cn5;
    hp0[k] = pp0;
    i2 = n5;
    for (G4int m = 0; m < i2; ++m) {
      cm =  m;
      pslowr[k + m * 20 - 20] = pp0 - delps[k] * cm;
      if (pslowr[k + m * 20 - 20] < 0.0) {
	pslowr[k + m * 20 - 20] = 0.0;
      }
      for (G4int ink = 0; ink < 10; ++ink) icoun[k + (m + ink * 50) * 20 - 1020] = 0;
    } 
  }
  istart = 10;
  amp = 0.94;
  ampi = 0.139;
  tamp = sqr(amp);
  tampi = sqr(ampi);
  amp2 = sqr(amp);
  ampi2 = sqr(ampi);
  goto L1;
 L100:
  if (i > 2) goto L102;
  if (itype > 2) goto L102;
  perem = sqrt(erem * (erem + tamp));
  goto L104;
 L102:
  perem = sqrt(erem * (erem + tampi)); // search primary table for momentum less than perem 
 L104:
  // compute high ps index for starting search, pslowr is descending 
  khi = static_cast<G4int>((hp0[jp0 - 1] - perem) / delps[jp0 - 1]); 
  if (khi < n5) ++khi;
  i1 = n5;
  for (k = khi; k < i1; ++k) {
    if (icoun[jp0 + (k + incp * 50) * 20 - 1020] > 0) goto L108;
  }
  w = 0.0;
  return;
 L108:
  w = 1.0;
  --icoun[jp0 + (k + incp * 50) * 20 - 1021];
  r = G4UniformRand();
  ps = pslowr[jp0 + k * 20 - 21] + delps[jp0 - 1] * r;
  if (ps < perem) goto L110;
  r = G4UniformRand();
  ps = pslowr[jp0 + k * 20 - 21] + delps[jp0 - 1] * r;
  if (ps < perem) goto L110;
  r = G4UniformRand();
  ps = pslowr[jp0 + k * 20 - 21] + delps[jp0 - 1] * r;
  if (ps < perem) goto L110;
  ps = pslowr[jp0 + k * 20 - 21];
  if (ps == 0.0) ps = perem * r;
 L110:
  if (i > 2) goto L112;
  e = sqrt(sqr(ps) + amp2) - amp;
  wmass = 0.0;
  if (itype > 2) goto L114;
  prot(p0, e, ps, a, b, g);
  return;
 L112:
  e = sqrt(sqr(ps) + ampi2) - ampi;
  wmass = ampi;
 L114:
  // pi(p0, i, e, ps, a, b, g); //::: defined in G4BertiniCascade
  return;
 L200:
  kps = static_cast<G4int>((hp0[jp0 - 1] - ps) / delps[jp0 - 1]); // score primary 
  if (kps < n5) ++kps;
  ++icoun[jp0 + (kps + incp * 50) * 20 - 1021];  // note: ++
  return;
 L300:
  kps = static_cast<G4int>((hp0[jp0 - 1] - ps) / delps[jp0 - 1]); // score secondary 
  if (kps < n5) ++kps;
  --icoun[jp0 + (kps + incp * 50) * 20 - 1021]; // note: --
  return;
 L400:
  kps = static_cast<G4int>((hp0[jp0 - 1] - ps) / delps[jp0 - 1]);  // search secondary table for momentum, ps 
  if (kps < n5) ++kps;
  if (icoun[jp0 + (kps + incp * 50) * 20 - 1021] < 0) goto L420;
  return;
 L420:
  ++icoun[jp0 + (kps + incp * 50) * 20 - 1021];
  i = 0;
}

void G4BertiniInelasticCollision::alp19() {
  // description:
  pt[1]  = 1.0;
  pt[25] = 1.0;
  pt[13] = 3.0;
  pt[15] = massPionZero;
  G4double r = G4UniformRand();
  if (isw[11] > 0) goto L70;
  if (r > 0.25) goto L55;
  if (isw[3] == 0) pt[1] = 2.0;
  if (G4UniformRand() > twoThirds) goto L40;
  pt[13] = 4.0;
 L26:
  if (isw[3] != 0) return;
  pt[25] = 2.0;
  return;
 L40:
  pt[15] = massPionCharged;
  if (isw[3] != 0) {
    pt[25] = 2.0;
    return;  
  }
 L50:
  pt[13] = 5.0;
  return;
 L55:
  if (isw[3] != 0) goto L65;
  pt[25] = 2.0;
  goto L40;
 L65:
  pt[1] = 2.0;
  pt[15] = massPionCharged;
  return;
 L70:
  if (r > 0.5) goto L95;
  if (isw[3] != 0) pt[1] = 2.0;
  r = G4UniformRand(); 
  pt[13] = 4.0;
  goto L26;
 L95:
  if (isw[3] != 0) goto L105;
  pt[1] = 2.0;
 L105:
  if (G4UniformRand() > twoThirds) goto L115;
  pt[13] = 4.0;
  if (isw[3] != 0){
    pt[25] = 2.0;
  } 
  return;
 L115:
  pt[15] = massPionCharged;
  if (isw[3] != 0) goto L50;
  pt[25] = 2.0;
  return;
}

void G4BertiniInelasticCollision::alp28() {
  G4double r = G4UniformRand();
  if (isw[12] != 0) goto L95;
  if (r > 0.6) goto L50;
  pt[3] = massPionCharged;
  r = G4UniformRand();
  if (isw[3] == 0) goto L35;
  if (r > oneThird) goto L30;
 L20:
  pt[25] = 5.0;
 L21:
  pt[27] = massPionCharged;
 L25:
  return;
 L30:
  pt[25] = 4.0;
 L31:
  pt[37] = 2.0;
  goto L25;
 L35:
  pt[1] = 5.0;
  pt[13] = 2.0;
  if (r <= oneThird) goto L40;
  else goto L45;
 L40:
  pt[27] = massPionCharged;
  goto L31;
 L45:
  pt[25] = 4.0;
  return;
 L50:
  r = G4UniformRand();
  if (isw[3] == 0) goto L75;
  if (r <= twoThirds) goto L60;
  else goto L70;
 L60:
  pt[1] = 4.0;
 L65:
  r = G4UniformRand();
  if (r <= twoThirds) goto L45;
  else goto L40;
 L70:
  pt[13] = 2.0;
 L71:
  pt[3] = massPionCharged;
  goto L65;
 L75:
  if (r > twoThirds) goto L90;
  pt[1] = 4.0;
 L81:
  pt[13] = 2.0;
 L85:
  r = G4UniformRand();
  if (r <= twoThirds) goto L30;
  else goto L20;
 L90:
  pt[1] = 5.0;
  pt[3] = massPionCharged;
  goto L85;
 L95:
  if (r > value1) goto L115;
  pt[3] = massPionCharged;
  if (isw[3] != 0) goto L110;
  pt[1] = 5.0;
  pt[13] = 2.0;
  goto L21;
 L110:
  pt[37] = 2.0;
  goto L20;
 L115:
  r = G4UniformRand();
  if (isw[3] == 0) goto L135;
  if (r > oneThird) goto L130;
  pt[3] = massPionCharged;
  goto L81;
 L130:
  pt[1] = 4.0;
  goto L85;
 L135:
  if (r > oneThird) goto L145;
  pt[1] = 5.0;
  goto L71;
 L145:
  pt[13] = 2.0;
  goto L60;
}

void G4BertiniInelasticCollision::alpha() {
  // description:
  G4double r = G4UniformRand();
  if (value3 < 0.0) goto L130;
  if (value3 > 0.0) goto L60;
  if (r > value1) goto L50;
  if (isw[10] != 0) goto L15;
  pt[1] = 5.0;
  pt[25] = 2.0;
 L15:
  pt[3] = massPionCharged;
  massParticle[3] = massPionCharged;
  r = G4UniformRand();
  if (r > value2) goto L30;
  pt[13] = 4.0;
  return;
 L30:
  if (isw[10] != 0) goto L45;
  pt[25] = 1.0;
 L36:
  pt[13] = 5.0;
 L40:
  pt[15] = massPionCharged;
  return;
 L45:
  pt[25] = 2.0;
  goto L40;
 L50:
  pt[1] = 4.0;
  if (isw[10] != 0) goto L40;
  pt[13] = 5.0;
  goto L45;
 L60:
  if (r > value1) goto L80;
  massParticle[3] = massPionCharged;
  if (isw[10] == 0) goto L75;
 L65:
  pt[1] = 5.0;
 L66:
  pt[15] = massPionCharged;
  pt[3]  = massPionCharged;
  return;
 L75:
  pt[13] = 5.0;
  pt[25] = 2.0;
  goto L66;
 L80:
  if (r > value2) goto L105;
  pt[1] = 4.0;
  r = G4UniformRand();
  if (r <= value3) goto L100;
  if (isw[10] != 0){
    pt[13] = 4.0;
    return;
  }
 L95:
  pt[25] = 2.0;
  pt[13] = 4.0;
  return;
 L100:
  if (isw[10] != 0) goto L45;
  else goto L36;
 L105:
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  r = G4UniformRand();
  if (r > twoThirds) goto L120;
  if (isw[10] != 0) goto L95;
  pt[1] = 5.0;
  pt[13] = 4.0;
  return;
 L120:
  if (isw[10] != 0) goto L36;
  pt[25] = 2.0;
  goto L65;
 L130:
  if (r > value1) goto L150;
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  r = G4UniformRand();
  if (value3 != -1.0) goto L145;
  if (r <= oneThird) goto L36;
  else goto L95;
 L145:
  pt[1] = 5.0;
  if (r <= oneThird) goto L45;
  else {
    pt[13] = 4.0;
    return; 
  }
 L150:
  if (r > value2) goto L170;
  pt[1] = 4.0;
  r = G4UniformRand();
  if (value3 != -1.0) goto L165;
  if (r <= twoThirds) {
    pt[13] = 4.0;
    return; 
  }
  else goto L45;
 L165:
  if (r <= twoThirds) goto L95;
  else goto L36;
 L170:
  massParticle[3] = massPionCharged;
  pt[3] = massPionCharged;
  if (value3 == -1.0) goto L65;
  pt[13] = 5.0;
  goto L45;
}

void G4BertiniInelasticCollision::angid() {
  // description: calculate  cos(theta) and sin(theta), sin(phi) and cos(phi)
  // uses: curr
  G4double r;
  G4int icurr;
  G4double tesiso;
  icurr = static_cast<G4int>(curr[0]);
  switch (icurr) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L30;
  case 4:  goto L30;
  case 5:  goto L30;
  }
 L10:
  // incident particle - nucleon 
  if (it == 21 || it == 22) goto L20; // single production 
  if (rlke > 3500.0) G4Exception("angid1");
  if (rlke < 500.0) goto L70;
  tesiso = 0.75;
  if (rlke < 1000.0) goto L60;
  tesiso = 0.5;
  if (rlke < 1300.0) goto L60;
  tesiso = 0.25;
  if (rlke < 2500.0) goto L60;
  goto L65; // double production 
 L20:
  if (rlke > 3500.0) G4Exception("angid2");
  goto L65;
 L30:
  // incident particle - pion 
  r = G4UniformRand();
  if (rlke > 2500.0) G4Exception("angid3");
  cst = -0.9999995;
  snt = 0.003162;
  if (it != 11) goto L40;
  if (r <= 0.75) goto L70;
  goto L80;
  // (pi+, p) (pi-, n) (pi0, n) (pi0, p) 
 L40:
  if (it != 12 && it != 28) G4Exception("angid4");
  if (rlke < 500.0) cst = -cst; // opposite didection
  if (r <= 0.8) goto L70;
  goto L80;
 L60:
  r = G4UniformRand();
  if (r <= tesiso) goto L70; // backward / forward 
 L65:
  r = G4UniformRand(); // test for direction 
  cst = 0.9999995;
  snt = 0.003162;
  if (r <= 0.5) goto L80;
  cst = -0.9999995;
  goto L80; // isotropic 
 L70:
  pol1(cst, snt); // calculate cos(phi), sin(phi) 
 L80:
  // azio(sopc, sops); //:::
  return;
} 

void G4BertiniInelasticCollision::bg6c(G4int int1) {
  // description:
  if (knot >= 7) goto L81; 
  abz = 0.0;
  if (knot >= 2) goto L15;
  if (int1 != 0) goto L40;
  else goto L50;
 L15:
  if (knot > 5) goto L30;
  if (knot == 5) goto L25;
  if (int1 != 0) goto L50;
  else goto L40;
 L25:
  it = 11;
  if (int1 != 0) goto L45;
  else goto L55;
 L30:
  it = 12;
  if (int1 != 0) goto L55;
  else goto L45;
 L40:
  it = 2 * knot - 1;
 L45:
  strkp = -1.0;
 L46:
  i1 = 0;
  goto L60;
 L50:
  it = 2 * knot;
 L55:
  strkp = -2.0;
 L56:
  i1 = 1;
 L60:
  i2 = static_cast<G4int>(clsm);
  bb();
 L65:
  isom();
  switch (knot) {
  case 1:  goto L70;
  case 2:  goto L70;
  case 3:  goto L70;
  case 4:  goto L74;
  case 5:  goto L70;
  case 6:  goto L70;
  case 7:  goto L80;
  case 8:  goto L80;
  case 9:  goto L80;
  case 10:  goto L80;
  case 11:  goto L80;
  case 12:  goto L80;
  case 13:  goto L70;
  case 14:  goto L70;
  case 15:  goto L70;
  case 16:  goto L70;
  case 17:  goto L74;
  case 18:  goto L70;
  case 19:  goto L70;
  }
 L70:
  if (rlke <= 2500.0) goto L75;
  else goto L65;
 L74:
  rlke = 0.0;
 L75:
  return;
 L80:
  if (rlke <= 3500.0) goto L75;
  else goto L65;
 L81:
  if (knot > 12) goto L135;
  if (in != 0) goto L155;
  if (knot > 8) goto L136;
  if (knot == 8) goto L115;
  if (int1 == 0) goto L105;
  it = 2 * (knot + 1);
 L100:
  strkp = -2.0;
  goto L46;
 L105:
  it = 2 * knot + 1;
 L110:
  strkp = -1.0;
  goto L56;
 L115:
  if (int1 != 0) {
    it = 2 * (knot + 1); 
    goto L45;
  }
  it = 2 * knot  + 1;
  goto L55;
 L135:
  if (knot >= 18) goto L195;
 L136:
  it = knot + 10; 
  if (knot > 10) goto L150;
  if (knot == 10) goto L145;
 L140:
  if (int1 != 0) goto L45;
  else goto L55;
 L145:
  if (int1 != 0) goto L100;
  else goto L110;
 L150:
  if (knot < 12) {
    goto L140;
  } else if (knot == 12) {
    goto L145;
  } else {
    goto L185;
  }
 L155:
  if (knot > 8) goto L136;
  if (knot == 8) goto L170;
  if (int1 == 0) goto L136;
  it = knot + 11;
  goto L45;
 L170:
  if (int1 != 0) goto L180;
  it = 2 * knot - 1;
  goto L110;
 L180:
  it = knot + knot;
  goto L100;
 L185:
  abz = 0.0;
  i1 = knot - 15;
  if (i1 < 0) {
    goto L45;
  } else if (i1 == 0) {
    goto L55;
  } 
  if (knot != 18) goto L55;
  else goto L45;
 L195:
  it = 28;
  goto L185;
}

G4int G4BertiniInelasticCollision::coll(G4int m) {
  // description:
  // uses:
  // changes:
  if (m < 0) a = sqr(massParticle[3]);
  else a = sqnm;
  col[14] = 0.0;
  medium = static_cast<G4int>(clsm);
  eco[0] = cfepn[medium - 1]; // proton  energy cut-off 
  eco[1] = cfepn[medium + 2]; // neutron energy cut-off 
  col[0] = energy[0] + energy[1]; // total energy particles 1 and 2 
  for (G4int i = 0; i < 3; ++i) col[i + 1] = sqr(massParticle[i]); 
  col[4] = col[2] + col[1] + (energy[0] * energy[1] - (pxyz[0] * pxyz[1] + pxyz[4] * pxyz[5] + pxyz[8] * pxyz[9])) * 2.0;
  col[5] = sqrt(col[4]);
  col[6] = col[5] / col[0]; // gamma 
  col[7] = col[5] * 2.0;
  col[8] = (col[3] + col[4] - a) / col[7];
  com2 = col[8] * col[8];
  if (col[3] > 2.9882156e27) goto L52; // pm(3) = isobar-- //::: fix to G4 constant 
  if (col[3] >= 2.9872156e27) goto L552; // test for roundoff range, (min) sqd + or - 5d23  // ::: compare to original
  else goto L53; // pion or nucleon mass = pm(3) 
 L552:
  col[3] = 2.9877156e27;
  massParticle[2] = 5.466005e13;
 L52:
  if (com2 >= col[3]) goto L56;
  else goto L54;
 L53:
  if (col[3] <= sqnm) goto L52; // have N or pi-
  else goto L11; // go to error 
 L54:
  if (com2 < (col[3] * 0.99)) goto L56;
  com2 = col[3];
  col[8] = massParticle[2];
 L56:
  col[9] = sqrt(com2 - col[3]); // p3 prime 
  col[10] = (col[4] + col[1] - col[2]) / col[7]; // e1 prime 
  col[11] = sqrt(sqr(col[10]) - col[1]); // p1 prime 
  col[12] = (col[6] * energy[0] - col[10]) / col[11]; // beta 
  com = 1.0 - (col[12] * col[12] + col[6] * col[6]);
  if (com >= 5e-06) goto L25;
  if (com >= -5e-06) goto L20;
 L11:
  col[14] = 1.0; // error 
 L12:
  return 0;
 L20:
  col[13] = 0.002236067977;
  goto L30;
 L25:
  col[13] = sqrt(com);
  // alpha 
 L30:
  energy[2] = (col[8] + col[9] * (col[12] * cst + col[13] * sopc * snt)) / col[6];
  energy[3] = col[0] - energy[2];
  goto L12;
}

G4int G4BertiniInelasticCollision::collm(G4int m) {
  // description:
  // uses:
  // changes:
  G4double b;
  G4int k;

  univ = energy[1] + col[5] - col[10];
  unive = energy[0] + col[10];
  univer = col[0] + col[5];
  k = 16;
  for (G4int i = 0; i < 9; i += 4) {
    col[k] = (pxyz[i] * univ - pxyz[i + 1] * unive) / univer;
    col[k + 3] = (pxyz[i] + pxyz[i + 1]) / col[0]; // vx 
    ++k;
  }
  col[21] = (pxyz[9] * pxyz[4] - pxyz[8] * pxyz[5]) / col[0]; // qx 
  col[22] = (pxyz[1] * pxyz[8] - pxyz[9] * pxyz[0]) / col[0]; // qy 
  col[23] = (pxyz[5] * pxyz[0] - pxyz[4] * pxyz[1]) / col[0];
  a = snt / col[13];
  b = a * col[9]; // (-beta * cos phi * sin theta / alpha + cos theta) / p1p * p3p 
  univ = col[9] * (cst - a * sopc * col[12]) / col[11];
  unive = b * sops / col[11]; // p3p * sin phi * sin theta / p1p * alpha 
  univer = sopc * b + (energy[2] + col[8]) / (col[6] + 1.0); // cos phi * sin theta * p3p / alpha  +  (e3 + e3p) / (1.0 + gamma) 
  k = 19;
  for (i = 1; i < 10; i += 4) {   
    pxyz[i] = col[k] * univer + col[k + 3] * unive + col[k - 3] * univ;  // :: index fixed
    ++k;
  }
  if (m == 0) goto L15;
  if (pt[14] == 0.0) goto L25;
 L15:
  for (i = 0; i < 9; i += 4) pxyz[i + 2] = pxyz[i - 1] + pxyz[i] - pxyz[i + 1];
  if (m == 0) goto L60;
 L25:
  if (pt[2] == 0.0) goto L45;
  pt[3] = massParticle[2];
  i1 = 3;
 L35:
  i2 = -1;
  pstor();
  if (i1 > 3) goto L55;
 L45:
  if (pt[14] == 0.0) goto L55;
  pt[15] = massNucleon;
  i1 = 4;
  goto L35;
 L55:
  pt[26] = 0.0;
  pt[38] = 0.0;
 L60:
  return 0;
}

void G4BertiniInelasticCollision::ecpl() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  i1 = 0;
  medium = static_cast<G4int>(clsm);
  if (pt[37] != 0.0) goto L240;
  if (curr[0] < 1.0) {
    goto L135;
  } else if (curr[0] == 1.0) {
    goto L215;
  }
  if (curr[0]< 3.0) {
    goto L170;
  } else if (curr[0] == 3.0) {
    goto L155;
  }
  if (curr[0] < 5.0) {
    goto L110;
  } else if (curr[0] > 5.0) goto L135;
  if (strkp == -1.0) goto L30;
  else goto L135;
  if (strkp != -2.0) goto L135;
  else goto L90;
 L30:
  if (pt[1] != 2.0) goto L32;
  pt[2] = vnvp[medium - 1];
  goto L33;
  // pi-, curr(1) = 5.0, pt6 + 1 = 0, strkp = 1.0 
 L32:
  pt[2] = 0.0;
 L33:
  pt[14] = hvp[medium - 1];
  if (pt[25] < 1.0) {
    goto L135;
  } else if (pt[25] == 1.0) {
    goto L75;
  }
  if (pt[25] != 2.0) goto L135;
  pt[26] = -ppan[medium - 1];
  if (curr[0] >= 2.0) goto L41;
  else goto L235;
  // ppan = - vnhp(neutron well depth - 1 / 2 proton well depth) 
 L41:
  if (curr[0] < 4.0) {
    goto L165;
  } else if (curr[0] == 4.0) {
    goto L135;
  }
 L45:
  if (pt[1] < 3.0) {
    goto L135;
  } else if (pt[1] == 3.0) {
    goto L55;
  }
  if (pt[1] < 5.0) {
    goto L65;
  } else if (pt[1] == 5.0) {
    goto L70;
  } else {
    goto L135;
  }
 L55:
  if (pt[13] != 5.0) goto L135;
 L60:
  return; 
 L65:
  if (pt[13] != 4.0) goto L135;
  else goto L60;
 L70:
  if (pt[13] != 3.0) goto L135;
  else goto L60; // 1 / 2 proton well depth 
 L75:
  pt[26] = hvp[medium - 1];
  if (curr[0] >= 2.0) goto L76;
  else goto L225;
 L76:
  if (curr[0] < 4.0) {
    goto L145;
  } else if (curr[0] == 4.0) {
    goto L135;
  } 
 L80:
  if (pt[1] < 4.0) {
    goto L135;
  } else if (pt[1] == 4.0) {
    goto L55;
  }
  if (pt[1] != 5.0) {
    goto L135;
  } else {
    goto L65;
  }
 L90:
  pt[2] = -vnvp[medium - 1];
 L91:
  pt[14] = -pmac[medium - 1];
 L92:
  if (pt[25] < 1.0) {
    goto L135;
  } else if (pt[25] == 1.0) {
    goto L100;
  }
  if (pt[25] != 2.0) {
    goto L135;
  } else {
    goto L105;
  }
 L100:
  pt[26] = -pmac[medium - 1];
  if (curr[0] <= 2.0) goto L55;
  if (curr[0] < 4.0) {
    goto L45;
  } else if (curr[0] == 4.0) {
    goto L80;
  }
  if (pt[1] != 5.0) {
    goto L135;
  } else {
    goto L55;
  }
 L105:
  pt[26] = hvn[medium - 1];
  if (curr[0] <= 2.0) goto L65;
  if (curr[0] < 4.0) {
    goto L145;
  } else if (curr[0] == 4.0) {
    goto L45;
  } else {
    goto L80;
  }
 L110:
  if (strkp == 1.0) goto L120;
  if (strkp > 1.0) goto L135;
  if (strkp != -2.0) goto L135;
  else goto L90;
 L120:
  pt[2] = 0.0; // pi0 
 L121:
  pt[14] = hvp[medium - 1];
  if (pt[25] < 1.0) {
    goto L135;
  } else if (pt[25] == 1.0) {
    goto L130;
  }
  if (pt[25] != 2.0) goto L135;
  else goto L140;
 L130:
  pt[26] = hvp[medium - 1];
  if (curr[0] < 4.0) {
    goto L65;
  } else if (curr[0] == 4.0) {
    goto L45;
  }
 L135:
  i1 = 1;
  goto L60;
 L140:
  pt[26] = -ppan[medium - 1];
  if (curr[0] < 4.0) goto L70;
  if (curr[0] > 4.0) goto L135;
 L145:
  if (pt[1] < 3.0) {
    goto L135;
  } else if (pt[1] == 3.0) {
    goto L65;
  }
  if (pt[1] != 4.0) goto L135;
  else goto L70;
 L155:
  if (strkp > -1.0) goto L135;
  if (strkp == -1.0) goto L30;
  if (strkp != 2.0) goto L135;   // pi+ 
  else goto L90;
 L165:
  if (pt[1] != 3.0) goto L135;
  else goto L70;
 L170:
  if (strkp > -1.0) goto L135;
  if (strkp == -1.0) goto L180;
  if (strkp != -2.0) goto L135; // neutron 
  else goto L190;
 L180:
  pt[2] = 0.0;
  if (pt[1] < 1.0) {
    goto L135;
  } else if (pt[1] == 1.0) {
    goto L91;
  }
  if (pt[1] != 2.0) goto L135;
  else goto L121;
 L190:
  pt[14] = -pmac[medium - 1];
  if (pt[1] < 1.0) {
    goto L135;
  } else if (pt[1] == 1.0) {
    goto L200;
  }
  if (pt[1] != 2.0) goto L135;
  else goto L210;
 L200:
  if (pt[25] != 2.0) goto L135;
  pt[2] = -vnvp[medium - 1];
  pt[26] = hvn[medium - 1];
  goto L55;
 L210:
  pt[2] = 0.0;
  goto L92;
 L215:
  if (strkp > -1.0) goto L135;
  if (strkp == -1.0) goto L30;
  if (strkp != -2.0) goto L135;
  else goto L180;
 L225:
  if (pt[1] < 1.0) { // proton 
    goto L135;
  } else if (pt[1] == 1.0) {
    goto L65;
  }
  if (pt[1] != 2.0) goto L135;
  else goto L70;
 L235:
  if (pt[1] != 1.0) goto L135;
  else goto L70;
 L240:
  if (curr[0] < 1.0) {
    goto L135;
  } else if (curr[0] == 1.0) {
    goto L250;
  }
  if (curr[0] != 2.0) goto L135;
  else goto L260;
 L250:
  if (strkp > -1.0) goto L135;
  if (strkp == -1.0) goto L360;
  if (strkp != 2.0) goto L135;
  else goto L270;
 L260:
  if (strkp > -1.0) goto L135;
  if (strkp == -1.0) goto L270;
  if (strkp != 2.0) goto L135;
  else goto L405;
 L270:
  if (pt[13] <= 0.0) goto L135;
  if (pt[37] <= 0.0) goto L135;
  if (pt[37] > 2.0) goto L135;
  if (pt[37] == 2.0) goto L314;
  if (pt[13] > 2.0) goto L135;
  if (pt[13] == 2.0) goto L340;
  pt[2] = tffn[medium - 1];
  pt[14] = tffn[medium - 1];
  pt[26] = tffn[medium - 1];
  pt[38] = tffn[medium - 1];
 L295:
  if (pt[1] <= 3.0) goto L135;
  if (pt[1] < 5.0) goto L310;
  if (pt[1] > 5.0) goto L135;
 L305:
  if (pt[25] != 4.0) goto L135;
  else goto L60;
 L310:
  if (pt[25] != 5.0) goto L135;
  else goto L60;
 L314:
  if (pt[13] < 2.0) {
    goto L315;
  } else if (pt[13] == 2.0) {
    goto L345;
  }
 L315:
  pt[26] = ffptfn[medium - 1];
  pt[2] = fvnp[medium - 1];
 L316:
  pt[38] = fvnp[medium - 1];
  pt[14] = fvnp[medium - 1];
 L317:
  if (pt[1] < 3.0) {
    goto L135;
  } else if (pt[1] == 3.0) {
    goto L310;
    if (pt[1] < 5.0) goto L305;
    if (pt[1] > 5.0) goto L135;
  L325:
    if (pt[25] != 3.0) goto L135;
    else goto L60;
  L340:
    pt[2] = ffptfn[medium - 1];
    pt[26] = fvnp[medium - 1];
    goto L316;
  L345:
    pt[2] = tffn[medium - 1];
    pt[26] = tffn[medium - 1];
    pt[14] = tffp[medium - 1];
    pt[38] = tffp[medium - 1];
  L350:
    if (pt[1] < 3.0) {
      goto L135;
    } else if (pt[1] == 3.0) {
      goto L305;
    }
    if (pt[1] != 4.0) {
      goto L135;
    } else {
      goto L325;
    }
  L360:
    if (pt[13] <= 0.0) goto L135;
    if (pt[37] < 1.0) {
      goto L135;
    } else if (pt[37] == 1.0) {
      goto L375;
    }
    if (pt[37] != 2.0) {
      goto L135;
    } else {
      goto L390;
    }
  L375:
    if (pt[13] > 2.0) goto L135;
    if (pt[13] == 2.0) goto L385;
    pt[2] = hvp[medium - 1];
  L381:
    pt[14] = hvp[medium - 1];
    pt[26] = hvp[medium - 1];
    pt[38] = hvp[medium - 1];
    if (pt[13] > 2.0) {
      goto L317;
    } else if (pt[13] == 2.0) {
      goto L295;
    } else {
      goto L135;
    }
  L385:
    pt[26] = hvn[medium - 1];
  L386:
    pt[2] = -pmac[medium - 1];
    pt[38] = hvn[medium - 1];
    pt[14] = hvn[medium - 1];
    if (pt[37] < 2.0) {
      goto L350;
    } else if (pt[37] == 2.0) {
      goto L317;
    } else {
      goto L135;
    }
  L390:
    if (pt[13] > 2.0) goto L135;
    if (pt[13] == 2.0)goto L400;
    pt[2] = hvn[medium - 1];
    pt[14] = hvn[medium - 1];
    pt[26] = -pmac[medium - 1];
    pt[38] = hvn[medium - 1];
    goto L350;
  L400:
    pt[14] = -ppan[medium - 1];
    pt[38] = -ppan[medium - 1];
    pt[2] = hvp[medium - 1];
    pt[26] = hvp[medium - 1];
    if (pt[1] != 3.0) {
      goto L135;
    } else {
      goto L325;
    }
  L405:
    if (pt[13] <= 0.0) goto L135;
    if (pt[37] < 1.0) {
      goto L135;
    } else if (pt[37] == 1.0) {
      goto L420;
    }
    if (pt[37] != 2.0) goto L135;
    else goto L435;
  L420:
    if (pt[13] > 2.0) goto L135;
    if (pt[13] == 2.0) goto L430;
    pt[38] = -pmac[medium - 1];
    pt[14] = -pmac[medium - 1];
    pt[26] = -pmac[medium - 1];
    pt[2] = -pmac[medium - 1];
    if (pt[1] != 5.0) {
      goto L135;
    } else {
      goto L310;
    }
  L430:
    pt[2] = thpn[medium - 1];
    goto L381;
  L435:
    if (pt[13] > 2.0) goto L135;
    if (pt[13] == 2.0) goto L445;
    pt[2] = hvp[medium - 1];
    pt[14] = hvp[medium - 1];
    pt[38] = hvp[medium - 1];
    pt[26] = thpn[medium - 1];
    goto L295;
  L445:
    pt[26] = -pmac[medium - 1];
    goto L386;
  }
}

void G4BertiniInelasticCollision::geo() {
  // decription: 
  // parameters:
  // uses: coordinate,
  // changes:
  G4double temp;
  G4double tempo;
  G4double t5; 
  G4double t6;
  i1 = 0;
  G4double t1 = out[1] * out[1]; // sqr(r1    ) 
  G4double t2 = out[2] * out[2]; // sqr(r1 + 1)  
  G4double t3 = out[3] * out[3]; // sqr(r1 + 2)
  G4double t4 = t3 * 2.0;        // sqr(r1 + 2) * 2 
 L95:
  t5 = sqr(coordinate[0]) + sqr(coordinate[1]) + sqr(coordinate[2]); // t5 = sqr(r) 
  switch (medium) {
  case 1:  goto L100;
  case 2:  goto L125;
  case 3:  goto L155;
  case 4:  goto L185;
  }
 L100:
  t6 = t5 - t1;
  if (t6 <= 0.0) goto L210; // medium = 1
  temp = t1;
 L110:
  if ((t6 / temp) > 5e-06) goto L205;
  G4int i;
  for (i = 0; i < 3; ++i) {
    coordinate[i] *= 0.999995; 
    curr[i + 1] = coordinate[i];
  }
  goto L95;
 L125:
  t6 = t5 - t1; // medium = 2 
  if (t6 > 0.0) goto L150;
  temp = t1;
 L135:
  if ((t6 / temp) < -5e-06) goto L205;
  for (i = 0; i < 3; ++i) {
    coordinate[i] *= 1.000005; 
    curr[i + 1] = coordinate[i];
  }
  goto L95;
 L150:
  t6 = t5 - t2;
  temp = t2;
  if (t6 <= 0.0) goto L210;
  else goto L110;
 L155:
  t6 = t5 - t2;
  // medium = 3 
  if (t6 > 0.0) goto L165;
  temp = t2;
  goto L135;
 L165:
  t6 = t5 - t3;
  if (t6 <= 0.0) goto L210;
  temp = t3; // dummy if st. follows to keep st. 175 
  if (temp != t3) goto L175;
  goto L110;
 L175:
  if (coordinate[1] != 0.0) goto L205; // medium = 4
  if (curr[4] != 0.0) goto L205;
 L185:
  t6 = t5 - t3;
  if (t6 > 0.0) goto L195;
  temp = t3;
  goto L135;
 L195:
  t6 = t5 - t4;
  if (t6 <= 0.0) goto L210;
  temp = t4;
  goto L110;
 L205:
  i1 = -1;
  goto L10;
 L210:
  t4 = coordinate[0] * dcos[0] + coordinate[1] * dcos[1] + coordinate[2] * dcos[2]; // t4 = -b = -rcos(theta) = sum of xi(i) * dcos(i) 
  t6 = sqr(t4); 
  t6 = t5 - t6; // t6 = sqr(r) - sqr(b) 
  if (t3 < t6) goto L205;
  t3 = sqrt(t3 - t6); // sqrt(sqr(b)-sqr(r) + sqr(radius3)) similar for t2 = a2 and t1 = a1 
  temp = t2 - t6;
  t2 = sqrt(fabs(temp));
  tempo = t1 - t6;
  t1 = sqrt(fabs(tempo));
  for (i = 0; i < 6; ++i) d[i] = 0.0;
  switch (medium) {
  case 1:  goto L5;
  case 2:  goto L15;
  case 3:  goto L39;
  case 4:  goto L65;
  }
 L5:
  if (temp < 0.0) goto L205;
  d[3] = t1 - t4; // b + a1 
  d[4] = t2 - t1; // a2 - a1 
  d[5] = t3 - t2; // a3 - a2 
 L10:
  return;
 L15:
  if (temp < 0.0) goto L205;
  d[5] = t3 - t2; 
  if (t4 < 0.0) goto L25;
 L20:
  d[2] = t2 - t4; // b + a2 
  goto L10;
 L25:
  if (tempo < 0.0) goto L20;
  d[2] = - (t4 + t1); // b - a1 
  d[3] = 2 * t1;      // 2 * a1 
  d[4] = t2 - t1;     // a2 - a1 
  goto L10;
 L39:
  if (t4 < 0.0) goto L45;
 L40:
  d[1] = t3 - t4;     // b + a3 
  goto L10;
 L45:
  if (temp < 0.0) goto L40;
  d[1] = - (t4 + t2); // b - a2 
  d[5] = t3 - t2;     // a3 - a2 
  if (tempo >= 0.0) goto L60;
  d[2] = 2 * t2;      // 2  * a2 
  goto L10;
 L60:
  d[2] = t2 - t1;     // a2 - a1 
  d[4] = d[2];
  d[3] = 2 * t1;      // 2 * a1 
  goto L10;
 L65:
  d[0] = - (t4 + t3); 
  if (temp >= 0.0) goto L75;
  d[1] = t3 + t3;
  goto L10;
 L75:
  d[1] = t3 - t2;
  d[5] = d[1];        // b - a3, a3 - a2, region 4 
  if (tempo < 0.0) goto L85;
  d[2] = t2 - t1;     // a2 - a1 
  d[4] = d[2];
  d[3] = t1 + t1;
  goto L10;
 L85:
  d[2] = 2 * t2;      // 2 * a2 
  goto L10;
} 

void G4BertiniInelasticCollision::pstor() {
  // decription: sets pt vector 
  // parameters: - 
  // uses: i1, massNucleon, pxyz, curr, clmsc, c
  // changes: pt
  G4int j;
  G4int k;
  G4int l;
  G4int m;
  G4int jj;
  l = i1 * 12 - 28;
  if (i2 == 0) goto L14;
  else goto L23;
  jj = 0;
  if (massParticle[2] > massNucleon) {
    ++i1;
    jj = 1;
  }
  univ = sqrt(sqr(pxyz[i1 - 1]) + sqr(pxyz[i1 + 3]) + sqr(pxyz[i1 + 7])); // x-y-z -coordinates of collision point 
  k = i1 + 8;
  i1 = k;
  G4int i;
  for (i = i1; i < i1; i += 4) {
    pt[l] = pxyz[i] / univ; 
    ++l;
  }
  i1 -= jj;
 L6:
  pt[l - 1] = clsm;
  pt[l    ] = curr[10];
  pt[l - 7] = c[0];
  pt[l - 6] = c[1];
  pt[l - 5] = c[2];
  return;
 L14:
  k = 14;
  goto L25;
 L23:
  if (i2 >= 2) goto L35;
  k = 17;
 L25:
  univ = sqrt(sqr(pnidk[k - 1]) + sqr(pnidk[k]) + sqr(pnidk[k + 1]));
  j = k + 2;
  i1 = j;
  for (i = k - 1; i < i1; ++i) {
    pt[l - 1] = pnidk[i] / univ;
    ++l;
  }
  goto L6;
 L35:
  univ = sqrt(sqr(pt[l - 4]) + sqr(pt[l - 3]) + sqr(pt[l - 2]));
  k = l - 1;
  m = l - 3;
  i1 = k;
  for (i = m; i <= i1; ++i) {
    pt[l - 1] = pt[i] / univ;
    ++l;
  }
  goto L6;
} 

void G4BertiniInelasticCollision::punp() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  G4int l;
  G4int m;
  G4double ac;
  G4double ar;
  G4double zc;
  G4double ze;
  if (pgvc[0] > 0.0) goto L14; // pgvc is particle with velocity greater than criterion 
  if (pgvc[0] == 0) goto L20; 
 L5:
  i1 = -1;
 L10:
  return;
 L14:
  pgvc[0] += -11.0;
  k = static_cast<G4int>(pgvc[0] + 2.005);
  G4int i;
  for (i = 0; i < 11; ++i) {
    curr[i] = pgvc[k];
    pgvc[k] = 0.0;
    ++k;
  }
 L16:
  i1 = 1;
  goto L10;
 L20:
  if (plvc[0] < 0.0) {
    goto L5;
  } else if (plvc[0] == 0) {
    goto L55;
  }
  univ = 0.0;
  l = static_cast<G4int>(plvc[0]);
  k = -10;
  i1 = l;
  for (i = 0; i < i1; ++i) {
  L30:
    k += 12;
    if (plvc[k] < 0.0) {
      goto L5;
    } else if (plvc[k] == 0) {
      goto L30;
    }
    if (plvc[k] >= univ){
      univ = plvc[k];
      m = k;
    }
  }
  plvc[m] = 0.0;
  for (i = 0; i < 11; ++i) {
    ++m;
    curr[i] = plvc[m];
    plvc[m] = 0.0;
  }
  plvc[0] += -1.0;
  goto L16;
 L55:
  i1 = 0;
  ac = amasno;
  if (no > 2) goto L150;
  zc = zee;
  ac = amasno + 1.0; // ac-compound nuc, zc = chg. compound, ar = mass cascade residual nucleus 
  if (no == 1) goto L160;
  goto L170;
 L150:
  if (no >= 5) {
    zc = zee - 1.0;
    goto L170;
  }
 L160:
  zc = zee + 1.0;
 L170:
  ar = ac - count[0] - count[1];
  ze = count[0] + count[2] - count[4];
  if (ar == 0) goto L190;
  if (ar > 0) goto L200;
  count[6] += 1.0;
  goto L220;
 L190:
  if (zc == ze) goto L10;
  count[7] += 1.0;
  goto L220;
 L200:
  if (ar >= (zc - ze)) goto L210;
  count[8] += 1.0;
  goto L220;
 L210:
  if (zc >= ze) goto L10;
  count[9] += 1.0;
 L220:
  if (begru == 1.0) goto L240;
  begru += -1.0;
 L222:
  for (i = 0; i < 12; ++i) {
    cc[i] = pcc[i];
    ipec[i] = nip[i];
  }
  for (i = 0; i < 5; ++i) pnbc[i] = ppnb[i];
  --nor;
  count[5] = -1.0;
  goto L10;
 L240:
  begru = -1.0;
  goto L222;
  // 1. picks up last 11 items in pgvc and stores 
  // them in curr.  stores zero in those 11 pgvc cells. 
  // group is greater than 1st in all other groups.  stores 
  // 2-12th items in curr--zeroes items 1-12 in plvc group
}

void G4BertiniInelasticCollision::spcn() {
  // decription: set vector fmax 
  // parameters: -
  // uses: curr fmax, i
  // changes:
  fmax[3] = 0.0;
  fmax[4] = 0.0;
  i6 = static_cast<G4int>(curr[0] - 1.95);
  G4double wk = wkrpn[i2 - 1];

  if (i6 == 1 || i6 ==3)   univ = massPionCharged;
  if (i6 == 2)   univ = massPionZero;
 
  unive = bovera(wk, univ);
  value1 = 0.0;
  if (wk > 220.0) goto L55;
  switch (i6) {
  case 1:  goto L30;
  case 2:  goto L45;
  case 3:  goto L30;
  }
 L30:
  fmax[4] = unive * 3.5e-28;  // (pi+, p) s.p. 
  fmax[5] = unive * 2.3e-27;  // (pi-, p) s.p. 
  fmax[0] = unive * 2e-25;    // (pi+, p) sc 
  fmax[2] = unive * 4.51e-26; // (pi-, p) ex 
 L35:
  fmax[1] = unive * 2.3e-26;  // (pi-, p) sc 
  return;
 L45:
  fmax[5] = unive * 1.35e-27; // (pi0, p) s.p. 
  fmax[0] = unive * 8.92e-26; // (pi0, p) sc 
  fmax[1] = unive * 4.51e-26; // (pi0, p) ex 
 L50:
  fmax[2] = fmax[0];          // (pi0, n) sc 
  fmax[3] = fmax[1];          // (pi0, n) ex 
  fmax[6] = fmax[5];          // (pi0, n) s.p. 
  return;
 L55:
  if (wk > 400.0) goto L75;
  switch (i6) {
  case 1:  goto L65;
  case 2:  goto L70;
  case 3:  goto L65;
  }
 L65:
  fmax[4] = unive * 4.5e-27;
  fmax[5] = unive * 2.06e-26;
  fmax[0] = 2.0e-25;
  fmax[2] = 4.51e-26;
  goto L35;
 L70:
  fmax[5] = unive * 1.25e-26;
  fmax[0] = 8.92e-26;
  fmax[1] = 4.51e-26;
  goto L50;
 L75:
  if (wk > 500.0) goto L100;
  switch (i6) {
  case 1:  goto L85;
  case 2:  goto L90;
  case 3:  goto L85;
  }
 L85:
  fmax[0] = 1.13e-25;
  fmax[1] = unive * 2.05e-26;
  fmax[2] = 2.77e-26;
  fmax[4] = unive * 1.17e-26;
  fmax[5] = unive * 2.17e-26;
  return;
 L90:
  fmax[0] = 5.08e-26;
  fmax[1] = 2.77e-26;
  fmax[5] = unive * 1.46e-26;
  goto L50;
 L100:
  value1 = 1.0;
  if (wk > 600.0) goto L120;
  switch (i6) {
  case 1:  goto L110;
  case 2:  goto L115;
  case 3:  goto L110;
  }
 L110:
  fmax[0] = 5.1e-26;
  fmax[1] = unive * 2.47e-26;
  fmax[2] = unive * 1.55e-26;
  fmax[4] = unive * 1.51e-26;
  fmax[5] = unive * 3.04e-26;
  return;
 L115:
  fmax[0] = unive * 2.3e-26;
  fmax[1] = unive * 1.55e-26;
  fmax[5] = unive * 2.26e-26;
  goto L50;
 L120:
  if (wk > 800.0) goto L140;
  switch (i6) {
  case 1:  goto L130;
  case 2:  goto L135;
  case 3:  goto L130;
  }
 L130:
  fmax[0] = 3.3e-26;
  fmax[1] = unive * 2.63e-26;
  fmax[2] = unive * 1.2e-26;
  fmax[4] = unive * 2.01e-26;
  fmax[5] = unive * 3.04e-26;
  return;
 L135:
  fmax[0] = unive * 1.6e-26;
  fmax[1] = unive * 1.2e-26;
  fmax[5] = unive * 2.26e-26;
  goto L50;
 L140:
  switch (i6) {
  case 1:  goto L145;
  case 2:  goto L150;
  case 3:  goto L145;
  }
 L145:
  fmax[0] = unive * 1.93e-26;
  fmax[1] = unive * 2.63e-26;
  fmax[2] = unive * 8.2e-27;
  fmax[4] = unive * 2.33e-26;
  fmax[5] = unive * 3.04e-26;
  return;
 L150:
  fmax[0] = unive * 1.6e-26;
  fmax[1] = unive * 8.2e-27;
  fmax[5] = unive * 2.49e-26;
  goto L50;
}

void G4BertiniInelasticCollision::store() {
  // decription: 
  // parameters:
  // uses: zero(), i3, i2, i5
  // changes:
  G4int k; 
  G4int n;
  G4int l;
  G4int m; 
  G4int id; 
  G4int mm;
  G4double sub1;
  G4double sub2;
  zero();
  n = i1;
  if (i3 != 0) goto L85;
  k = i2 + 6;
  l = i2;
  switch (i5) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L10;
  case 4:  goto L15;
  case 5:  goto L10;
  }
 L15:
  if (i1 == 5) goto L100;
  mm = -6;
 L19:
  l = k;
  --n;
 L10:
  i1 = n;
  G4int i;
  for (i = 0; i < i1; i += 2) {
    id = i; // ::: verify for new indexes (was id = i)
    ce[i + 1] = fmax[i] * 1e30 * space[k];
    ce[i + 2] = fmax[i + 1] * 1e30 * space[l];
    switch (i5) {
    case 1:  continue;
    case 2:  continue;
    case 3:  goto L11;
    case 4:  goto L11;
    case 5:  goto L11;
    }
  L11:
    m = k;
    k = l;
    switch (i5) {
    case 1:  continue;
    case 2:  continue;
    case 3:  goto L30;
    case 4:  goto L40;
    case 5:  goto L35;
    }
  L35:
    if (id <= 2) continue;
    else  goto L36;
  L40:
    k = m + mm;
    mm = -mm;
    l = k;
    continue;
  L36:
    k = i2;
  L30:
    l = m;
  }
  switch (i5) {
  case 1:  goto L26;
  case 2:  goto L26;
  case 3:  goto L26;
  case 4:  goto L95;
  case 5:  goto L26;
  }
 L26:
  sign = 0.0;
  for (i = 1; i < 8; ++i) sign += ce[i];
  switch (i5) {
  case 1:  goto L50;
  case 2:  goto L50;
  case 3:  goto L60;
  case 4:  goto L110;
  case 5:  goto L60;
  }
 L50:
  sign *= 0.999999;
  return;
 L60:
  if (i1 > 4) goto L50;
  if (i3 > 0) goto L75;
  space[i2 + 86] = sign;
  goto L80;
 L75:
  space[i2 + 105] = sign;
 L80:
  fmax[4] = 0.0;
  fmax[5] = 0.0;
  goto L50;
 L85:
  k = i2;
  l = i2 + 6;
  switch (i5) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L10;
  case 4:  goto L90;
  case 5:  goto L10;
  }
 L90:
  mm = 6;
  goto L19;
 L95:
  ce[7] = fmax[6] * 1e30 * space[l - 1];
  goto L26;
 L100:
  sub1 = fmax[0] * 1e30;
  sub2 = fmax[1] * 1e30;
  i1 = n;
  for (i = 0; i < i1; i += 2) {
    ce[i    ] = sub1 * space[k];
    ce[i + 1] = sub2 * space[k];
    k = l;
    l = k + 6;
  }
  ce[n + 1] = 1e30 * space[k] * fmax[n];
  goto L26;
 L110:
  if (i1 == 5) {
    space[i2 + 67] = sign;
    fmax[5] = 0.0;
    fmax[6] = 0.0;
  }
  goto L50;
}

void G4BertiniInelasticCollision::stpr() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int k;
  i1 = 0;
  medium = static_cast<G4int>(clsm);
  for (G4int i = 2; i < 39; i += 12) {
    k = i;
    if (pt[i] == 0.0) return;
    if (pt[k - 1] == 2.0) goto L25;
    pt[k - 2] = pt[k] - space[medium + 9];
    goto L40;
  L25:
    pt[k - 2] = pt[k] - space[medium + 3];
  L30:
    if (pt[k - 2] > 500.0) goto L50;
    stpl(&pt[k - 2]); // velocity less than criterion 
    if (i1 != 0) return;
    else return;
  L40:
    if (pt[k - 2] == 1.0) goto L30;
    pt[k - 2] = massNucleon * pt[k - 2] / pt[k + 1];
    goto L30;
  L50:
    stph(&pt[k - 1]); // velocity greater than criterion 
    if (i1 != 0) return; // ::: verify
    else return;
  }
}

void G4BertiniInelasticCollision::undis() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4double r1;
  G4double r2;

  if (begru < 0.0) {
    goto L45;
  } else if (begru == 0) {
    goto L50;
  } else {
    goto L40;
  }
 L5:
  begru = 1.0;
  G4int i;
  for (i = 0; i < 4; ++i) rands[i] = randi[i];
 L10:
  r1 = 2.0 * G4UniformRand() - 1.0;
  r2 = 2.0 * G4UniformRand() - 1.0;
  if ((sqr(r1) + sqr(r2)) >= 1.0) goto L10;
  coordinate[0] = r1 * out[3]; // x, y and z coordinates, also alpha, beta and gamma direction cosines. 
  coordinate[1] = r2 * out[3];
  x = coordinate[0];
  y = coordinate[1];
  curr[3] = coordinate[0];
  curr[4] = coordinate[1];
  coordinate[2] = -out[3];
  dcos[0] = 0.0;
  dcos[1] = 0.0;
  dcos[2] = 1.0;
  curr[5] = coordinate[2];
  curr[6] = 0.0; 
  curr[7] = 0.0;
  curr[8] = 1.0;
  medium = 4; // no. of geom 
  curr[9] = medium; 
  return;
 L19:
  begru += 1.0;
 L20:
  begru += 1.0;
  if (casesn < begru) goto L30;
  frand = G4UniformRand();
  goto L10;
 L30:
  begru = 0.0;
  for (i = 0; i < 4; ++i) {
    randi[i] = rands[i];
    erand[i] = rands[i];
  }
  return; // final random in erand. run completed 
 L40:
  if (count[5] == 0.0) goto L50;
 L45:
  count[5] = 0.0;
 L47:
  if (begru < 0.0) {
    goto L19;
  } else if (begru == 0) {
    goto L5;
  } else {
    goto L20;
  }
 L50:
  for (i = 0; i < 12; ++i) {
    pcc[i] = cc[i]; 
    nip[i] = ipec[i];
  }
  for (i = 0; i < 5; ++i) ppnb[i] = pnbc[i];
  goto L47;
  // when begru = 0, munpu common is zeroed 
}

void G4BertiniInelasticCollision::xyi(G4int ii, G4int jj, G4int kk) {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int    ll; 
  G4int    nm; 
  G4int    mm;
  G4int    nn;
  G4double w1;
  G4double w2;
  G4double w3;
  G4double w4;
  G4double w5;
  G4double w6;
  w1 = s[ii - 1] * 1e30;
  w2 = s[ii] * 1e30;
  if (abs(iv) == 1) goto L20;
  w3 = space[jj - 1] * 1e30;
  w4 = space[jj    ] * 1e30;
  w5 = space[jj + 1] * 1e30;
  w6 = space[jj + 2] * 1e30;
  ll = ii + 2;
  if (iv != 0) goto L15;
  else goto L10;
 L20:
  w3 =     s[ii + 1] * 1e30;
  w4 =     s[ii + 2] * 1e30;
  w5 = space[jj - 1] * 1e30;
  w6 = space[jj    ] * 1e30;
  ll = ii + 4;
  if (iv < 0) goto L35;
  nm = 7;
  i7 = 1;
 L15:
  mm = 7;
  nn = 1;
  goto L40;
 L35:
  nm = 1;
  i7 = 7;
 L10:
  mm = 1;
  nn = 7;
 L40:
  for (G4int i = 0; i < 3; ++i) {
    s[ll - 1] = w1 * space[mm - 1];
    s[ll    ] = w2 * space[nn - 1];
    if (abs(iv) == 1) {
      s[ll + 1] = w3 * space[nm - 1];
      s[ll + 2] = w4 * space[i7 - 1];
      ll += 2;
      ++nm;
      ++i7;
    }
    ++mm;
    ++nn; 
    ll += 2;
  }
  if (iv > 0) goto L70;
  if (iv == 0) goto L65;
  nm = 7;
  i7 = 1;
  goto L80;
 L65:
  mm = 1;
  nn = 7;
  nm = 7;
  i7 = 7;
  goto L72;
 L70:
  if (iv <= 1) goto L75;
  mm = 7;
  nn = 1;
  nm = 1;
  i7 = 7;
 L72:
  ll = jj + 4;
  goto L85;
 L75:
  nm = 1;
  i7 = 7;
 L80:
  ll = jj;
 L85:
  for (i = 0; i < 3; ++i) {
    space[ll + 1] = w5 * space[nm - 1];
    space[ll + 2] = w6 * space[i7 - 1];
    if (abs(iv) != 1) {
      space[ll - 1] = w3 * space[mm - 1];
      space[ll    ] = w4 * space[nn - 1];
      ll += 2;
      ++mm;
      ++nn;
    }
    ++nm;
    ++i7;
    ll += 2;
  }
  ll = kk + 2;
  if (abs(iv) == 1) goto L120;
  mm = jj + 4;
  nn = ii + 2;
  i1 = ll;
  for (i = kk-1; i < i1; ++i) {
    space[i] = space[mm] + space[mm + 1] + space[mm + 2] + space[mm + 3] + s[nn] + s[nn + 1];
    mm += 4; 
    nn += 2;
  }
 L115:
  return;
 L120:
  mm = ii + 4;
  nn = jj + 2;
  i1 = ll;
  for (i = kk - 1; i < i1; ++i) {
    space[i] = space[nn] + space[nn + 1] + s[mm] + s[mm + 1] + s[mm + 2] + s[mm + 3];
    mm += 4;
    nn += 2;
  }
  goto L115;
}

//-----------------------------------------------------------------------------------------
// methods ready for testing


void G4BertiniInelasticCollision::abran(G4int k1) {
  // description :
  // parameters  : k1
  // uses        : sign, ce
  // changes     : value1, value2, not
  value1 = sign * G4UniformRand();
  value2 = 0.0; // sum of coulomb energy for a particular region - sum f(i1) mass 
  not = 1;
  for (G4int i = 1; i < k1; ++i) {
    value2 += ce[i];
    if (value2 >= value1) return;
    ++not;
  }
} 

void G4BertiniInelasticCollision::bb() { //::: verify
  // description :
  // parameters  : -
  // uses        : i2, cfepn, wkrpn, rcpmv
  // changes     : energy
  G4int i = i2;          // collision cut-off energy (proton) 
  if (i1 != 0) i += 3;   // collision cut-off energy (neutron) 
  if (knot <= 6) {
    clcfe = cfepn[i];
  } else if (knot > 12) clcfe = cfepn[i];
  energy[0] = wkrpn[i] * rcpmv + massParticle[0]; // total energy particle 1  
  if (in != 0) p1clc();  // p1oe1 = current = moment / total 
  //  else g4cascade->p1cli(); //:::
}

void G4BertiniInelasticCollision::bg6b() {
  // description:
  // uses variables i1, i2, i3, i4
  // changes i1, i3
  zero();
  G4int j = i2;
  for (G4int i = 1; i < i4; ++i) {
    ce[i] = space[j];
    ++j;
  }
  j = i4;
  for (i = j; i < 6; ++i) {
    ce[i] = s[i3];
    ++i3;
  }
}

void G4BertiniInelasticCollision::bg6ca(G4int k, G4int l) {
  // description : collision medium stored in region of initial collision and medium where particle was born 
  clsm = k;
  if (in == 0) {
    curr[9 ] = clsm;
    curr[10] = clsm;
  }
  efrn = space[k + 2] - bindingEnergy;
  efrp = space[k + 8] - bindingEnergy;  // proton fermi energy (MeV) = proton well depth + binding energy
  massParticle[0] = massNucleon;
  massParticle[1] = massNucleon; 
  if (k < l) {
    massParticle[0] = massPionZero;
  } else if (i1 == 0) {
    massParticle[0] = massPionCharged;
  }
}

void G4BertiniInelasticCollision::ccpes() {
  // description :
  // uses        : curr
  // changes     : esps
  G4int i; 
  G4int k;
  G4int l;
  G4int m; 
  G4int n;
  i1 = 0;
  i = static_cast<G4int>(curr[0] + 0.05);
  count[i] += 1.0; // count the number of times each type of particle escapes 
  if (curr[0] != 2.0) k = static_cast<G4int>(curr[9] + 9.05);
  else k = static_cast<G4int>(curr[9] + 3.05);
  if (esps[0] >= 60.0) {
    i1 = 1; // storage already filled
  } else {
    l = static_cast<G4int>(esps[0] * 8.0 + 2.05);
    esps[l] = curr[0];
    esps[l + 1] = curr[1] - space[k];
    m = static_cast<G4int>(13.05 - curr[10]);
    cc[m] += 1.0;
    m = 4;
    l += 2;
    n = l + 2;
    for (G4int i = l-1; i < n; ++i) {
      esps[i    ] = curr[m + 2];
      esps[i + 3] = curr[m - 1];
      ++m;
    }
    esps[0] += 1.0;
  }
}

void G4BertiniInelasticCollision::crjab(G4int    k1, 
                      G4double pp) {
  // decription : 
  // parameters : k1, pp cross-section data
  // uses       : k1, pp, value1
  // changes    : value1, value2
  //  crdet(k1, pp, value1); //:::
  value1 = (pxyz[0] * pxyz[1] + pxyz[4] * pxyz[5] + pxyz[8] * pxyz[9]) / energy[0];  // p1 . p2 / e(1) 
  value2 = value1 / (p2 * p2) * (energy[1] / massNucleon - 1.0) - 1.0 / massNucleon; // s=((p1 . p2) / (e1 * sqr(p2))) * ((e2 / m) - 1.0) - 1.0 / m 
  value2 = massNucleon * crdt[0] * sqrt(sqr(p1oe1) + value1 * 2.0 * value2 + sqr(p2) * sqr(value2)) / (energy[1] * p1oe1 * any); // (m)(c.s) (j^2 + 2 s (p1 . p2) / e1 + (p2) (p2) (s) (s) 
  if (value2 <= 1.0) { // tests sampling technique to ensure fmaxs were selected like value2 <= 1 
    value1 = G4UniformRand(); 
  } else {
    G4Exception("crjab-1"); 
    value1 = G4UniformRand();
  }
}

void G4BertiniInelasticCollision::dfmax() {
  // decription : 
  // parameters :
  // uses       :
  // changes    :
  G4int i;
  G4double wk;
  i = i2;
  if (curr[0] >= 2.0) i += 3;
  wk = wkrpn[i1];
  univ = bovera(wk, massNucleon);
  if (wk < 560.0) {
    fmax[0] = univ * 2.72e-26; // 820  MeV (p, p) s 
    fmax[1] = univ * 3.80e-26; // 230  MeV (p, n) s 
    fmax[2] = univ * 2.26e-26; // 1020 MeV (p, p) s.p. 
    fmax[3] = univ * 1.40e-26; // 750  MeV (p, n) s.p. 
    fmax[4] =        0.0;      //          (p, p) d.p. 
    fmax[5] =        0.0;      //          (p, n) d.p. 
    return;
  }
  if (wk < 800.0) {
    fmax[1] = univ * 3.70e-26; // 250  MeV 
    fmax[4] = univ * 1.90e-27; // 5 and 6 at 1380 MeV 
    fmax[5] = univ * 9.00e-27;
    fmax[0] = univ * 2.72e-26; // 820  MeV
    fmax[2] = univ * 2.26e-26; // 1020 MeV
    fmax[3] = univ * 1.40e-26; // 750  MeV
    return;
  }
  if (wk < 1680.0) {
    fmax[1] = univ * 3.30e-26; // 400  MeV 
    fmax[4] = univ * 1.08e-26; // 5 and 6 at 2600 MeV
    fmax[5] = univ * 1.74e-26;
    fmax[0] = univ * 2.72e-26; // 820  MeV
    fmax[2] = univ * 2.26e-26; // 1020 MeV
    fmax[3] = univ * 1.40e-26; // 750  MeV
    return;
  }
  fmax[0] = univ * 2.63e-26; // 1000 
  fmax[1] = univ * 2.47e-26; // 1000 
  fmax[2] = univ * 2.26e-26; // 1020 
  fmax[3] = univ * 1.35e-26; // 1000 
  fmax[4] = univ * 1.80e-26;
  fmax[5] = univ * 2.36e-26; // 3500 
  return;
}

void G4BertiniInelasticCollision::frmic(G4double &gpart) {
  // decription : find the largest of 3 random numbers
  // parameters : gpart
  // uses       : -
  // changes    : gpart
  G4double r[3];
  for (G4int i = 0; i < 3; ++i) r[i] = G4UniformRand();
  G4double r1 = max(r[0], r[1]);
  G4double r2 = max(r[1], r[2]);
  gpart = max(r1, r2);
}  

void G4BertiniInelasticCollision::idk() {
  // decription : 
  // parameters : -
  // uses       :
  // changes    :
  univ = sqr(pnidk[5]);                         // m(p1) squared decay pion mass squared 
  pnidk[6]  = (sqr(pnidk[0]) + univ - sqnm) / (2 * pnidk[0]); // e(pi) prime decay pion energy prime 
  pnidk[7]  = sqrt(sqr(pnidk[6]) - univ);       // decay pion momentum prime p(d) 
  pol1(pnidk[19], pnidk[20]);                   // cos(theta), sin(theta) 
  azio(pnidk[21], pnidk[22]);                   // cos(phi), sin(phi) 
  pnidk[8]  =  pnidk[21] * pnidk[20] * pnidk[7]; // decay pion x momentum component prime 
  pnidk[9]  =  pnidk[20] * pnidk[22] * pnidk[7]; //            y         
  pnidk[10] =  pnidk[19] * pnidk[7 ];            //            z 
  univ      =  pnidk[8 ] * pnidk[1 ] + pnidk[9] * pnidk[2] + pnidk[10] * pnidk[3]; // p p1 prime dot p 
  pnidk[11] = (pnidk[6 ] * pnidk[4 ] + univ) / pnidk[0]; // decay pion energy e(pi) 
  pnidk[12] =  pnidk[4 ] - pnidk[11];
  univ      = (pnidk[4 ] / pnidk[0 ] - 1.0) * univ / (sqr(pnidk[1]) + sqr(pnidk[2]) + sqr(pnidk[3])); // (e / m - 1.0) * p(p1) prime dot p / p squared 
  unive     =  pnidk[6 ] / pnidk[0 ];              // e pi prime over mass    
  for (G4int i = 0; i < 3; ++i) {
    pnidk[i + 12] = pnidk[i] * (univ + unive) + pnidk[i + 7 ]; // pion momentum components 
    pnidk[i + 15] = pnidk[i]                  - pnidk[i + 12]; // nucleon momentum components 
  }
}

void G4BertiniInelasticCollision::isom() {
  // decription: 
  // parameters:
  // uses:
  // changes:
  G4int m;
  G4double fermn;
  pol1(polc, pols);
  azio(sopc, sops);
  m =  static_cast<G4int>(clsm + 0.05);
  if (strkp + 2.0 > 0.0) fermn = fmpn[m - 1]; // struck proton 
  else fermn = fmpn[m + 2];                   // struck neutron 
  frmic(p2);                                  // select largest of 3 random numbers 
  p2 *= fermn;                                // momentum of particle selected from proper fermi distribution 
  a = p2 * pols;                              // p2 sin(theta) 
  pxyz[1] = a * sopc;                         // p2 sin(theta) cos(phi) 
  pxyz[5] = a * sops;                         // p2 sin(phi) 
  pxyz[9] = p2 * polc;                        // p2 cos(theta) 
  energy[1] = sqrt(sqr(p2) + sqnm);           // sqrt(sqr(momentum struck particle) + sqr(nucleon mass)) 
  rlke = ((energy[0] * energy[1] - pxyz[0] * pxyz[1] - pxyz[4] * pxyz[5] - pxyz[8] * pxyz[9]) / massNucleon - massParticle[0]) / rcpmv; // relative kinetic energy (MeV) - constant = nucleon mass / cm second = MeV / cm. 
}

void G4BertiniInelasticCollision::nn() {
  // decription : initialize vector fmax
  // parameters : -
  // uses       : -
  // changes    : fmax
  fmax[0] = 4.6e-25;
  fmax[1] = 1.4e-24;
  for (G4int i = 1; i < 5; ++i) fmax[i] = 0.0;
}

void G4BertiniInelasticCollision::p1clc() {
  // decription : 
  // parameters : -
  // uses       : energy[0], massParticle[0], p1oe1, curr[6, 7, 8]
  // changes    : p1oe1, pxyz[0, 4, 8]
  p1oe1 = sqrt(sqr(energy[0]) - sqr(massParticle[0]));
  pxyz[0] = p1oe1 * curr[6];
  pxyz[4] = p1oe1 * curr[7];
  pxyz[8] = p1oe1 * curr[8];
  p1oe1 /= energy[0];
}

void G4BertiniInelasticCollision::partin() {
  // decription: keep track of particle numbers in different nucleon regions 
  // parameters: -
  // uses: d[3, 2], ipec[9, 5, 0]
  // changes: ipec[9, 5, 0], isw[0-2]
  //                               // isw[0] = 0 if start in region 3 or 4, 
  //                               // isw[1] = 0 if in rg. 3 not passing through region 1, and 
  //                               // isw[2] = 0 if in rg. 2 not passing through reqion 1 
  for (G4int i = 0; i < 3; ++i) isw[i] = 0;
  if (d[3] != 0.0) {
    ++ipec[9];                     // 1 (number of inc. particles on region 1 only)
  } else if (d[2] != 0.0) {
    ++ipec[5];                     // 2  
  } else ++ipec[0];                // 3
}

void G4BertiniInelasticCollision::pfmax() {
  // decription : set fmax[0-3] 
  // parameters : -
  // uses       : curr[0], wkrpn
  // changes    : fmax[0-3]
  G4int i = i2;
  if (curr[0] >= 2.0) i += 3;
  G4double  wk = wkrpn[i - 1];
  if (wk < 160.0) {
    fmax[0] = 1.76e-25;
    fmax[1] = 5.20e-25;
    fmax[2] = 0.0;
    fmax[3] = 0.0;
    return;
  }
  univ = bovera(wk, massNucleon);
  if (wk < 200.0) {
    fmax[2] = univ * 1.95e-27; // 3 and 4 at 465 MeV 
    fmax[3] = univ * 1.70e-27;
    fmax[0] =        1.03e-25;
    fmax[1] =        3.13e-25;
  }
  if (wk < 300.0) {
    fmax[0] =        9.00e-26; // 1 and 2 at 35 MeV
    fmax[1] =        2.60e-25;
    fmax[2] = univ * 1.14e-26; // 3 and 4 at 630 MeV 
    fmax[3] = univ * 1.12e-26;
  }
  if (wk < 400.0){
    fmax[0] = univ * 2.80e-26; // 1 and 2 at 100 MeV 
    fmax[1] =        7.30e-26;
    fmax[2] = univ * 2.00e-26; // 3 and 4 at 780 MeV 
    fmax[3] = univ * 1.40e-26;
  }
  if (wk >= 400.0) {
    fmax[0] = univ * 2.72e-26; // 1 and 2 at 155 MeV 
    fmax[1] = univ * 4.80e-26; // fmax(1) = (p, p) s - (2) = (p, n) s - (3) = (p, p) s.p. - (4) = (p, n) s.p. //::: index 0/1? 
    fmax[2] = univ * 2.26e-26; // 3 and 4 at 1020 MeV 
    fmax[3] = univ * 1.40e-26;
  }
}

void G4BertiniInelasticCollision::spac32(G4int i) {  // ::: verify
  // decription : 
  // parameters :
  // uses       :
  // changes    :
  ex = 0.0;
  i4 = 5;
  sign = space[i] * 0.999999;
  if (i < 31) {
    i2 = 18;
    i3 = 3; 
  } else {
    if (i == 31) {
      i2 = 22;
      i3 = 5;
    } else {
      if (i < 41) {
	i2 = 26;
	i3 = 7;
      } else {
	i4 = 3;
	if (i < 42) {
	  i2 = 35;
	  i3 = 13;
	} else if (i == 42) {
	  i2 = 37;
	  i3 = 17;
	} else {
	  i2 = 39;
	  i3 = 21;
	}
      }
    }
  }
  bg6b();
  signex();
} 


void G4BertiniInelasticCollision::spisom() {
  // decription : set space, e and pxyz vectors 
  // parameters : -
  // uses       : isom(), pxyz, space
  // changes    : space[175-177], pxyz[1, 5, 9], energy[1]
  while(1) {
    isom();
    if (i1 <= 0) {
      space[175] = pxyz[1];
      space[176] = pxyz[5];
      space[177] = pxyz[9];
      if (i1 >= 0) return;
      else i1 = 1;
    } else {
      pxyz[1] += space[175];
      pxyz[5] += space[176];
      pxyz[9] += space[177];
      energy[1] = (sqr(pxyz[9]) + sqr(pxyz[5]) + sqr(pxyz[1])) / 1.9032e14;
      return; 
    }
  }
}

void G4BertiniInelasticCollision::stph(G4double *w) { // ::: verify
  // decription : 
  // parameters :
  // uses       :
  // changes    : 
  if (pgvc[0] >= 440.0) { // pgvc(0) = number of times velocity greater than criterion entered 
    i1 = 1;
  } else {
  G4int i = static_cast<G4int>(pgvc[0] + 1.005);
  G4int j = i + 10;
  G4int l = 2;
  for (G4int k = i; k < j; ++k) { 
    pgvc[k] = w[l]; 
    ++l;
  }
  pgvc[0] += 11.0;
  }
  return;
}

void G4BertiniInelasticCollision::stpl(G4double *w) { 
  // decription:  set plvc
  // parameters:
  // uses:
  // changes:
  G4int i;
  G4int j;
  G4int l;
  for (G4int k = 1; k < 950; k += 12) {
    if (plvc[k] == 0.0) {
      i = k;
      j = i + 11;
      l = 1;
      for (G4int k = i; k < j; ++k) {
	plvc[k] = w[l]; // plvc[0] is the number of times entered for storage of velocity less than criterion 
	++l;
      }
      plvc[0] += 1.0;
      return;
    }
  }
  i1 = 1;
}

void G4BertiniInelasticCollision::zero() {
  // decription : initialize vector ce to 0.0
  // parameters : -
  // uses       : -
  // changes    : ce
  for (G4int i = 0; i < 20; ++i) ce[i] = 0.0;
}

void G4BertiniInelasticCollision::signex() {
  // decription : update ex = distance in sampling routine exponential random divided by sigma ci region i  
  // parameters : -
  // uses       : exprn(), sign
  // changes    : ex
  ex += exprn() / sign;
}




