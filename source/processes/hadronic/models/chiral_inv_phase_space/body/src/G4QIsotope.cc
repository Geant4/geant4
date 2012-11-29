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
//
// $Id$
//
//      ---------------- G4QIsotope class ----------------
//             by Mikhail Kossov, December 2003.
//  class G4QIsotope gives Isotopes for Elements (CHIPS branch of G4)
// ------------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
//
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
// ------------------------------------------------------------------------------
// Short description: containes the natural abundance DB for isotops and permits
// new artificial materials with unnatural abundance (e.g. enreached Deuterium).
// Uses isotope cross-sections for calculation of the mean cross-sections for the
// Element (fixed Z).
// ------------------------------------------------------------------------------


//#define debug
//#define cdebug
//#define pdebug
//#define ppdebug
//#define sdebug

#include "G4QIsotope.hh"
#include <cmath>
using namespace std;

vector<vector<pair<G4int,G4double>*>*>G4QIsotope::natElements;//NaturalElems
vector<vector<pair<G4int,G4double>*>*>G4QIsotope::natSumAbund;//NaturalSumAb
vector<vector<pair<G4int,G4double>*>*>G4QIsotope::natIsoCrosS;//CSOfNatElems
vector<pair<G4int,vector<pair<G4int,G4double>*>*>*>G4QIsotope::newElems;
vector<pair<G4int,vector<pair<G4int,G4double>*>*>*>G4QIsotope::newSumAb;
vector<pair<G4int,vector<pair<G4int,G4double>*>*>*>G4QIsotope::newIsoCS;

G4QIsotope::G4QIsotope() 
{
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor is called"<<G4endl;
#endif
  vector<vector<pair<G4int,G4double> >*> natEl;
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor natEl is booked"<<G4endl;
#endif
  vector<pair<G4int,G4double> >*a0=new vector<pair<G4int,G4double> >;
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor a0 is booked"<<G4endl;
#endif
  a0->push_back(make_pair(1,1.));
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor a0 is filled by a pair"<<G4endl;
#endif
  natEl.push_back(a0);
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor a0 is filled in natEl"<<G4endl;
#endif
  // If an error is found in this initialization, please, correct the simple tree below
  vector<pair<G4int,G4double> >*a1=new vector<pair<G4int,G4double> >; // 1-H
  a1->push_back(make_pair(0,.99985));
  a1->push_back(make_pair(1,1.));
  natEl.push_back(a1);
  vector<pair<G4int,G4double> >*a2=new vector<pair<G4int,G4double> >; // 2-He
  a2->push_back(make_pair(2,.999999863));
  a2->push_back(make_pair(1,1.));
  natEl.push_back(a2);
  vector<pair<G4int,G4double> >*a3=new vector<pair<G4int,G4double> >; // 3-Li
  a3->push_back(make_pair(4,.925));
  a3->push_back(make_pair(3,1.));
  natEl.push_back(a3);
  vector<pair<G4int,G4double> >*a4=new vector<pair<G4int,G4double> >; // 4-Be
  a4->push_back(make_pair(5,1.));
  natEl.push_back(a4);
  vector<pair<G4int,G4double> >*a5=new vector<pair<G4int,G4double> >; // 5-B
  a5->push_back(make_pair(6,.801));
  a5->push_back(make_pair(5,1.));
  natEl.push_back(a5);
  vector<pair<G4int,G4double> >*a6=new vector<pair<G4int,G4double> >; // 6-C
  a6->push_back(make_pair(6,.989));
  a6->push_back(make_pair(7,1.));
  natEl.push_back(a6);
  vector<pair<G4int,G4double> >*a7=new vector<pair<G4int,G4double> >; // 7-N
  a7->push_back(make_pair(7,.9963));
  a7->push_back(make_pair(8,1.));
  natEl.push_back(a7);
  vector<pair<G4int,G4double> >*a8=new vector<pair<G4int,G4double> >; // 8-O
  a8->push_back(make_pair(8,.9976));
  a8->push_back(make_pair(10,.9996));
  a8->push_back(make_pair(9,1.));
  natEl.push_back(a8);
  vector<pair<G4int,G4double> >*a9=new vector<pair<G4int,G4double> >; // 9-F
  a9->push_back(make_pair(10,1.));
  natEl.push_back(a9);
  vector<pair<G4int,G4double> >*b0=new vector<pair<G4int,G4double> >; // 10-Ne
  b0->push_back(make_pair(10,.9948));
  b0->push_back(make_pair(11,.9975));
  b0->push_back(make_pair(12,1.));
  natEl.push_back(b0);
  vector<pair<G4int,G4double> >*b1=new vector<pair<G4int,G4double> >; // 11-Na
  b1->push_back(make_pair(12,1.));
  natEl.push_back(b1);
  vector<pair<G4int,G4double> >*b2=new vector<pair<G4int,G4double> >; // 12-Mg
  b2->push_back(make_pair(12,.7899));
  b2->push_back(make_pair(13,.8899));
  b2->push_back(make_pair(14,1.));
  natEl.push_back(b2);
  vector<pair<G4int,G4double> >*b3=new vector<pair<G4int,G4double> >; // 13-Al
  b3->push_back(make_pair(14,1.));
  natEl.push_back(b3);
  vector<pair<G4int,G4double> >*b4=new vector<pair<G4int,G4double> >; // 14-Si
  b4->push_back(make_pair(14,.9223));
  b4->push_back(make_pair(15,.969));
  b4->push_back(make_pair(16,1.));
  natEl.push_back(b4);
  vector<pair<G4int,G4double> >*b5=new vector<pair<G4int,G4double> >; // 15-P
  b5->push_back(make_pair(16,1.));
  natEl.push_back(b5);
  vector<pair<G4int,G4double> >*b6=new vector<pair<G4int,G4double> >; // 16-S
  b6->push_back(make_pair(16,.9502));
  b6->push_back(make_pair(18,.9923));
  b6->push_back(make_pair(17,.9998));
  b6->push_back(make_pair(20,1.));
  natEl.push_back(b6);
  vector<pair<G4int,G4double> >*b7=new vector<pair<G4int,G4double> >; // 17-Cl
  b7->push_back(make_pair(18,.7577));
  b7->push_back(make_pair(20,1.));
  natEl.push_back(b7);
  vector<pair<G4int,G4double> >*b8=new vector<pair<G4int,G4double> >; // 18-Ar
  b8->push_back(make_pair(22,.996));
  b8->push_back(make_pair(18,.99937));
  b8->push_back(make_pair(20,1.));
  natEl.push_back(b8);
  vector<pair<G4int,G4double> >*b9=new vector<pair<G4int,G4double> >; // 19-K
  b9->push_back(make_pair(20,.932581));
  b9->push_back(make_pair(22,.999883));
  b9->push_back(make_pair(21,1.));
  natEl.push_back(b9);
  vector<pair<G4int,G4double> >*c0=new vector<pair<G4int,G4double> >; // 20-Ca
  c0->push_back(make_pair(20,.96941));
  c0->push_back(make_pair(24,.99027));
  c0->push_back(make_pair(22,.99674));
  c0->push_back(make_pair(28,.99861));
  c0->push_back(make_pair(23,.99996));
  c0->push_back(make_pair(26,1.));
  natEl.push_back(c0);
  vector<pair<G4int,G4double> >*c1=new vector<pair<G4int,G4double> >; // 21-Sc
  c1->push_back(make_pair(24,1.));
  natEl.push_back(c1);
  vector<pair<G4int,G4double> >*c2=new vector<pair<G4int,G4double> >; // 22-Ti
  c2->push_back(make_pair(26,.738));
  c2->push_back(make_pair(24,.818));
  c2->push_back(make_pair(25,.891));
  c2->push_back(make_pair(27,.946));
  c2->push_back(make_pair(28,1.));
  natEl.push_back(c2);
  vector<pair<G4int,G4double> >*c3=new vector<pair<G4int,G4double> >; // 23-V
  c3->push_back(make_pair(28,.9975));
  c3->push_back(make_pair(27,1.));
  natEl.push_back(c3);
  vector<pair<G4int,G4double> >*c4=new vector<pair<G4int,G4double> >; // 24-Cr
  c4->push_back(make_pair(28,.8379));
  c4->push_back(make_pair(29,.9329));
  c4->push_back(make_pair(26,.97635));
  c4->push_back(make_pair(30,1.));
  natEl.push_back(c4);
  vector<pair<G4int,G4double> >*c5=new vector<pair<G4int,G4double> >; // 25-Mn
  c5->push_back(make_pair(30,1.));
  natEl.push_back(c5);
  vector<pair<G4int,G4double> >*c6=new vector<pair<G4int,G4double> >; // 26-Fe
  c6->push_back(make_pair(30,.9172));
  c6->push_back(make_pair(28,.9762));
  c6->push_back(make_pair(31,.9972));
  c6->push_back(make_pair(32,1.));
  natEl.push_back(c6);
  vector<pair<G4int,G4double> >*c7=new vector<pair<G4int,G4double> >; // 27-Co
  c7->push_back(make_pair(32,1.));
  natEl.push_back(c7);
  vector<pair<G4int,G4double> >*c8=new vector<pair<G4int,G4double> >; // 28-Ni
  c8->push_back(make_pair(30,.68077));
  c8->push_back(make_pair(32,.943));
  c8->push_back(make_pair(34,.97934));
  c8->push_back(make_pair(33,.99074));
  c8->push_back(make_pair(36,1.));
  natEl.push_back(c8);
  vector<pair<G4int,G4double> >*c9=new vector<pair<G4int,G4double> >; // 29-Cu
  c9->push_back(make_pair(34,.6917));
  c9->push_back(make_pair(36,1.));
  natEl.push_back(c9);
  vector<pair<G4int,G4double> >*d0=new vector<pair<G4int,G4double> >; // 30-Zn
  d0->push_back(make_pair(34,.486));
  d0->push_back(make_pair(36,.765));
  d0->push_back(make_pair(38,.953));
  d0->push_back(make_pair(37,.994));
  d0->push_back(make_pair(40,1.));
  natEl.push_back(d0);
  vector<pair<G4int,G4double> >*d1=new vector<pair<G4int,G4double> >; // 31-Ga
  d1->push_back(make_pair(38,.60108));
  d1->push_back(make_pair(40,1.));
  natEl.push_back(d1);
  vector<pair<G4int,G4double> >*d2=new vector<pair<G4int,G4double> >; // 32-Ge
  d2->push_back(make_pair(42,.3594));
  d2->push_back(make_pair(40,.6360));
  d2->push_back(make_pair(38,.8484));
  d2->push_back(make_pair(41,.9256));
  d2->push_back(make_pair(44,1.));
  natEl.push_back(d2);
  vector<pair<G4int,G4double> >*d3=new vector<pair<G4int,G4double> >; // 33-As
  d3->push_back(make_pair(42,1.));
  natEl.push_back(d3);
  vector<pair<G4int,G4double> >*d4=new vector<pair<G4int,G4double> >; // 34-Se
  d4->push_back(make_pair(46,.4961));
  d4->push_back(make_pair(44,.7378));
  d4->push_back(make_pair(42,.8274));
  d4->push_back(make_pair(48,.9148));
  d4->push_back(make_pair(43,.9911));
  d4->push_back(make_pair(40,1.));
  natEl.push_back(d4);
  vector<pair<G4int,G4double> >*d5=new vector<pair<G4int,G4double> >; // 35-Br
  d5->push_back(make_pair(44,.5069));
  d5->push_back(make_pair(46,1.));
  natEl.push_back(d5);
  vector<pair<G4int,G4double> >*d6=new vector<pair<G4int,G4double> >; // 36-Kr
  d6->push_back(make_pair(48,.57));
  d6->push_back(make_pair(50,.743));
  d6->push_back(make_pair(46,.859));
  d6->push_back(make_pair(47,.974));
  d6->push_back(make_pair(44,.9965));
  d6->push_back(make_pair(42,1.));
  natEl.push_back(d6);
  vector<pair<G4int,G4double> >*d7=new vector<pair<G4int,G4double> >; // 37-Rb
  d7->push_back(make_pair(48,.7217));
  d7->push_back(make_pair(50,1.));
  natEl.push_back(d7);
  vector<pair<G4int,G4double> >*d8=new vector<pair<G4int,G4double> >; // 38-sr
  d8->push_back(make_pair(50,.8258));
  d8->push_back(make_pair(48,.9244));
  d8->push_back(make_pair(49,.9944));
  d8->push_back(make_pair(46,1.));
  natEl.push_back(d8);
  vector<pair<G4int,G4double> >*d9=new vector<pair<G4int,G4double> >; // 39-Y
  d9->push_back(make_pair(50,1.));
  natEl.push_back(d9);
  vector<pair<G4int,G4double> >*e0=new vector<pair<G4int,G4double> >; // 40-Zr
  e0->push_back(make_pair(50,.5145));
  e0->push_back(make_pair(54,.6883));
  e0->push_back(make_pair(52,.8598));
  e0->push_back(make_pair(51,.972));
  e0->push_back(make_pair(56,1.));
  natEl.push_back(e0);
  vector<pair<G4int,G4double> >*e1=new vector<pair<G4int,G4double> >; // 41-Nb
  e1->push_back(make_pair(52,1.));
  natEl.push_back(e1);
  vector<pair<G4int,G4double> >*e2=new vector<pair<G4int,G4double> >; // 42-Mo
  e2->push_back(make_pair(56,.2413));
  e2->push_back(make_pair(54,.4081));
  e2->push_back(make_pair(53,.5673));
  e2->push_back(make_pair(50,.7157));
  e2->push_back(make_pair(58,.8120));
  e2->push_back(make_pair(55,.9075));
  e2->push_back(make_pair(52,1.));
  natEl.push_back(e2);
  vector<pair<G4int,G4double> >*e3=new vector<pair<G4int,G4double> >; // 43-Tc
  e3->push_back(make_pair(55,1.));
  natEl.push_back(e3);
  vector<pair<G4int,G4double> >*e4=new vector<pair<G4int,G4double> >; // 44-Ru
  e4->push_back(make_pair(58,.316));
  e4->push_back(make_pair(60,.502));
  e4->push_back(make_pair(57,.673));
  e4->push_back(make_pair(55,.8));
  e4->push_back(make_pair(56,.926));
  e4->push_back(make_pair(52,.9814));
  e4->push_back(make_pair(54,1.));
  natEl.push_back(e4);
  vector<pair<G4int,G4double> >*e5=new vector<pair<G4int,G4double> >; // 45-Rh
  e5->push_back(make_pair(58,1.));
  natEl.push_back(e5);
  vector<pair<G4int,G4double> >*e6=new vector<pair<G4int,G4double> >; // 46-Pd
  e6->push_back(make_pair(60,.2733));
  e6->push_back(make_pair(62,.5379));
  e6->push_back(make_pair(59,.7612));
  e6->push_back(make_pair(55,.8784));
  e6->push_back(make_pair(58,.9898));
  e6->push_back(make_pair(56,1.));
  natEl.push_back(e6);
  vector<pair<G4int,G4double> >*e7=new vector<pair<G4int,G4double> >; // 47-Ag
  e7->push_back(make_pair(60,.51839));
  e7->push_back(make_pair(62,1.));
  natEl.push_back(e7);
  vector<pair<G4int,G4double> >*e8=new vector<pair<G4int,G4double> >; // 48-Cd
  e8->push_back(make_pair(66,.2873));
  e8->push_back(make_pair(64,.5286));
  e8->push_back(make_pair(59,.6566));
  e8->push_back(make_pair(62,.7815));
  e8->push_back(make_pair(65,.9037));
  e8->push_back(make_pair(68,.9786));
  e8->push_back(make_pair(58,.9911));
  e8->push_back(make_pair(60,1.));
  natEl.push_back(e8);
  vector<pair<G4int,G4double> >*e9=new vector<pair<G4int,G4double> >; // 49-In
  e9->push_back(make_pair(66,.9577));
  e9->push_back(make_pair(64,1.));
  natEl.push_back(e9);
  vector<pair<G4int,G4double> >*f0=new vector<pair<G4int,G4double> >; // 50-Sn
  f0->push_back(make_pair(70,.3259));
  f0->push_back(make_pair(68,.5681));
  f0->push_back(make_pair(66,.7134));
  f0->push_back(make_pair(69,.7992));
  f0->push_back(make_pair(67,.8760));
  f0->push_back(make_pair(74,.9339));
  f0->push_back(make_pair(72,.9802));
  f0->push_back(make_pair(62,.9899));
  f0->push_back(make_pair(64,1.));
  //f0->push_back(make_pair(64,.9964));
  //f0->push_back(make_pair(65,1.)); // Nine isotopes is the maximum, so Sn115 is out
  natEl.push_back(f0);
  vector<pair<G4int,G4double> >*f1=new vector<pair<G4int,G4double> >; // 51-Sb
  f1->push_back(make_pair(70,.5736));
  f1->push_back(make_pair(72,1.));
  natEl.push_back(f1);
  vector<pair<G4int,G4double> >*f2=new vector<pair<G4int,G4double> >; // 52-Te
  f2->push_back(make_pair(78,.3387));
  f2->push_back(make_pair(76,.6557));
  f2->push_back(make_pair(74,.8450));
  f2->push_back(make_pair(73,.9162));
  f2->push_back(make_pair(72,.9641));
  f2->push_back(make_pair(70,.9900));
  f2->push_back(make_pair(71,.99905));
  f2->push_back(make_pair(68,1.));
  natEl.push_back(f2);
  vector<pair<G4int,G4double> >*f3=new vector<pair<G4int,G4double> >; // 53-I
  f3->push_back(make_pair(74,1.));
  natEl.push_back(f3);
  vector<pair<G4int,G4double> >*f4=new vector<pair<G4int,G4double> >; // 54-Xe
  f4->push_back(make_pair(78,.269));
  f4->push_back(make_pair(75,.533));
  f4->push_back(make_pair(77,.745));
  f4->push_back(make_pair(80,.849));
  f4->push_back(make_pair(82,.938));
  f4->push_back(make_pair(76,.979));
  f4->push_back(make_pair(74,.9981));
  f4->push_back(make_pair(70,.9991));
  f4->push_back(make_pair(72,1.));
  natEl.push_back(f4);
  vector<pair<G4int,G4double> >*f5=new vector<pair<G4int,G4double> >; // 55-Cs
  f5->push_back(make_pair(78,1.));
  natEl.push_back(f5);
  vector<pair<G4int,G4double> >*f6=new vector<pair<G4int,G4double> >; // 56-Ba
  f6->push_back(make_pair(82,.717));
  f6->push_back(make_pair(81,.8293));
  f6->push_back(make_pair(80,.9078));
  f6->push_back(make_pair(79,.97373));
  f6->push_back(make_pair(78,.99793));
  f6->push_back(make_pair(74,.99899));
  f6->push_back(make_pair(76,1.));
  natEl.push_back(f6);
  vector<pair<G4int,G4double> >*f7=new vector<pair<G4int,G4double> >; // 57-La
  f7->push_back(make_pair(82,.999098));
  f7->push_back(make_pair(81,1.));
  natEl.push_back(f7);
  vector<pair<G4int,G4double> >*f8=new vector<pair<G4int,G4double> >; // 58-Ce
  f8->push_back(make_pair(82,.8843));
  f8->push_back(make_pair(84,.9956));
  f8->push_back(make_pair(80,.9981));
  f8->push_back(make_pair(78,1.));
  natEl.push_back(f8);
  vector<pair<G4int,G4double> >*f9=new vector<pair<G4int,G4double> >; // 59-Pr
  f9->push_back(make_pair(82,1.));
  natEl.push_back(f9);
  vector<pair<G4int,G4double> >*g0=new vector<pair<G4int,G4double> >; // 60-Nd
  g0->push_back(make_pair(82,.2713));
  g0->push_back(make_pair(84,.5093));
  g0->push_back(make_pair(86,.6812));
  g0->push_back(make_pair(83,.8030));
  g0->push_back(make_pair(85,.8860));
  g0->push_back(make_pair(88,.9436));
  g0->push_back(make_pair(90,1.));
  natEl.push_back(g0);
  vector<pair<G4int,G4double> >*g1=new vector<pair<G4int,G4double> >; // 61-Pm
  g1->push_back(make_pair(85,1.));
  natEl.push_back(g1);
  vector<pair<G4int,G4double> >*g2=new vector<pair<G4int,G4double> >; // 62-Sm
  g2->push_back(make_pair(90,.267));
  g2->push_back(make_pair(92,.494));
  g2->push_back(make_pair(85,.644));
  g2->push_back(make_pair(87,.782));
  g2->push_back(make_pair(86,.895));
  g2->push_back(make_pair(88,.969));
  g2->push_back(make_pair(82,1.));
  natEl.push_back(g2);
  vector<pair<G4int,G4double> >*g3=new vector<pair<G4int,G4double> >; // 63-Eu
  g3->push_back(make_pair(90,.522));
  g3->push_back(make_pair(89,1.));
  natEl.push_back(g3);
  vector<pair<G4int,G4double> >*g4=new vector<pair<G4int,G4double> >; // 64-Gd
  g4->push_back(make_pair(94,.2484));
  g4->push_back(make_pair(96,.4670));
  g4->push_back(make_pair(92,.6717));
  g4->push_back(make_pair(93,.8282));
  g4->push_back(make_pair(91,.9762));
  g4->push_back(make_pair(90,.9980));
  g4->push_back(make_pair(88,1.));
  natEl.push_back(g4);
  vector<pair<G4int,G4double> >*g5=new vector<pair<G4int,G4double> >; // 65-Tb
  g5->push_back(make_pair(94,1.));
  natEl.push_back(g5);
  vector<pair<G4int,G4double> >*g6=new vector<pair<G4int,G4double> >; // 66-Dy
  g6->push_back(make_pair(98,.282));
  g6->push_back(make_pair(96,.537));
  g6->push_back(make_pair(97,.786));
  g6->push_back(make_pair(95,.975));
  g6->push_back(make_pair(94,.9984));
  g6->push_back(make_pair(92,.9994));
  g6->push_back(make_pair(90,1.));
  natEl.push_back(g6);
  vector<pair<G4int,G4double> >*g7=new vector<pair<G4int,G4double> >; // 67-Ho
  g7->push_back(make_pair(98,1.));
  natEl.push_back(g7);
  vector<pair<G4int,G4double> >*g8=new vector<pair<G4int,G4double> >; // 68-Er
  g8->push_back(make_pair( 98,.3360));
  g8->push_back(make_pair(100,.6040));
  g8->push_back(make_pair( 99,.8335));
  g8->push_back(make_pair(102,.9825));
  g8->push_back(make_pair( 96,.9986));
  g8->push_back(make_pair( 94,1.));
  natEl.push_back(g8);
  vector<pair<G4int,G4double> >*g9=new vector<pair<G4int,G4double> >; // 69-Tm
  g9->push_back(make_pair(100,1.));
  natEl.push_back(g9);
  vector<pair<G4int,G4double> >*h0=new vector<pair<G4int,G4double> >; // 70-Yb
  h0->push_back(make_pair(104,.3180));
  h0->push_back(make_pair(102,.5370));
  h0->push_back(make_pair(103,.6982));
  h0->push_back(make_pair(101,.8412));
  h0->push_back(make_pair(106,.9682));
  h0->push_back(make_pair(100,.9987));
  h0->push_back(make_pair( 98,1.));
  natEl.push_back(h0);
  vector<pair<G4int,G4double> >*h1=new vector<pair<G4int,G4double> >; // 71-Lu
  h1->push_back(make_pair(104,.9741));
  h1->push_back(make_pair(105,1.));
  natEl.push_back(h1);
  vector<pair<G4int,G4double> >*h2=new vector<pair<G4int,G4double> >; // 72-Hf
  h2->push_back(make_pair(108,.35100));
  h2->push_back(make_pair(106,.62397));
  h2->push_back(make_pair(105,.81003));
  h2->push_back(make_pair(107,.94632));
  h2->push_back(make_pair(104,.99838));
  h2->push_back(make_pair(102,1.));
  natEl.push_back(h2);
  vector<pair<G4int,G4double> >*h3=new vector<pair<G4int,G4double> >; // 73-Ta
  h3->push_back(make_pair(108,.99988));
  h3->push_back(make_pair(107,1.));
  natEl.push_back(h3);
  vector<pair<G4int,G4double> >*h4=new vector<pair<G4int,G4double> >; // 74-W
  h4->push_back(make_pair(110,.307));
  h4->push_back(make_pair(112,.593));
  h4->push_back(make_pair(108,.856));
  h4->push_back(make_pair(109,.9988));
  h4->push_back(make_pair(106,1.));
  natEl.push_back(h4);
  vector<pair<G4int,G4double> >*h5=new vector<pair<G4int,G4double> >; // 75-Re
  h5->push_back(make_pair(112,.626));
  h5->push_back(make_pair(110,1.));
  natEl.push_back(h5);
  vector<pair<G4int,G4double> >*h6=new vector<pair<G4int,G4double> >; // 78-Os
  h6->push_back(make_pair(116,.410));
  h6->push_back(make_pair(114,.674));
  h6->push_back(make_pair(113,.835));
  h6->push_back(make_pair(112,.968));
  h6->push_back(make_pair(111,.984));
  h6->push_back(make_pair(110,.9998));
  h6->push_back(make_pair(108,1.));
  natEl.push_back(h6);
  vector<pair<G4int,G4double> >*h7=new vector<pair<G4int,G4double> >; // 77-Ir
  h7->push_back(make_pair(116,.627));
  h7->push_back(make_pair(114,1.));
  natEl.push_back(h7);
  vector<pair<G4int,G4double> >*h8=new vector<pair<G4int,G4double> >; // 78-Pt
  h8->push_back(make_pair(117,.338));
  h8->push_back(make_pair(116,.667));
  h8->push_back(make_pair(118,.920));
  h8->push_back(make_pair(120,.992));
  h8->push_back(make_pair(114,.9999));
  h8->push_back(make_pair(112,1.));
  natEl.push_back(h8);
  vector<pair<G4int,G4double> >*h9=new vector<pair<G4int,G4double> >; // 79-Au
  h9->push_back(make_pair(118,1.));
  natEl.push_back(h9);
  vector<pair<G4int,G4double> >*i0=new vector<pair<G4int,G4double> >; // 80-Hg
  i0->push_back(make_pair(122,.2986));
  i0->push_back(make_pair(120,.5296));
  i0->push_back(make_pair(119,.6983));
  i0->push_back(make_pair(121,.8301));
  i0->push_back(make_pair(118,.9298));
  i0->push_back(make_pair(124,.9985));
  i0->push_back(make_pair(116,1.));
  natEl.push_back(i0);
  vector<pair<G4int,G4double> >*i1=new vector<pair<G4int,G4double> >; // 81-Tl
  i1->push_back(make_pair(124,.70476));
  i1->push_back(make_pair(122,1.));
  natEl.push_back(i1);
  vector<pair<G4int,G4double> >*i2=new vector<pair<G4int,G4double> >; // 82-Pb
  i2->push_back(make_pair(126,.524));
  i2->push_back(make_pair(124,.765));
  i2->push_back(make_pair(125,.986));
  i2->push_back(make_pair(122,1.));
  natEl.push_back(i2);
  vector<pair<G4int,G4double> >*i3=new vector<pair<G4int,G4double> >; // 83-Bi
  i3->push_back(make_pair(126,1.));
  natEl.push_back(i3);
  vector<pair<G4int,G4double> >*i4=new vector<pair<G4int,G4double> >; // 84-Po
  i4->push_back(make_pair(125,1.));
  natEl.push_back(i4);
  vector<pair<G4int,G4double> >*i5=new vector<pair<G4int,G4double> >; // 85-At
  i5->push_back(make_pair(136,1.));
  natEl.push_back(i5);
  vector<pair<G4int,G4double> >*i6=new vector<pair<G4int,G4double> >; // 86-Ru
  i6->push_back(make_pair(136,1.));
  natEl.push_back(i6);
  vector<pair<G4int,G4double> >*i7=new vector<pair<G4int,G4double> >; // 87-Fr
  i7->push_back(make_pair(138,1.));
  natEl.push_back(i7);
  vector<pair<G4int,G4double> >*i8=new vector<pair<G4int,G4double> >; // 88-Ra
  i8->push_back(make_pair(138,1.));
  natEl.push_back(i8);
  vector<pair<G4int,G4double> >*i9=new vector<pair<G4int,G4double> >; // 89-Ac
  i9->push_back(make_pair(142,1.));
  natEl.push_back(i9);
  vector<pair<G4int,G4double> >*j0=new vector<pair<G4int,G4double> >; // 90-Th
  j0->push_back(make_pair(142,1.));
  natEl.push_back(j0);
  vector<pair<G4int,G4double> >*j1=new vector<pair<G4int,G4double> >; // 91-Pa
  j1->push_back(make_pair(140,1.));
  natEl.push_back(j1);
  vector<pair<G4int,G4double> >*j2=new vector<pair<G4int,G4double> >; // 92-U
  j2->push_back(make_pair(146,.992745));
  j2->push_back(make_pair(143,.999945));
  j2->push_back(make_pair(142,1.));
  natEl.push_back(j2);
  vector<pair<G4int,G4double> >*j3=new vector<pair<G4int,G4double> >; // 93-Np
  j3->push_back(make_pair(144,1.));
  natEl.push_back(j3);
  vector<pair<G4int,G4double> >*j4=new vector<pair<G4int,G4double> >; // 94-Pu
  j4->push_back(make_pair(150,1.));
  natEl.push_back(j4);
  vector<pair<G4int,G4double> >*j5=new vector<pair<G4int,G4double> >; // 95-Am
  j5->push_back(make_pair(148,1.));
  natEl.push_back(j5);
  vector<pair<G4int,G4double> >*j6=new vector<pair<G4int,G4double> >; // 96-Cm
  j6->push_back(make_pair(151,1.));
  natEl.push_back(j6);
  vector<pair<G4int,G4double> >*j7=new vector<pair<G4int,G4double> >; // 97-Bk
  j7->push_back(make_pair(150,1.));
  natEl.push_back(j7);
  vector<pair<G4int,G4double> >*j8=new vector<pair<G4int,G4double> >; // 98-Cf
  j8->push_back(make_pair(153,1.));
  natEl.push_back(j8);
  vector<pair<G4int,G4double> >*j9=new vector<pair<G4int,G4double> >; // 99-Es
  j9->push_back(make_pair(157,1.));
  natEl.push_back(j9);
  vector<pair<G4int,G4double> >*k0=new vector<pair<G4int,G4double> >; // 100-Fm
  k0->push_back(make_pair(157,1.));
  natEl.push_back(k0);
  vector<pair<G4int,G4double> >*k1=new vector<pair<G4int,G4double> >; // 101-Md
  k1->push_back(make_pair(157,1.));
  natEl.push_back(k1);
  vector<pair<G4int,G4double> >*k2=new vector<pair<G4int,G4double> >; // 102-No
  k2->push_back(make_pair(157,1.));
  natEl.push_back(k2);
  vector<pair<G4int,G4double> >*k3=new vector<pair<G4int,G4double> >; // 103-Lr
  k3->push_back(make_pair(157,1.));
  natEl.push_back(k3);
  vector<pair<G4int,G4double> >*k4=new vector<pair<G4int,G4double> >; // 104-Rf
  k4->push_back(make_pair(157,1.));
  natEl.push_back(k4);
  vector<pair<G4int,G4double> >*k5=new vector<pair<G4int,G4double> >; // 105-Db
  k5->push_back(make_pair(157,1.));
  natEl.push_back(k5);
  vector<pair<G4int,G4double> >*k6=new vector<pair<G4int,G4double> >; // 106-Sg
  k6->push_back(make_pair(157,1.));
  natEl.push_back(k6);
  vector<pair<G4int,G4double> >*k7=new vector<pair<G4int,G4double> >; // 107-Bh
  k7->push_back(make_pair(155,1.));
  natEl.push_back(k7);
  vector<pair<G4int,G4double> >*k8=new vector<pair<G4int,G4double> >; // 108-Hs
  k8->push_back(make_pair(157,1.));
  natEl.push_back(k8);
  vector<pair<G4int,G4double> >*k9=new vector<pair<G4int,G4double> >; // 109-Mt
  k9->push_back(make_pair(157,1.));
  natEl.push_back(k9);
  // Now fill natElements and natIsoCrossS
  G4int nona=natEl.size();
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor natEl filling is finished nE="<<nona<<G4endl;
#endif
  for(G4int i=0; i<nona; i++)
  {
    vector<pair<G4int,G4double> >* is=natEl[i]; // Pointer to theElement
    G4int n=is->size();
#ifdef cdebug
    G4cout<<"G4QIsotope::Constructor: Element # "<<i<<", nOfIsotopes="<<n<<G4endl;
#endif
    vector<pair<G4int,G4double>*>*a=new vector<pair<G4int,G4double>*>;
    vector<pair<G4int,G4double>*>*s_vec=new vector<pair<G4int,G4double>*>;
    G4double last=0.;
    if(n) for(G4int j=0; j<n; j++)
    {
      G4int    nn =(*is)[j].first; // #ofNeutrons in the isotope
      G4double cur=(*is)[j].second;// value of the summed abundancy
      pair<G4int,G4double>* aP = new pair<G4int,G4double>(nn,cur-last);
      last=cur;                     // Update the summed value
      pair<G4int,G4double>* sP = new pair<G4int,G4double>((*is)[j]);
      a->push_back(aP);
      s_vec->push_back(sP);
#ifdef cdebug
      G4cout<<"G4QIsotope::Constructor:Element# "<<i<<", Pair # "<<j<<" is filled"<<G4endl;
#endif
    }
    natElements.push_back(a);       // Fill abundancies for the particular isotope
    natSumAbund.push_back(s_vec);   // Fill summes abundancies up to this isotope
#ifdef cdebug
    G4cout<<"G4QIsotope::Constructor: natElements is filled"<<G4endl;
#endif
    vector<pair<G4int,G4double>*>*c=new vector<pair<G4int,G4double>*>;
    if(n) for(G4int j=0; j<n; j++)  // Cross sections are 0. by default
    {
      pair<G4int,G4double>* cP = new pair<G4int,G4double>((*is)[j].first,0.);
      c->push_back(cP);
#ifdef cdebug
      G4cout<<"G4QIsotope::Constructor:CrosSecPair i="<<i<<", j="<<j<<" is filled"<<G4endl;
#endif
    }
    natIsoCrosS.push_back(c);       // FillPrototypeCrossSec's (0) for theParticularIsotope
#ifdef cdebug
    G4cout<<"G4QIsotope::Constructor: natIsoCrosS is filled"<<G4endl;
#endif
    delete is;
  }
#ifdef cdebug
  G4cout<<"G4QIsotope::Constructor: is finished"<<G4endl;
#endif
}

G4QIsotope::~G4QIsotope()          // The QIsotopes are destructed only in theEnd of theJob
{
#ifdef debug
  G4cout<<"G4QIsotope::Destructor is called"<<G4endl;
#endif
  G4int uP=natElements.size();
  if(uP) for(G4int i=0; i<uP; i++)
  {
    vector<pair<G4int,G4double>*>* curA=natElements[i];
    G4int nn=curA->size();         // Can not be 0 by definition
    if(nn) for(G4int n=0; n<nn; n++) delete (*curA)[n]; // Delete pair(N,Ab)
    delete curA;                   // Delet abundancy vector
    vector<pair<G4int,G4double>*>* curS=natSumAbund[i];
    G4int ns_value=curS->size();   // Can not be 0 by definition
    if(ns_value) for(G4int n=0; n<ns_value; n++) delete (*curS)[n]; // Delete pair(N,Ab)
    delete curS;                   // Delet abundancy vector
    vector<pair<G4int,G4double>*>* curC=natIsoCrosS[i];
    G4int nc=curC->size();         // Can not be 0 by definition
    if(nc) for(G4int k=0; k<nc; k++) delete (*curC)[k]; // Delete pair(N,CS)
    delete curC;                   // Delete cross section vector
  }
  G4int nP=newElems.size();
  if(nP) for(G4int j=0; j<nP; j++) // LOOP over new UserDefinedElements
  {
    pair<G4int, vector<pair<G4int,G4double>*>* >* nEl= newElems[j];
    G4int nEn=nEl->second->size();
    if(nEn) for(G4int k=0; k<nEn; k++) delete (*(nEl->second))[k]; // Del vect<pair(N,A)*>
    delete nEl->second;            // Delete the vector
    delete nEl;                    // Delete vect<IndZ,vect<pair(N,Ab)*>*> newElementVector
    //
    pair<G4int, vector<pair<G4int,G4double>*>* >* nSA= newSumAb[j];
    G4int nSn=nSA->second->size();
    if(nSn) for(G4int n=0; n<nSn; n++) delete (*(nSA->second))[n]; // Del vect<pair(N,S)*>
    delete nSA->second;            // Delete the vector
    delete nSA;                    // Delete vect<IndZ,vect<pair(N,SA)*>*> newSumAbunVector
    //
    pair<G4int, vector<pair<G4int,G4double>*>* >* nCS= newIsoCS[j];
    G4int nCn=nCS->second->size();
    if(nCn) for(G4int n=0; n<nCn; n++) delete (*(nCS->second))[n]; // Del vect<pair(N,C)*>
    delete nCS->second;            // Delete the vector
    delete nCS;                    // Delete vect<IndZ,vect<pair(N,CS)*>*> newIsoCroSVector
    //
    if(nEn!=nCn) G4cerr<<"*G4QIsotope-WORNING-:#El="<<j<<":nE="<<nEn<<"!=nC="<<nCn<<G4endl;
    if(nEn!=nSn) G4cerr<<"*G4QIsotope-WORNING-:#El="<<j<<":nE="<<nEn<<"!=nS="<<nSn<<G4endl;
  }
}

// Returns Pointer to the G4QIsotope singletone
G4QIsotope* G4QIsotope::Get()
{
#ifdef pdebug
  G4cout<<"G4QIsotope::Get is called"<<G4endl;
#endif
  static G4QIsotope theIsotopes;             // *** Static body of the G4QIsotope class ***
  return &theIsotopes;
}

// #ofProtons in stable isotopes with fixed A=Z+N. Returns length and fils VectOfIsotopes
G4int G4QIsotope::GetProtons(G4int A, vector<G4int>& isoV) 
{
  const G4int nAZ=270;  // Dimension of the table
  // Best Z for the given A
  const G4int bestZ[nAZ] = {
     0,  1,  1,  2,  2,  0,  3,  3,  4,  4,   //0
     5,  5,  6,  6,  7,  7,  8,  8,  8,  9,   //10
    10, 10, 10, 11, 12, 12, 12, 13, 14, 14,   //20
    14, 15, 16, 16, 16, 17, 18, 17, 18, 19,   //30
    18, 19, 20, 20, 20, 21, 22, 22, 23, 23,   //40
    22, 23, 24, 24, 26, 25, 26, 26, 28, 27,   //50
    28, 28, 28, 29, 30, 29, 30, 30, 30, 31,   //60
    32, 31, 32, 32, 32, 33, 34, 34, 34, 35,   //70
    34, 35, 36, 36, 36, 37, 39, 36, 38, 39,   //80
    40, 40, 41, 40, 40, 42, 42, 42, 42, 44,   //90
    44, 44, 44, 45, 44, 46, 46, 47, 46, 47,   //100
    48, 48, 48, 48, 48, 49, 50, 50, 50, 50,   //110
    50, 51, 50, 51, 50, 52, 52, 53, 52, 54,   //120
    52, 54, 54, 55, 54, 56, 54, 56, 56, 57,   //130
    58, 59, 60, 60, 60, 60, 60, 62, 62, 62,   //140
    62, 63, 62, 63, 62, 64, 64, 64, 64, 65,   //150
    64, 66, 66, 66, 66, 67, 68, 68, 68, 69,   //160
    68, 70, 70, 70, 70, 71, 70, 72, 72, 72,   //170
    72, 73, 74, 74, 74, 75, 74, 75, 76, 76,   //180
    76, 77, 76, 77, 78, 78, 78, 79, 80, 80,   //190
    80, 80, 80, 81, 80, 81, 82, 82, 82, 83,   //200
    82,  0, 82,  0, 82,  0, 84,  0,  0,  0,   //210
    86,  0, 86, 87, 88,  0, 88, 89, 88, 89,   //220
    89, 91, 90,  0, 92, 92,  0, 93, 92, 94,   //230
     0,  0,  0, 95, 94,  0,  0, 96,  0,  0,   //240
     0, 98, 99,  0,  0,  0,  0,100,101,102,   //250
   103,104,105,106,  0,108,109,  0,  0,  0};  //260
  // 0   1   2   3   4   5   6   7   8   9
  // Second candidate
  const G4int secoZ[nAZ] = {
     0,  0,  0,  1,  0,  0,  0,  4,  0,  0,   //0
     4,  6,  5,  7,  0,  8,  0,  0,  0,  8,   //10
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   //20
     0,  0, 15,  0,  0, 16,  0,  0,  0,  0,   //30
    20, 20,  0,  0,  0,  0, 20,  0,  0,  0,   //40
    24,  0,  0,  0, 24,  0,  0,  0, 26,  0,   //50
    27,  0,  0,  0, 28,  0,  0,  0,  0,  0,   //60
    30,  0,  0,  0, 34,  0, 32,  0,  0,  0,   //70
    36,  0, 34,  0, 38,  0, 38, 38,  0,  0,   //80
     0,  0, 42,  0, 42,  0, 44,  0, 44,  0,   //90
    42,  0, 46,  0, 46,  0, 48,  0, 48,  0,   //100
    46,  0, 50, 49, 50, 50, 48,  0,  0,  0,   //110
    52,  0, 52,  0, 52,  0, 54,  0, 54,  0,   //120
    54,  0, 56,  0, 56,  0, 56,  0, 58,  0,   //130
    54,  0, 58,  0, 62, 61,  0,  0, 60,  0,   //140
    60,  0, 64,  0, 64,  0, 66,  0, 66,  0,   //150
    66,  0, 68,  0, 68,  0,  0,  0, 70,  0,   //160
    70,  0,  0,  0, 72,  0, 72,  0,  0,  0,   //170
    74,  0,  0,  0, 76,  0, 76, 76,  0,  0,   //180
    78,  0, 78,  0,  0,  0, 80,  0, 78,  0,   //190
     0,  0,  0,  0, 82,  0,  0,  0,  0, 84,   //200
    84,  0, 83,  0, 83,  0,  0,  0,  0,  0,   //210
     0,  0,  0,  0,  0,  0,  0,  0, 89,  0,   //220
     0,  0,  0,  0, 93,  0,  0,  0, 93,  0,   //230
     0,  0,  0,  0,  0,  0,  0, 97,  0,  0,   //240
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   //250
     0,  0,107,  0,  0,  0,  0,  0,  0,  0};  //260
  // 0   1   2   3   4   5   6   7   8   9
  // Third candidate
  const G4int thrdZ[nAZ] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //0
    0, 0, 7, 0, 0, 0, 0, 0, 0, 0,   //10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //20
    0, 0, 0, 0, 0,20, 0, 0, 0, 0,   //30
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //40
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //50
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //60
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //70
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //80
    0, 0,36, 0,38, 0,40, 0, 0, 0,   //90
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //100
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //110
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //120
   56, 0,50, 0, 0, 0,58, 0,57, 0,   //130
    0, 0, 0, 0, 0, 0, 0, 0,65, 0,   //140
    0, 0,66, 0, 0, 0, 0, 0, 0, 0,   //150
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //160
    0, 0, 0, 0, 0, 0,71, 0, 0, 0,   //170
   73, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //180
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //190
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //200
   83, 0,84, 0,84, 0, 0, 0, 0, 0,   //210
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //220
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //230
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //240
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //250
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //260
  //0  1  2  3  4  5  6  7  8  9
  // Fourth candidate (only two isotopes)
  const G4int quadZ[nAZ] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //20
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //30
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //40
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //50
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //60
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //70
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //80
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //90
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //100
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //110
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //120
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //130
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //140
    0, 0,67, 0, 0, 0, 0, 0, 0, 0,   //150
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //160
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //170
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //180
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //190
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //200
   85, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //210
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //220
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //230
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //240
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //250
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //260
  //0  1  2  3  4  5  6  7  8  9
  isoV.clear();                     // Always cleans up before filling
  if(A>=nAZ) return 0;
  G4int fZ=bestZ[A];
  if(fZ)
  {
    isoV.push_back(fZ);
    G4int sZ=secoZ[A];
    if(sZ)
    {
      isoV.push_back(sZ);
      G4int tZ=thrdZ[A];
      if(tZ)
      {
        isoV.push_back(tZ);
        G4int qZ=quadZ[A];
        if(qZ)
        {
          isoV.push_back(qZ);
          return 4;
        }
        else return 3;
      }
      else return 2;
    }
    else return 1;
  }
  else return 0;
}

// Randomize a#ofNeutrons in the Element (independent function @@ can be used for initial)
G4int G4QIsotope::RandomizeNeutrons(G4int i) 
{
  static const G4int nElements = 110; // Max=Meitnerium(Mt)Z=99(starts with Z=0 - neuteron)
  // If an error is found in this simple tree, please, correct the initialization above
  static G4int dN[nElements]=
  {
 //   n              Be                   F              Al       P
      1, -1, -1, -1,  5, -1, -1, -1, -1, 10, -1, 12, -1, 14, -1, 16, -1, -1, -1, -1,
 //      Sc              Mn      Co                      As                       Y
     -1, 24, -1, -1, -1, 30, -1, 32, -1, -1, -1, -1, -1, 42, -1, -1, -1, -1, -1, 50,
 //      Nb      Tc      Rh                               I      Cs              Pr
     -1, 52, -1, 55, -1, 58, -1, -1, -1, -1, -1, -1, -1, 74, -1, 78, -1, -1, -1, 82,
 //      Pm              Tb      Ho      Tm                                      Au
     -1,-85, -1, -1, -1, 94, -1, 98, -1,100, -1, -1, -1, -1, -1, -1, -1, -1, -1,118,
 //                      Po   At   Rn   Fr   Ra   Ac
     -1,  -1,  -1, 126,-125,-136,-136,-138,-138,-142,
 //  Th   Pa    U   Np   Pu   Am   Cm   Bk   Cf   Es
    142,-140,  -1,-144,-150,-148,-151,-150,-153,-153
 //  Fm   Md   No   Lr   Rf   Db   Sg   Bh   Hs   Mt
   -157,-157,-157,-157,-157,-157,-157,-155,-157,-157
  };
#ifdef debug
  G4cout<<"G4QIsotope::RandomizeNeutrons is called Z="<<i<<G4endl;
#endif
  G4int N=dN[i];
  if(N==-1)
  {
    G4double rnd=G4UniformRand();
    if          (i<44)        // =----= H - Mo
    {
      if        (i<23)        // ------ H - Ti
      {
        if      (i<12)        // ______ H - Ne
        {
          if     (i<6)        // ...... H - B
          {
            if   (i<3)
            {
              if(i==1)        // H
              {
                if(rnd>.00015)       N=0;
                else                 N=1;
              }
              else            // He (2)
              {
                if(rnd>1.37e-6)      N=2;
                else                 N=1;
              }
            }
            else
            {
              if(i==3)        // Li
              {
                if(rnd>.075)         N=4;
                else                 N=3;
              }
              else            // B (5)
              {
                if(rnd>.199)         N=6;
                else                 N=5;
              }
            }
          }
          else                // ...... C - Ne
          {
            if   (i<8)
            {
              if(i==6)        // C
              {
                if(rnd>.011)         N=6;
                else                 N=7;
              }
              else            // N (7)
              {
                if(rnd>.0037)        N=7;
                else                 N=8;
              }
            }
            else
            {
              if(i==8)        // O
              {
                if     (rnd<.9976)   N=8;
                else if(rnd<.9996)   N=10;
                else                 N=9;
              }
              else            // Ne (10)
              {
                if     (rnd<.9948)   N=10;
                else if(rnd<.9975)   N=11;
                else                 N=12;
              }
            }
          }
        }
        else                  // ______ Mg - Ti
        {
          if     (i<18)       // ...... Mg - Cl
          {
            if   (i<16)
            {
              if(i==12)       // Mg
              {
                if     (rnd<.7899)   N=12;
                else if(rnd<.8899)   N=13;
                else                 N=14;
              }
              else            // Si (14)
              {
                if    (rnd<.9223)    N=14;
                else if(rnd<.969)    N=15;
                else                 N=16;
              }
            }
            else
            {
              if(i==16)       // S
              {
                if     (rnd<.9502)   N=16;
                else if(rnd<.9923)   N=18;
                else if(rnd<.9998)   N=17;
                else                 N=20;
              }
              else            // Cl (17)
              {
                if     (rnd>.7577)   N=18;
                else                 N=20;
              }
            }
          }
          else                // ...... Ar - Ti
          {
            if   (i<20)
            {
              if(i==18)       // Ar
              {
                if     (rnd<.996)    N=22;
                else if(rnd<.99937)  N=18;
                else                 N=20;
              }
              else            // K (19)
              {
                if     (rnd<.932581) N=20;
                else if(rnd<.999883) N=22;
                else                 N=21;
              }
            }
            else
            {
              if(i==20)       // Ca
              {
                if     (rnd<.96941)  N=20;
                else if(rnd<.99027)  N=24;
                else if(rnd<.99674)  N=22;
                else if(rnd<.99861)  N=28;
                else if(rnd<.99996)  N=23;
                else                 N=26;
              }
              else            // Ti (22)
              {
                if     (rnd<.738)    N=26;
                else if(rnd<.818)    N=24;
                else if(rnd<.891)    N=25;
                else if(rnd<.946)    N=27;
                else                 N=28;
              }
            }
          }
        }
      }
      else                    // ------ V - Mo
      {
        if     (i<32)         // ______ V - Ga
        {
          if     (i<28)       // ...... H - Fe
          {
            if   (i<26)
            {
              if(i==23)       // V
              {
                if     (rnd<.9975)   N=28;
                else                 N=27;
              }
              else            // Cr (24)
              {
                if     (rnd<.8379)   N=28;
                else if(rnd<.9329)   N=29;
                else if(rnd<.97635)  N=26;
                else                 N=30;
              }
            }
            else              // Fe (26)
            {
                if     (rnd<.9172)   N=30;
                else if(rnd<.9762)   N=28;
                else if(rnd<.9972)   N=31;
                else                 N=32;
            }
          }
          else                // ...... Ni - Ga
          {
            if   (i<30)
            {
              if(i==28)       // Ni
              {
                if     (rnd<.68077)  N=30;
                else if(rnd<.943)    N=32;
                else if(rnd<.97934)  N=34;
                else if(rnd<.99074)  N=33;
                else                 N=36;
              }
              else            // Cu (29)
              {
                if     (rnd<.6917)   N=34;
                else                 N=36;
              }
            }
            else
            {
              if(i==30)       // Zn
              {
                if     (rnd<.486)    N=34;
                else if(rnd<.765)    N=36;
                else if(rnd<.953)    N=38;
                else if(rnd<.994)    N=37;
                else                 N=40;
              }
              else            // Ga (31)
              {
                if     (rnd<.60108)  N=38;
                else                 N=40;
              }
            }
          }
        }
        else                  // ______ Ge - Mo
        {
          if     (i<37)       // ...... H - B
          {
            if   (i<35)
            {
              if(i==32)       // Ge
              {
                if     (rnd<.3594)  N=42;
                else if(rnd<.6360)  N=40;
                else if(rnd<.8484)  N=38;
                else if(rnd<.9256)  N=41;
                else                N=44;
              }
              else            // Se (34)
              {
                if     (rnd>.4961)  N=46;
                else if(rnd<.7378)  N=44;
                else if(rnd<.8274)  N=42;
                else if(rnd<.9148)  N=48;
                else if(rnd<.9911)  N=43;
                else                N=40;
              }
            }
            else
            {
              if(i==35)       // Br
              {
                if     (rnd<.5069)  N=44;
                else                N=46;
              }
              else            // Kr (36)
              {
                if     (rnd<.57)    N=48;
                else if(rnd<.743)   N=50;
                else if(rnd<.859)   N=46;
                else if(rnd<.974)   N=47;
                else if(rnd<.9965)  N=44;
                else                N=42;
              }
            }
          }
          else                // ...... Rb - Mo
          {
            if     (i<40)
            {
              if(i==37)       // Rb
              {
                if     (rnd<.7217)  N=48;
                else                N=50;
              }
              else            // SR (38)
              {
                if     (rnd<.8258)  N=50;
                else if(rnd<.9244)  N=48;
                else if(rnd<.9944)  N=49;
                else                N=46;
              }
            }
            else
            {
              if(i==40)       // Zr
              {
                if     (rnd<.5145)  N=50;
                else if(rnd<.6883)  N=54;
                else if(rnd<.8598)  N=53;
                else if(rnd<.972)   N=51;
                else                N=56;
              }
              else            // Mo (42)
              {
                if     (rnd<.2413)  N=56;
                else if(rnd<.4081)  N=54;
                else if(rnd<.5673)  N=53;
                else if(rnd<.7157)  N=50;
                else if(rnd<.8120)  N=58;
                else if(rnd<.9075)  N=55;
                else                N=52;
              }
            }
          }
        }
      }
    }
    else                      // =----= Ru - U
    {
      if         (i<66)       // ------ Ru - Gd
      {
        if       (i<54)       // ______ Ru - Te
        {
          if     (i<49)       // ...... Ru - Cd
          {
            if   (i<47)
            {
              if(i==44)       // Ru
              {
                if     (rnd<.316)   N=58;
                else if(rnd<.502)   N=60;
                else if(rnd<.673)   N=57;
                else if(rnd<.8)     N=55;
                else if(rnd<.926)   N=56;
                else if(rnd<.9814)  N=52;
                else                N=54;
              }
              else            // Pd (46)
              {
                if     (rnd<.2733)  N=60;
                else if(rnd<.5379)  N=62;
                else if(rnd<.7612)  N=59;
                else if(rnd<.8784)  N=55;
                else if(rnd<.9898)  N=58;
                else                N=56;
              }
            }
            else
            {
              if(i==47)       // Ag
              {
                if(rnd<.51839)      N=60;
                else                N=62;
              }
              else            // Cd (48)
              {
                if     (rnd<.2873)  N=66;
                else if(rnd<.5286)  N=64;
                else if(rnd<.6566)  N=59;
                else if(rnd<.7815)  N=62;
                else if(rnd<.9037)  N=65;
                else if(rnd<.9786)  N=68;
                else if(rnd<.9911)  N=58;
                else                N=60;
              }
            }
          }
          else                // ...... In - Te
          {
            if   (i<51)
            {
              if(i==49)       // In
              {
                if     (rnd<.9577)  N=66;
                else                N=64;
              }
              else            // Sn (50)
              {
                if     (rnd<.3259)  N=70;
                else if(rnd<.5681)  N=68;
                else if(rnd<.7134)  N=66;
                else if(rnd<.7992)  N=69;
                else if(rnd<.8760)  N=67;
                else if(rnd<.9339)  N=74;
                else if(rnd<.9802)  N=72;
                else if(rnd<.9899)  N=62;
                else                N=64;
                //else if(rnd<.9964)  N=64;
                //else                N=65;
              }
            }
            else
            {
              if(i==51)       // Sb
              {
                if     (rnd<.5736)  N=70;
                else                N=72;
              }
              else            // Te (52)
              {
                if     (rnd<.3387)  N=78;
                else if(rnd<.6557)  N=76;
                else if(rnd<.8450)  N=74;
                else if(rnd<.9162)  N=73;
                else if(rnd<.9641)  N=72;
                else if(rnd<.9900)  N=70;
                else if(rnd<.99905) N=71;
                else                N=68;
              }
            }
          }
        }
        else                // ______ Xe - Gd
        {
          if     (i<60)     // ...... Xe - B
          {
            if   (i<57)
            {
              if(i==54)       // Xe
              {
                if     (rnd<.269)   N=78;
                else if(rnd<.533)   N=75;
                else if(rnd<.745)   N=77;
                else if(rnd<.849)   N=80;
                else if(rnd<.938)   N=82;
                else if(rnd<.979)   N=76;
                else if(rnd<.9981)  N=74;
                else if(rnd<.9991)  N=70;
                else                N=72;
              }
              else            // Ba (56)
              {
                if     (rnd<.717)   N=82;
                else if(rnd<.8293)  N=81;
                else if(rnd<.9078)  N=80;
                else if(rnd<.97373) N=79;
                else if(rnd<.99793) N=78;
                else if(rnd<.99899) N=74;
                else                N=76;
              }
            }
            else
            {
              if(i==57)       // La
              {
                if     (rnd<.999098)N=82;
                else                N=81;
              }
              else            // Ce (58)
              {
                if     (rnd<.8843)  N=82;
                else if(rnd<.9956)  N=84;
                else if(rnd<.9981)  N=80;
                else                N=78;
              }
            }
          }
          else                // ...... Nd - Gd
          {
            if   (i<63)
            {
              if(i==60)       // Nd
              {
                if     (rnd<.2713)  N=82;
                else if(rnd<.5093)  N=84;
                else if(rnd<.6812)  N=86;
                else if(rnd<.8030)  N=83;
                else if(rnd<.8860)  N=85;
                else if(rnd<.9436)  N=88;
                else                N=90;
              }
              else            // Sm (62)
              {
                if     (rnd<.267)   N=90;
                else if(rnd<.494)   N=92;
                else if(rnd<.644)   N=85;
                else if(rnd<.782)   N=87;
                else if(rnd<.895)   N=86;
                else if(rnd<.969)   N=88;
                else                N=82;
              }
            }
            else
            {
              if(i==63)       // Eu
              {
                if     (rnd<.522)   N=90;
                else                N=89;
              }
              else            // Gd (64)
              {
                if     (rnd<.2484)  N=94;
                else if(rnd<.4670)  N=96;
                else if(rnd<.6717)  N=92;
                else if(rnd<.8282)  N=93;
                else if(rnd<.9762)  N=91;
                else if(rnd<.9980)  N=90;
                else                N=88;
              }
            }
          }
        }
      }
      else                    // ------ Dy - U
      {
        if       (i<76)       // ______ Dy - Re
        {
          if     (i<72)       // ...... Dy - Lu
          {
            if   (i<70)
            {
              if(i==66)       // Dy
              {
                if     (rnd<.282)   N=98;
                else if(rnd<.537)   N=96;
                else if(rnd<.786)   N=97;
                else if(rnd<.975)   N=95;
                else if(rnd<.9984)  N=94;
                else if(rnd<.9994)  N=92;
                else                N=90;
              }
              else            // Er (68)
              {
                if     (rnd<.3360)  N= 98;
                else if(rnd<.6040)  N=100;
                else if(rnd<.8335)  N= 99;
                else if(rnd<.9825)  N=102;
                else if(rnd<.9986)  N= 96;
                else                N= 94;
              }
            }
            else
            {
              if(i==70)       // Yb
              {
                if     (rnd<.3180)  N=104;
                else if(rnd<.5370)  N=102;
                else if(rnd<.6982)  N=103;
                else if(rnd<.8412)  N=101;
                else if(rnd<.9682)  N=106;
                else if(rnd<.9987)  N=100;
                else                N= 98;
              }
              else            // Lu (71)
              {
                if     (rnd<.9741)  N=104;
                else                N=105;
              }
            }
          }
          else                // ...... Hf - Re
          {
            if   (i<74)
            {
              if(i==72)       // Hf
              {
                if     (rnd<.35100) N=108;
                else if(rnd<.62397) N=106;
                else if(rnd<.81003) N=105;
                else if(rnd<.94632) N=107;
                else if(rnd<.99838) N=104;
                else                N=102;
              }
              else            // Ta (73)
              {
                if(rnd<.99988) N=108;
                else           N=107;
              }
            }
            else
            {
              if(i==74)       // W
              {
                if     (rnd<.307)   N=110;
                else if(rnd<.593)   N=112;
                else if(rnd<.856)   N=108;
                else if(rnd<.9988)  N=109;
                else                N=106;
              }
              else            // Re (75)
              {
                if     (rnd<.626)   N=112;
                else                N=110;
              }
            }
          }
        }
        else                  // ______ Os - U
        {
          if     (i<81)       // ...... Os - Hg
          {
            if     (i<78)
            {
              if(i==76)       // Os
              {
                if     (rnd<.410)   N=116;
                else if(rnd<.674)   N=114;
                else if(rnd<.835)   N=113;
                else if(rnd<.968)   N=112;
                else if(rnd<.984)   N=111;
                else if(rnd<.9998)  N=110;
                else                N=108;
              }
              else            // Ir (77)
              {
                if     (rnd<.627)   N=116;
                else                N=114;
              }
            }
            else
            {
              if(i==78)       // Pt
              {
                if     (rnd<.338)   N=117;
                else if(rnd<.667)   N=116;
                else if(rnd<.920)   N=118;
                else if(rnd<.992)   N=120;
                else if(rnd<.9999)  N=114;
                else                N=112;
              }
              else            // Hg (80)
              {
                if     (rnd<.2986)  N=122;
                else if(rnd<.5296)  N=120;
                else if(rnd<.6983)  N=119;
                else if(rnd<.8301)  N=121;
                else if(rnd<.9298)  N=118;
                else if(rnd<.9985)  N=124;
                else                N=116;
              }
            }
          }
          else                // ...... Tl - U
          {
            if        (i<92)
            {
              if     (i==81)  // Tl
              {
                if     (rnd<.70476) N=124;
                else                N=122;
              }
              else            // Pb (82)
              {
                if     (rnd<.524)   N=126;
                else if(rnd<.765)   N=124;
                else if(rnd<.986)   N=125;
                else                N=122;
              }
            }
            else              // U (92)
            {
                if     (rnd<.992745)N=146;
                else if(rnd<.999945)N=143;
                else                N=142;
            }
          }
        }
      }
    }
  }
  else if(N<0)
  {
    N=-N;
    G4cerr<<"---Worn---G4QIsotope::RandomizeNeutrons:UnstableElement's used Z="<<i<<G4endl;
  }
#ifdef debug
  G4cout<<"G4QIsotope::RandomizeNeutrons: Z="<<i<<", N="<<N<<G4endl;
#endif
  return N;
}

// Returns the input index (if it is >0 & unique) or theFirstFreeIndex (<=0 or nonunique)
G4int G4QIsotope::InitElement(G4int Z, G4int index, // Ret: -1 - Empty, -2 - Wrong (sum>1)
                              vector<pair<G4int,G4double>*>* abund)
{
  G4int I=abund->size();
#ifdef debug
  G4cout<<"G4QIsotope::InitElement: called with I="<<I<<" pairs,Z="<<Z<<",i="<<indexG4endl;
#endif
  if(I<=0)
  {
    G4cerr<<"--Worning--G4QIsotope::InitEl:(-1)0VectorOfNewEl,Z="<<Z<<",i="<<index<<G4endl;
    return -2;
  }
  if(IsDefined(Z,index))              // This index is already defined
  {
    G4cerr<<"-Worning-G4QIsotope::InitEl:VONewEl,Z="<<Z<<",ind="<<index<<" exists"<<G4endl;
    return index;
  }
  G4int ZInd=1000*index+Z;            // Fake Z increased by the UserDefinedIndex
  vector<pair<G4int,G4double>*>*A = new vector<pair<G4int,G4double>*>;
  vector<pair<G4int,G4double>*>*S = new vector<pair<G4int,G4double>*>;
  vector<pair<G4int,G4double>*>*C = new vector<pair<G4int,G4double>*>;
#ifdef debug
  G4cout<<"G4QIsotope::InitElement: A & S & C vectors are alocated"<<G4endl;
#endif
  G4double sumAbu=0;                  // Summ of abbundancies
  for(G4int j=0; j<I; j++)
  {
    G4int N=(*abund)[j]->first;
    G4double abu=(*abund)[j]->second;
#ifdef debug
    G4cout<<"G4QIsotope::InitElement: pair#"<<j<<", N="<<N<<", abund="<<abu<<G4endl;
#endif
    sumAbu+=abu;
    if(j==I-1.)
    {
      if(fabs(sumAbu-1.)>.00001)
      {
        G4cerr<<"--Worning--G4QIsotope::InitEl:maxSum="<<sumAbu<<" is fixed to 1."<<G4endl;
        abu+=1.-sumAbu;
        sumAbu=1.;
      }
      else if(sumAbu-abu>1.)
      {
        G4cerr<<"--Worning--G4QIsotope::InitEl:(-2)WrongAbund,Z="<<Z<<",i="<<index<<G4endl;
        for(G4int k=0; k<I-1; k++)
        {
          delete (*A)[k];
          delete (*S)[k];
          delete (*C)[k];
        }
        delete A; 
        delete S;
        delete C;
        return -2;
      }
#ifdef debug
    G4cout<<"G4QIsotope::InitElement:TheLastIsChecked it isOK or coredTo "<<sumAbu<<G4endl;
#endif
    }
    pair<G4int,G4double>* abP= new pair<G4int,G4double>(N,abu);
    A->push_back(abP); // @@ Valgrind thinks that it is not deleted (?) (Line 703)
    pair<G4int,G4double>* saP= new pair<G4int,G4double>(N,sumAbu);
    S->push_back(saP); // @@ Valgrind thinks that it is not deleted (?) (Line 713)
    pair<G4int,G4double>* csP= new pair<G4int,G4double>(N,0.);
    C->push_back(csP); // @@ Valgrind thinks that it is not deleted (?) (Line 723)
#ifdef debug
    G4cout<<"G4QIsotope::InitElement: A & S & C are filled nP="<<C->size()<<G4endl;
#endif
  }
  pair<G4int,vector<pair<G4int,G4double>*>*>* newAP=
                                    new pair<G4int,vector<pair<G4int,G4double>*>*>(ZInd,A);
  newElems.push_back(newAP);
  pair<G4int,vector<pair<G4int,G4double>*>*>* newSA=
                                    new pair<G4int,vector<pair<G4int,G4double>*>*>(ZInd,S);
  newSumAb.push_back(newSA);
  pair<G4int,vector<pair<G4int,G4double>*>*>* newCP=
                                    new pair<G4int,vector<pair<G4int,G4double>*>*>(ZInd,C);
  newIsoCS.push_back(newCP);
#ifdef debug
  G4cout<<"G4QIsotope::InitElement: newElems & newSumAb & newIsoCS are filled "<<G4endl;
#endif
  return index;
}

// The highest index defined for Element with Z (Index>0 correspondToUserDefinedElements)
G4int G4QIsotope::GetLastIndex(G4int Z) // Returns theLastDefinedIndex (if onlyNatural: =0)
{
#ifdef debug
  G4cout<<"G4QIsotope::GetLastIndex is called Z="<<Z<<G4endl;
#endif
  G4int mind=0;                       // Prototype of the maximum existing index for this Z
  G4int nE=newElems.size();           // A number of definitions in the newElements vector
  if(nE) for(G4int i=0; i<nE; i++)    // LOOP over new UserDefinedElements
  {
    G4int zin=newElems[i]->first;
    G4int zi=zin%1000;                // Existing Z
    G4int in=zin/1000;                // Existing index
    if(Z==zi && in>mind) mind=in;     // maximum index for this Z
  }
  return mind;
}

// Indices can have differen numbers (not 1,2,3,...) & in different sequences (9,3,7,...)
G4bool G4QIsotope::IsDefined(G4int Z, G4int Ind) // Ind is an index to be found (true)
{
#ifdef debug
  G4cout<<"G4QIsotope::IsDefined is called Z="<<Z<<", I="<<Ind<<G4endl;
#endif
  if(Ind<=0)
  {
    if(Ind<0) G4cerr<<"-W-G4QIsotope::IsDefined: Z="<<Z<<", Ind="<<Ind<<" < 0->=0"<<G4endl;
    return true;                      // to avoid definition with the negative index
  }
  G4int nE=newElems.size();           // A number of definitions in the newElements vector
  if(nE) for(G4int i=0; i<nE; i++)    // LOOP over new UserDefinedElements
  {
    G4int zin=newElems[i]->first;
    G4int zi=zin%1000;                // Existing Z
    G4int in=zin/1000;                // Existing index
    if(Z==zi && Ind==in) return true;  // The index for the element Z is found
  }
  return false;                       // The index for the element Z is not found
}

// A#ofNeutrons in theElement with Z & UserDefIndex. Universal for Nat(index=0) & UserDefEl
G4int G4QIsotope::GetNeutrons(G4int Z, G4int index) // If theElem doesn't exist, returns <0
{
#ifdef debug
  G4cout<<"G4QIsotope::GetNeutrons is called Z="<<Z<<", index="<<index<<G4endl;
#endif
  // To reduce the code, but make the member function a bit slower, one can use for natural
  // isotopes the same algorithm as for the newElements, splitting the natElements Vector
  if(!index) return RandomizeNeutrons(Z); // @@ Fast decision for the natural isotopes
  else if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetNeutrons:(-2) Negative Index i="<<index<<G4endl;
    return -2;
  }
  // For the positive index tries to randomize the newUserDefinedElement
  G4bool found=false;                 // Prototype of the"ZWithTheSameIndex is found" event
  G4int nE=newElems.size();           // A number of definitions in the newElements Vector
  G4int i=0;
  if(nE) for(i=0; i<nE; i++)
  {
    G4int zin=newElems[i]->first;
    G4int in=zin/1000;                // Existing index
    G4int zi=zin%1000;                // Existing Z
    if(Z==zi && in==index)
    {
      found=true;                     // The newElement with the same Z & index is found
      break;                          // Finish the search and quit the loop
    }
  }
  if(!found)
  {
    G4cerr<<"--Worning--G4QIsotope::GetNeutrons:(-1) NotFound Z="<<Z<<",i="<<index<<G4endl;
    return -1;
  }
  vector<pair<G4int,G4double>*>* abu = newSumAb[i]->second;
  G4int nn = abu->size();             // A#Of UserDefinedIsotopes for the newElement
  if(nn>0)
  {
    if(nn==1) return (*abu)[0]->first;
    else
    {
      G4double rnd=G4UniformRand();
      G4int j=0;
      for(j=0; j<nn; j++) if ( rnd < (*abu)[j]->second ) break;
      if(j>=nn) j=nn-1;
      return (*abu)[j]->first;
    }
  }
  else
  {
    G4cerr<<"--Worning--G4QIsotope::GetNeutrons:(-3) Empty Z="<<Z<<",i="<<index<<G4endl;
    return -3;
  }
}

// Get a pointer to the vector of pairs(N,CrosSec), where N is used to calculate CrosSec
vector<pair<G4int,G4double>*>* G4QIsotope::GetCSVector(G4int Z, G4int index)
{
#ifdef debug
  G4cout<<"G4QIsotope::GetCSVector is called"<<G4endl;
#endif
  if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetSCVector:(-1) Negative Index i="<<index<<G4endl;
    return 0;
  }
  else if(!index) return natIsoCrosS[Z];
  // For the positive index tries to find the newUserDefinedElement
  G4bool found=false;                 // Prototype of the"ZWithTheSameIndex is found" event
  G4int nE=newIsoCS.size();           // A number of definitions in the newElements Vector
  G4int i=0;
  if(nE) for(i=0; i<nE; i++)
  {
    G4int zin=newIsoCS[i]->first;
    G4int in=zin/1000;                // Existing index
    G4int zi=zin%1000;                // Existing Z
    if(Z==zi && in==index)
    {
      found=true;                     // The newElement with the same Z & index is found
      break;                          // Finish the search and quit the loop
    }
  }
  if(!found)
  {
    G4cerr<<"--Worning--G4QIsotope::GetSCVector:(-2) NotFound Z="<<Z<<",i="<<index<<G4endl;
    return 0;
  }
  return newIsoCS[i]->second;
}

// Get a pointer to the vector of pairs(N,IntAbundancy) for the element with Z
vector<pair<G4int,G4double>*>* G4QIsotope::GetAbuVector(G4int Z, G4int index)
{
#ifdef debug
  G4cout<<"G4QIsotope::GetAbuVector is called"<<G4endl;
#endif
  if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetAbuVector:(-1) Negative Index i="<<index<<G4endl;
    return 0;
  }
  else if(!index) return natElements[Z];
  // For the positive index tries to find the newUserDefinedElement
  G4bool found=false;                 // Prototype of the"ZWithTheSameIndex is found" event
  G4int nE=newElems.size();           // A number of definitions in the newElements Vector
  G4int i=0;
  if(nE) for(i=0; i<nE; i++)
  {
    G4int zin=newElems[i]->first;
    G4int in=zin/1000;                // Existing index
    G4int zi=zin%1000;                // Existing Z
    if(Z==zi && in==index)
    {
      found=true;                     // The newElement with the same Z & index is found
      break;                          // Finish the search and quit the loop
    }
  }
  if(!found)
  {
    G4cerr<<"--Worning--G4QIsotope::GetAbuVector:(-2)NotFound Z="<<Z<<",i="<<index<<G4endl;
    return 0;
  }
  return newElems[i]->second;
}

// Get a pointer to the vector of pairs(N,SumAbundancy) for the element with Z
vector<pair<G4int,G4double>*>* G4QIsotope::GetSumAVector(G4int Z, G4int index)
{
#ifdef debug
  G4cout<<"G4QIsotope::GetSumAVector is called"<<G4endl;
#endif
  if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetSumAVector:(-1) Negative Index i="<<index<<G4endl;
    return 0;
  }
  else if(!index) return natSumAbund[Z];
  // For the positive index tries to find the newUserDefinedElement
  G4bool found=false;                 // Prototype of the"ZWithTheSameIndex is found" event
  G4int nE=newSumAb.size();           // A number of definitions in the newElements Vector
  G4int i=0;
  if(nE) for(i=0; i<nE; i++)
  {
    G4int zin=newSumAb[i]->first;
    G4int in=zin/1000;                // Existing index
    G4int zi=zin%1000;                // Existing Z
    if(Z==zi && in==index)
    {
      found=true;                     // The newElement with the same Z & index is found
      break;                          // Finish the search and quit the loop
    }
  }
  if(!found)
  {
    G4cerr<<"-Worning-G4QIsotope::GetSumAVector:(-2)Not Found Z="<<Z<<",i="<<index<<G4endl;
    return 0;
  }
  return newSumAb[i]->second;
}

// Calculates the mean Cross Section for the initialized Element(ind=0 Nat,ind>0 UserDef)
G4double G4QIsotope::GetMeanCrossSection(G4int Z, G4int index)
{
  vector<pair<G4int,G4double>*>* ab;
  vector<pair<G4int,G4double>*>* cs;
#ifdef ppdebug
  G4cout<<"G4QIsotope::GetMeanCrossSection is called"<<G4endl;
#endif
  if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetMeanCS:(-1) Negative Index i="<<index<<G4endl;
    return -1.;
  }
  else if(!index)           // =-------=> Natural Abundancies for Isotopes of the Element
  {
#ifdef ppdebug
    G4cout<<"G4QIsotope::GetMeanCrossSection: Nat Abundance, Z="<<Z<<G4endl;
#endif
    ab=natElements[Z];
    cs=natIsoCrosS[Z];
  }
  else                      // =------=> UserDefinedAbundancies for Isotopes of theElement
  {
#ifdef ppdebug
    G4cout<<"G4QIsotope::GetMeanCrossSection: Art Abund, Z="<<Z<<",ind="<<index<<G4endl;
#endif
    // For the positive index tries to find the newUserDefinedElement
    G4bool found=false;               // Prototype of the"ZWithTheSameIndex is found" event
    G4int nE=newIsoCS.size();         // A number of definitions in the newElements Vector
    G4int i=0;
    if(nE) for(i=0; i<nE; i++)
    {
      G4int zin=newIsoCS[i]->first;
      G4int in=zin/1000;              // Existing index
      G4int zi=zin%1000;              // Existing Z
      if(Z==zi && in==index)
      {
        found=true;                   // The newElement with the same Z & index is found
        break;                        // Finish the search and quit the loop
      }
    }
    if(!found)
    {
      G4cerr<<"--Worning--G4QIsotope::GetMeanCS:(-2) NotFound Z="<<Z<<",i="<<index<<G4endl;
      return -2.;
    }
    ab=newElems[i]->second;
    cs=newIsoCS[i]->second;
  }
  G4int nis=ab->size();
  //G4double last=0.;
  if(!nis)
  {
    G4cerr<<"--Worning--G4QIsotope::GetMeanCS:(-3) Empty Z="<<Z<<",i="<<index<<G4endl;
    return -3.;
  }
  else
  {
    G4double sum=0.;
    for(G4int j=0; j<nis; j++)
    {
      G4double cur=(*ab)[j]->second;
      //G4double abunda=cur-last;
      //last=cur;
#ifdef ppdebug
      G4cout<<"G4QIsot::GetMeanCS:j="<<j<<",ab="<<cur<<",CS="<<(*cs)[j]->second<<G4endl;
#endif
      //sum+=abunda * (*cs)[j]->second;
      sum+=cur * (*cs)[j]->second;
    }
    return sum;
  }
}

// Randomize A#OfNeutrons in the Isotope weighted by theAbubdancies and theCrossSections
G4int G4QIsotope::GetCSNeutrons(G4int Z, G4int index)
{
  vector<pair<G4int,G4double>*>* ab;
  vector<pair<G4int,G4double>*>* cs;
#ifdef debug
  G4cout<<"G4QIsotope::GetCSNeutrons is called"<<G4endl;
#endif
  if(index<0)
  {
    G4cerr<<"---Worning---G4QIsotope::GetCSNeutrons:(-1) Negative Index i="<<index<<G4endl;
    return -1;
  }
  else if(!index)           // =---------=> Natural Abundancies for Isotopes of the Element
  {
    ab=natElements[Z];
    cs=natIsoCrosS[Z];
  }
  else                      // =-------=> UserDefinedAbundancies for Isotopes of theElement
  {
    // For the positive index tries to find the newUserDefinedElement
    G4bool found=false;               // Prototype of the"ZWithTheSameIndex is found" event
    G4int nE=newIsoCS.size();         // A number of definitions in the newElements Vector
    G4int i=0;
    if(nE) for(i=0; i<nE; i++)
    {
      G4int zin=newIsoCS[i]->first;
      G4int in=zin/1000;              // Existing index
      G4int zi=zin%1000;              // Existing Z
      if(Z==zi && in==index)
      {
        found=true;                   // The newElement with the same Z & index is found
        break;                        // Finish the search and quit the loop
      }
    }
    if(!found)
    {
      G4cerr<<"--Worning--G4QIsotope::GetCSNeut:(-2) NotFound Z="<<Z<<",i="<<index<<G4endl;
      return -2;
    }
    ab=newElems[i]->second;
    cs=newIsoCS[i]->second;
  }
  G4int nis=ab->size();
  G4double last=0.;
  if(!nis)
  {
    G4cerr<<"--Worning--G4QIsotope::GetCSNeutrons:(-3) Empty Z="<<Z<<",i="<<index<<G4endl;
    return -3;
  }
  else
  {
    G4double sum=0.;
    vector<G4double> scs(nis);
    for(G4int j=0; j<nis; j++)
    {
      G4double cur=(*ab)[j]->second;
      G4double abunda=cur-last;
      last=cur;
      sum+=abunda * (*cs)[j]->second;;
      scs.push_back(sum);
    }
    G4double rnd=sum*G4UniformRand();
    sum=0;
    G4int k=0;
    if(nis>1) for(k=0; k<nis; k++) if(rnd<scs[k]) break;
    return (*ab)[k]->first;
  }
}
