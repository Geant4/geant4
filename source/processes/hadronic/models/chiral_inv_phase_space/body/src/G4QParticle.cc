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
// $Id: G4QParticle.cc,v 1.36 2010-06-25 14:03:44 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QParticle ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Particles in the CHIPS Model
// -------------------------------------------------------------------
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
// -----------------------------------------------------------------------
// Short description: The G4QParticle is a part of the CHIPS World. It is
// characterized by the quark content, spin, mass, width and a vector of
// the decay channels (G4QDecayCannelVector).
// -----------------------------------------------------------------------
//#define debug
//#define pdebug

#include "G4QParticleVector.hh"

G4QParticle::G4QParticle()
{
#ifdef debug
  G4cout<<"G4QParticle::Constr: Default constructor is called"<<G4endl;
#endif
}

G4QParticle::G4QParticle(G4int thePDG)
{
  aQPDG      = G4QPDGCode(thePDG);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(aQPDG.GetQCode());
}

G4QParticle::G4QParticle(G4bool f, G4int theQCode)
{
  aQPDG      = G4QPDGCode(f,theQCode);
#ifdef debug
  G4cout<<"G4QParticle::Constr: PDG="<<aQPDG.GetPDGCode()<<G4endl;
#endif
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(theQCode);
}

G4QParticle::G4QParticle(const G4QParticle& right)
{
  aQPDG                = right.aQPDG;
  //aDecay (Vector)
  G4int nD             = right.aDecay.size(); 
  if(nD) for(G4int id=0; id<nD; id++)
  {
    G4QDecayChan* curD = new G4QDecayChan(right.aDecay[id]);
    aDecay.push_back(curD);
  }

  aQuarkCont           = right.aQuarkCont;
}

G4QParticle::G4QParticle(G4QParticle* right)
{
  aQPDG                = right->aQPDG;
  //aDecay (Vector)
  G4int nD             = right->aDecay.size();
  if(nD) for(G4int id=0; id<nD; id++)
  {
    G4QDecayChan* curD = new G4QDecayChan(right->aDecay[id]);
    aDecay.push_back(curD);
  }

  aQuarkCont           = right->aQuarkCont;
}

G4QParticle::~G4QParticle() 
{
  G4int nDC=aDecay.size();
  //G4cout<<"G4QParticle::Destructor: Before nDC="<<nDC<<G4endl; // TMP
  if(nDC) std::for_each(aDecay.begin(), aDecay.end(), DeleteQDecayChan());
  //G4cout<<"G4QParticle::Destructor: After"<<G4endl; // TMP
  aDecay.clear();
}

// Assignment operator
const G4QParticle& G4QParticle::operator=(const G4QParticle &right)
{
  if(this != &right)                          // Beware of self assignment
  {
    aQPDG                = right.aQPDG;
    //aDecay (Vector)
    G4int iD             = aDecay.size();
    if(iD) for(G4int jd=0; jd<iD; jd++) delete aDecay[jd];
    aDecay.clear();
    G4int nD             = right.aDecay.size();
    if(nD) for(G4int id=0; id<nD; id++)
    {
      G4QDecayChan* curD = new G4QDecayChan(right.aDecay[id]);
      aDecay.push_back(curD);
    }

    aQuarkCont           = right.aQuarkCont;
  }
  return *this;
}

// Standard output for QParticle
std::ostream& operator<<(std::ostream& lhs, G4QParticle& rhs)
//       =========================================
{
  G4QPDGCode rhsQPDG = rhs.GetQPDG();
  lhs << G4endl << "Particle with PDG=" << rhsQPDG << ", Spin=" << rhs.GetSpin()
      << ", mass=" << rhs.GetMass() << ", width=" << rhs.GetWidth() << G4endl;
  lhs<<" Quark Content of the Particle="<<rhs.GetQContent()<<", Decay Channels:"<<G4endl;
  G4QDecayChanVector DCV = rhs.GetDecayVector();
  G4int n = DCV.size();
  for (int i=0; i<n; i++)
  {
    lhs << DCV[i]->GetDecayChanLimit() << "PDG codes";
    G4QPDGCodeVector PCV=DCV[i]->GetVecOfSecHadrons();
    G4int m = PCV.size();
    for (int j=0; j<m; j++)
    {
      if(!j) lhs << ":";
      else   lhs << ",";
      lhs << PCV[j]->GetPDGCode() ;
    }
  }
  return lhs;
}

// Initialize the PDG-Particle by QCode @@ Can be improved, using PDG.DATA file
G4QDecayChanVector G4QParticle::InitDecayVector(G4int nQ)
//    ===================================================
{
  //static G4int nP = 486;                  // Up to A=80
  //static const G4int nP = 494;              // Up to A=80(?) "Isonuclear revision"
  static const G4int nP = 512; // A<57 "Leptons/Hypernuclei" G4QCHIPSWorld::GetParticles(!)
  //static G4QDecayChanVector* DecayDB = new G4QDecayChanVector[nP];
  static G4QDecayChanVector DecayDB[nP];
  static int limit= 0;
  if(nQ>=limit && nQ<nP)
  {
    //*** Secondary PDG-particles should be ordered in a Channel by increasing width***!***
    //** Channels should be ordered by increasing minimum mass of the secondary particles**
    //if(limit<=  0 && nQ>=  0)DecayDB[  0] = 0;    // e-     (11)
    //if(limit<=  1 && nQ>=  1)DecayDB[  1] = 0;    // nu_e   (12)
    //if(limit<=  2 && nQ>=  2)DecayDB[  2] = 0;    // mu-    (13)
    //if(limit<=  3 && nQ>=  3)DecayDB[  3] = 0;    // mu_e   (14)
    //if(limit<=  4 && nQ>=  4)DecayDB[  4] = 0;    // tau-   (15)
    //if(limit<=  5 && nQ>=  5)DecayDB[  5] = 0;    // nu_tau (16)
    //if(limit<=  6 && nQ>=  6)DecayDB[  6] = 0;    // gamma  (22)
    if(limit<=  7 && nQ>=  7)                       // Z0     (23)
    {
      DecayDB[  7].push_back(new G4QDecayChan(.036,  11, -11));
      DecayDB[  7].push_back(new G4QDecayChan(.073,  13, -13));
      DecayDB[  7].push_back(new G4QDecayChan(.107,  15, -15));
      DecayDB[  7].push_back(new G4QDecayChan(.174,  12, -12)); // @@ Fake invisible decays
      DecayDB[  7].push_back(new G4QDecayChan(.240,  14, -14));
      DecayDB[  7].push_back(new G4QDecayChan(.307,  16, -16));
      DecayDB[  7].push_back(new G4QDecayChan(.400,2112,-2112)); // @@ Fake Hadronic decays
      DecayDB[  7].push_back(new G4QDecayChan(.500,2212,-2212)); // @@ Need heavy quarks
      DecayDB[  7].push_back(new G4QDecayChan(.600,2212,-2212, 111));
      DecayDB[  7].push_back(new G4QDecayChan(.700,2112,-2112, 111));
      DecayDB[  7].push_back(new G4QDecayChan(.800,2212,-2112,-211));
      DecayDB[  7].push_back(new G4QDecayChan(.990,2112,-2212, 211));
      DecayDB[  7].push_back(new G4QDecayChan(1.00,2112,-3122, 111));
    }
    if(limit<=  8 && nQ>=  8)                       // W-     (24) @@ Update HadronicDecays
    {
      DecayDB[  8].push_back(new G4QDecayChan(.107,  11, -12));
      DecayDB[  8].push_back(new G4QDecayChan(.214,  13, -14));
      DecayDB[  8].push_back(new G4QDecayChan(.321,  15, -16));
      DecayDB[  8].push_back(new G4QDecayChan(.421,2112,-2212)); // @@ Fake Hadronic decays
      DecayDB[  8].push_back(new G4QDecayChan(.521,2112,-2112,-211));
      DecayDB[  8].push_back(new G4QDecayChan(.621,2212,-2212,-211));
      DecayDB[  8].push_back(new G4QDecayChan(.721,3122,-3122,-211));
      DecayDB[  8].push_back(new G4QDecayChan(.821,2112,-2212, 111));
      DecayDB[  8].push_back(new G4QDecayChan(.921,3122,-2212, 111));
      DecayDB[  8].push_back(new G4QDecayChan(1.00,2112,-3122,-211));
    }
    //if(limit<=  9 && nQ>=  9)DecayDB[  9] = 0;    // H0     (25)
    //if(limit<= 10 && nQ>= 10)DecayDB[ 10] = 0;    // H-     (37)
    if(limit<= 11 && nQ>= 11)    // Low sigma=pi,pi S-wave : f_0 (800)
    {
      DecayDB[ 11].push_back(new G4QDecayChan(.333,211,-211));
      DecayDB[ 11].push_back(new G4QDecayChan(1.00,111, 111));
    }
    if(limit<= 12 && nQ>= 12)    // Midle Regeon-Pomeron   : f_0 (980)
    {
      DecayDB[ 12].push_back(new G4QDecayChan(.333,211,-211));
      DecayDB[ 12].push_back(new G4QDecayChan(1.00,111, 111));
    }
    if(limit<= 13 && nQ>= 13)    // High Regeon-Pomeron    : f_0 (1500)
    {
      DecayDB[ 13].push_back(new G4QDecayChan(.019,221, 331));
      DecayDB[ 13].push_back(new G4QDecayChan(.070,221, 221));
      DecayDB[ 13].push_back(new G4QDecayChan(.113,311,-311));
      DecayDB[ 13].push_back(new G4QDecayChan(.156,321,-321));
      DecayDB[ 13].push_back(new G4QDecayChan(.578,211,-211)); //@@ include 4pi decays
      DecayDB[ 13].push_back(new G4QDecayChan(1.00,111, 111));
    }
    //if(limit<= 14 && nQ>= 14)DecayDB[ 14].push_back(new G4QDecayChan(1.00,22,22));//Pi0
    //if(limit<= 15 && nQ>= 15)DecayDB[ 15] = 0;    // Pi +
    if(limit<= 16 && nQ>= 16)    // eta
    {
      DecayDB[ 16].push_back(new G4QDecayChan(.226,211,-211,111));
      DecayDB[ 16].push_back(new G4QDecayChan(.551,111, 111,111));
      DecayDB[ 16].push_back(new G4QDecayChan(.598,211,-211, 22));
      DecayDB[ 16].push_back(new G4QDecayChan(.606, 11, -11, 22)); //@@ .002 (pi+)(pi-)2gam
      DecayDB[ 16].push_back(new G4QDecayChan(1.00, 22,  22));
    }
    //if(limit<= 17 && nQ>= 17)    // K 0 (K_short - probab 1/2) @@@@@@@@@@@@
    //{
    //  DecayDB[ 17].push_back(new G4QDecayChan(.6861,211,-211));
    //  DecayDB[ 17].push_back(new G4QDecayChan(1.00, 111, 111));
    //}
    //if(limit<= 18 && nQ>= 18)DecayDB[  8] = 0;    // K +
    if(limit<= 19 && nQ>= 19)    // eta'
    {
      DecayDB[ 19].push_back(new G4QDecayChan(.443,211,-211,221));
      DecayDB[ 19].push_back(new G4QDecayChan(.652,111, 111,221));
      DecayDB[ 19].push_back(new G4QDecayChan(.947, 22, 223));
      DecayDB[ 19].push_back(new G4QDecayChan(.949,111, 111,111));
      DecayDB[ 19].push_back(new G4QDecayChan(.979, 22, 113));
      DecayDB[ 19].push_back(new G4QDecayChan(1.00, 22,  22));
    }
    //if(limit<= 20 && nQ>= 20)DecayDB[ 20] = 0;    // n
    //if(limit<= 21 && nQ>= 21)DecayDB[ 21] = 0;    // p
    //if(limit<= 22 && nQ>= 22)    // Lambda ===>>> all week decays are closed at this time
    //{
    //  DecayDB[ 22].push_back(new G4QDecayChan(.640,2212,-211));
    //  DecayDB[ 22].push_back(new G4QDecayChan(1.00,2112, 111));
    //}
    //if(limit<= 23 &&nQ>=23)DecayDB[23].push_back(new G4QDecayChan(1.,2112,-211));//Sigma-
    if(limit<= 24 &&nQ>=24)DecayDB[24].push_back(new G4QDecayChan(1.,3122,22));//Sigma0(EM)
    //if(limit<= 25 && nQ>= 25)    // Sigma +
    //{
    //  DecayDB[ 25].push_back(new G4QDecayChan(.484,2112, 211));
    //  DecayDB[ 25].push_back(new G4QDecayChan(1.00,2212, 111));
    //}
    //if(limit<= 26 && nQ>=26)DecayDB[26].push_back(new G4QDecayChan(1.,3122,-211));// Ksi-
    //if(limit<= 27 && nQ>=27)DecayDB[27].push_back(new G4QDecayChan(1.,3122, 111));// Ksi0
    if(limit<= 28 && nQ>= 28)DecayDB[ 28].push_back(new G4QDecayChan(1., 211,-211));// rho0
    if(limit<= 29 && nQ>= 29)DecayDB[ 29].push_back(new G4QDecayChan(1., 211, 111));// rho+
    if(limit<= 30 && nQ>= 30)    // omega
    {
      DecayDB[ 30].push_back(new G4QDecayChan(.891, 211,-211,111));
      DecayDB[ 30].push_back(new G4QDecayChan(.908, 211,-211));
      DecayDB[ 30].push_back(new G4QDecayChan(.997,  22, 111));
      DecayDB[ 30].push_back(new G4QDecayChan(.998,  11, -11, 111)); //@@NeedsMoreAccurate
      DecayDB[ 30].push_back(new G4QDecayChan(.998,  13, -13, 111));
      DecayDB[ 30].push_back(new G4QDecayChan(1.00,  22, 221));
    }
    if(limit<= 31 && nQ>= 31)    // K* 0
    {
      DecayDB[ 31].push_back(new G4QDecayChan(.667,-211, 321));
      DecayDB[ 31].push_back(new G4QDecayChan(1.00, 111, 311));
    }
    if(limit<= 32 && nQ>= 32)    // K* +
    {
      DecayDB[ 32].push_back(new G4QDecayChan(.667, 211, 311));
      DecayDB[ 32].push_back(new G4QDecayChan(1.00, 111, 321));
    }
    if(limit<= 33 && nQ>= 33)    // phi
    {
      DecayDB[ 33].push_back(new G4QDecayChan(.491, 311,-311));
      DecayDB[ 33].push_back(new G4QDecayChan(.831, 321,-321));
      DecayDB[ 33].push_back(new G4QDecayChan(.844,  22, 221));
      DecayDB[ 33].push_back(new G4QDecayChan(.846,  22, 111));
      DecayDB[ 33].push_back(new G4QDecayChan(.897, 211,-213));
      DecayDB[ 33].push_back(new G4QDecayChan(.948,-211, 213));
      DecayDB[ 33].push_back(new G4QDecayChan(1.00, 111, 113));
    }
    if(limit<= 34 && nQ>= 34)DecayDB[34].push_back(new G4QDecayChan(1.,2112,-211));//Delta-
    if(limit<= 35 && nQ>= 35)    // Delta 0
    {
      DecayDB[ 35].push_back(new G4QDecayChan(.333,2212,-211));
      DecayDB[ 35].push_back(new G4QDecayChan(1.00,2112, 111));
    }
    if(limit<= 36 && nQ>= 36)    // Delta +
    {
      DecayDB[ 36].push_back(new G4QDecayChan(.333,2112, 211));
      DecayDB[ 36].push_back(new G4QDecayChan(1.00,2212, 111));
    }
    if(limit<= 37 && nQ>= 37)DecayDB[37].push_back(new G4QDecayChan(1.,2212,211));//Delta++
    if(limit<= 38 && nQ>= 38)    // Lambda* (1520)
    {
      DecayDB[ 38].push_back(new G4QDecayChan(.225,3112,-311));
      DecayDB[ 38].push_back(new G4QDecayChan(.450,3222,-321));
      DecayDB[ 38].push_back(new G4QDecayChan(.453,3112,211,111));
      DecayDB[ 38].push_back(new G4QDecayChan(.456,3212,211,-211));
      DecayDB[ 38].push_back(new G4QDecayChan(.459,3212,111,111));
      DecayDB[ 38].push_back(new G4QDecayChan(.462,3222,-211,111));
      DecayDB[ 38].push_back(new G4QDecayChan(.512,3122,211,-211));
      DecayDB[ 38].push_back(new G4QDecayChan(.562,3122,111,111));
      DecayDB[ 38].push_back(new G4QDecayChan(.702,3222,-211));
      DecayDB[ 38].push_back(new G4QDecayChan(.842,3212, 111));
      DecayDB[ 38].push_back(new G4QDecayChan(.982,3112, 211));
      DecayDB[ 38].push_back(new G4QDecayChan(1.00,3122,  22));
    }
    if(limit<= 39 && nQ>= 39)    // Sigma* -
    {
      DecayDB[ 39].push_back(new G4QDecayChan(.060,3112, 111));
      DecayDB[ 39].push_back(new G4QDecayChan(.120,3212,-211));
      DecayDB[ 39].push_back(new G4QDecayChan(1.00,3122,-211));
    }
    if(limit<= 40 && nQ>= 40)    // Sigma* 0
    {
      DecayDB[ 40].push_back(new G4QDecayChan(.040,3112, 211));
      DecayDB[ 40].push_back(new G4QDecayChan(.080,3222,-211));
      DecayDB[ 40].push_back(new G4QDecayChan(.120,3212, 111));
      DecayDB[ 40].push_back(new G4QDecayChan(1.00,3122, 111));
    }
    if(limit<= 41 && nQ>= 41)    // Sigma* +
    {
      DecayDB[ 41].push_back(new G4QDecayChan(.060,3212, 211));
      DecayDB[ 41].push_back(new G4QDecayChan(.120,3222, 111));
      DecayDB[ 41].push_back(new G4QDecayChan(1.00,3122, 211));
    }
    if(limit<= 42 && nQ>= 42)    // Ksi* -
    {
      DecayDB[ 42].push_back(new G4QDecayChan(.667,3322,-211));
      DecayDB[ 42].push_back(new G4QDecayChan(1.00,3312, 111));
    }
    if(limit<= 43 && nQ>= 43)    // Ksi* 0
    {
      DecayDB[ 43].push_back(new G4QDecayChan(.667,3312, 211));
      DecayDB[ 43].push_back(new G4QDecayChan(1.00,3322, 111));
    }
    //if(limit<= 44 && nQ>= 44)    // OMEGA - (Weak)
    //{
    //  DecayDB[ 44].push_back(new G4QDecayChan(.678,3122, 321));
    //  DecayDB[ 44].push_back(new G4QDecayChan(.914,3322,-211));
    //  DecayDB[ 44].push_back(new G4QDecayChan(1.00,3312, 111));
    //}
    if(limit<= 45 && nQ>= 45)    // a_2 0
    {
      DecayDB[ 45].push_back(new G4QDecayChan(.070, 211,-211,223));
      DecayDB[ 45].push_back(new G4QDecayChan(.106, 111, 111,223));
      DecayDB[ 45].push_back(new G4QDecayChan(.131, 321,-321));
      DecayDB[ 45].push_back(new G4QDecayChan(.156, 311,-311));
      DecayDB[ 45].push_back(new G4QDecayChan(.301, 111, 221));
      DecayDB[ 45].push_back(new G4QDecayChan(.534,-211, 213));
      DecayDB[ 45].push_back(new G4QDecayChan(.767, 211,-213));
      DecayDB[ 45].push_back(new G4QDecayChan(1.00, 111, 113));
    }
    if(limit<= 46 && nQ>= 46)    // a_2 +
    {
      DecayDB[ 46].push_back(new G4QDecayChan(.106,111,211,223));
      DecayDB[ 46].push_back(new G4QDecayChan(.156, 321,-311));
      DecayDB[ 46].push_back(new G4QDecayChan(.301, 211, 221));
      DecayDB[ 46].push_back(new G4QDecayChan(.651, 211, 113));
      DecayDB[ 46].push_back(new G4QDecayChan(1.00, 111, 213));
    }
    if(limit<= 47 && nQ>= 47)    // f_2 0
    {
      DecayDB[ 47].push_back(new G4QDecayChan(.005, 221, 221));
      DecayDB[ 47].push_back(new G4QDecayChan(.028, 311,-311));
      DecayDB[ 47].push_back(new G4QDecayChan(.051, 321,-321));
      DecayDB[ 47].push_back(new G4QDecayChan(.123, 111, 113));
      DecayDB[ 47].push_back(new G4QDecayChan(.126, 111, 221));
      DecayDB[ 47].push_back(new G4QDecayChan(.152, 211,-211,113));
      DecayDB[ 47].push_back(new G4QDecayChan(.717, 211,-211));
      DecayDB[ 47].push_back(new G4QDecayChan(1.00, 111, 111));
    }
    if(limit<= 48 && nQ>= 48)    // K_2 0
    {
      DecayDB[ 48].push_back(new G4QDecayChan(.028, 311, 223));
      DecayDB[ 48].push_back(new G4QDecayChan(.074, 211,-211,313));
      DecayDB[ 48].push_back(new G4QDecayChan(.143,111,-211,323));
      DecayDB[ 48].push_back(new G4QDecayChan(.166,111, 111,313));
      DecayDB[ 48].push_back(new G4QDecayChan(.190,-211, 323));
      DecayDB[ 48].push_back(new G4QDecayChan(.314, 111, 313));
      DecayDB[ 48].push_back(new G4QDecayChan(.357, 311, 113));
      DecayDB[ 48].push_back(new G4QDecayChan(.500, 321,-213));
      DecayDB[ 48].push_back(new G4QDecayChan(.750,-211, 321));
      DecayDB[ 48].push_back(new G4QDecayChan(1.00, 111, 311));
    }
    if(limit<= 49 && nQ>= 49)    // K_2 +
    {
      DecayDB[ 49].push_back(new G4QDecayChan(.028, 321, 223));
      DecayDB[ 49].push_back(new G4QDecayChan(.074,211,-211,323));
      DecayDB[ 49].push_back(new G4QDecayChan(.143,111, 211,313));
      DecayDB[ 49].push_back(new G4QDecayChan(.166,111, 111,323));
      DecayDB[ 49].push_back(new G4QDecayChan(.190, 211, 313));
      DecayDB[ 49].push_back(new G4QDecayChan(.314, 111, 323));
      DecayDB[ 49].push_back(new G4QDecayChan(.357, 311, 213));
      DecayDB[ 49].push_back(new G4QDecayChan(.500, 321, 113));
      DecayDB[ 49].push_back(new G4QDecayChan(.750, 211, 311));
      DecayDB[ 49].push_back(new G4QDecayChan(1.00, 111, 321));
    }
    if(limit<= 50 && nQ>= 50)    // f_2' 0
    {
      DecayDB[ 50].push_back(new G4QDecayChan(.103, 221, 221));
      DecayDB[ 50].push_back(new G4QDecayChan(.547, 311,-311));
      DecayDB[ 50].push_back(new G4QDecayChan(.991, 321,-321));
      DecayDB[ 50].push_back(new G4QDecayChan(.997, 211,-211));
      DecayDB[ 50].push_back(new G4QDecayChan(1.00, 111, 111));
    }
    if(limit<= 51 && nQ>= 51)    // N_5/2 0
    {
      DecayDB[ 51].push_back(new G4QDecayChan(.040, 211, 1114));
      DecayDB[ 51].push_back(new G4QDecayChan(.080, 111, 2114));
      DecayDB[ 51].push_back(new G4QDecayChan(.120,-211, 2214));
      DecayDB[ 51].push_back(new G4QDecayChan(.180, 2112, 113));
      DecayDB[ 51].push_back(new G4QDecayChan(.210, 2212,-213));
      DecayDB[ 51].push_back(new G4QDecayChan(.340, 2112, 110));
      DecayDB[ 51].push_back(new G4QDecayChan(.780, 2212,-211));
      DecayDB[ 51].push_back(new G4QDecayChan(1.00, 2112, 111));
    }
    if(limit<= 52 && nQ>= 52)    // N_5/2 +
    {
      DecayDB[ 52].push_back(new G4QDecayChan(.040,-211, 2224));
      DecayDB[ 52].push_back(new G4QDecayChan(.080, 211, 2114));
      DecayDB[ 52].push_back(new G4QDecayChan(.120, 111, 2214));
      DecayDB[ 52].push_back(new G4QDecayChan(.180, 2112, 213));
      DecayDB[ 52].push_back(new G4QDecayChan(.210, 2212, 113));
      DecayDB[ 52].push_back(new G4QDecayChan(.340, 2212, 229));
      DecayDB[ 52].push_back(new G4QDecayChan(.780, 2112, 211));
      DecayDB[ 52].push_back(new G4QDecayChan(1.00, 2212, 111));
    }
    if(limit<= 53 && nQ>= 53)    // LAMBDA_5/2
    {
      DecayDB[ 53].push_back(new G4QDecayChan(.350, 2112,-311));
      DecayDB[ 53].push_back(new G4QDecayChan(.700, 2212,-321));
      DecayDB[ 53].push_back(new G4QDecayChan(.740, 211, 3114));
      DecayDB[ 53].push_back(new G4QDecayChan(.780,-211, 3224));
      DecayDB[ 53].push_back(new G4QDecayChan(.820, 111, 3214));
      DecayDB[ 53].push_back(new G4QDecayChan(.880, 3112, 211));
      DecayDB[ 53].push_back(new G4QDecayChan(.940, 3222,-211));
      DecayDB[ 53].push_back(new G4QDecayChan(1.00, 3212, 111));
    }
    if(limit<= 54 && nQ>= 54)    // SIGMA_5/2 -
    {
      DecayDB[ 54].push_back(new G4QDecayChan(.600, 2112,-321));
      DecayDB[ 54].push_back(new G4QDecayChan(.660,-211, 3214));
      DecayDB[ 54].push_back(new G4QDecayChan(.720, 111, 3114));
      DecayDB[ 54].push_back(new G4QDecayChan(.810, 3212,-211));
      DecayDB[ 54].push_back(new G4QDecayChan(.900, 3112, 111));
      DecayDB[ 54].push_back(new G4QDecayChan(1.00, 3122,-211));
    }
    if(limit<= 55 && nQ>= 55)    // SIGMA_5/2 0
    {
      DecayDB[ 55].push_back(new G4QDecayChan(.300, 2112,-311));
      DecayDB[ 55].push_back(new G4QDecayChan(.600, 2212,-321));
      DecayDB[ 55].push_back(new G4QDecayChan(.640, 211, 3114));
      DecayDB[ 55].push_back(new G4QDecayChan(.680,-211, 3224));
      DecayDB[ 55].push_back(new G4QDecayChan(.720, 111, 3214));
      DecayDB[ 55].push_back(new G4QDecayChan(.780, 3112, 211));
      DecayDB[ 55].push_back(new G4QDecayChan(.840, 3222,-211));
      DecayDB[ 55].push_back(new G4QDecayChan(.900, 3212, 111));
      DecayDB[ 55].push_back(new G4QDecayChan(1.00, 3122, 111));
    }
    if(limit<= 56 && nQ>= 56)    // SIGMA_5/2 +
    {
      DecayDB[ 56].push_back(new G4QDecayChan(.600, 2212,-311));
      DecayDB[ 56].push_back(new G4QDecayChan(.660, 211, 3214));
      DecayDB[ 56].push_back(new G4QDecayChan(.720, 111, 3224));
      DecayDB[ 56].push_back(new G4QDecayChan(.810, 3212, 211));
      DecayDB[ 56].push_back(new G4QDecayChan(.900, 3222, 111));
      DecayDB[ 56].push_back(new G4QDecayChan(1.00, 3122, 211));
    }
    if(limit<= 57 && nQ>= 57)    // KSI_5/2 -
    {
      DecayDB[ 57].push_back(new G4QDecayChan(.400, 3112,-311));
      DecayDB[ 57].push_back(new G4QDecayChan(.800, 3212,-321));
      DecayDB[ 57].push_back(new G4QDecayChan(1.00, 3122,-321));
    }
    if(limit<= 58 && nQ>= 58)    // KSI_5/2 0
    {
      DecayDB[ 58].push_back(new G4QDecayChan(.400, 3212,-311));
      DecayDB[ 58].push_back(new G4QDecayChan(.800, 3222,-321));
      DecayDB[ 58].push_back(new G4QDecayChan(1.00, 3122,-311));
    }
    if(limit<= 59 && nQ>= 59)    // rho_3 0
    {
      DecayDB[ 59].push_back(new G4QDecayChan(.019,311,-313));
      DecayDB[ 59].push_back(new G4QDecayChan(.038,321,-323));
      DecayDB[ 59].push_back(new G4QDecayChan(.046,311,-311));
      DecayDB[ 59].push_back(new G4QDecayChan(.054,321,-321));
      DecayDB[ 59].push_back(new G4QDecayChan(.224,111, 223));
      DecayDB[ 59].push_back(new G4QDecayChan(.404,111,-211,213));
      DecayDB[ 59].push_back(new G4QDecayChan(.584,111, 211,-213));
      DecayDB[ 59].push_back(new G4QDecayChan(.764,111, 111,113));
      DecayDB[ 59].push_back(new G4QDecayChan(1.00,211,-211));
    }
    if(limit<= 60 && nQ>= 60)    // rho_3 +
    {
      DecayDB[ 60].push_back(new G4QDecayChan(.019, 321,-313));
      DecayDB[ 60].push_back(new G4QDecayChan(.038,-311, 323));
      DecayDB[ 60].push_back(new G4QDecayChan(.054, 321,-311));
      DecayDB[ 60].push_back(new G4QDecayChan(.224, 211, 223));
      DecayDB[ 60].push_back(new G4QDecayChan(.404,211,-211,213));
      DecayDB[ 60].push_back(new G4QDecayChan(.584,211,211,-213));
      DecayDB[ 60].push_back(new G4QDecayChan(.764,211,111,113));
      DecayDB[ 60].push_back(new G4QDecayChan(1.00, 211, 111));
    }
    if(limit<= 61 && nQ>= 61)    // omega_3
    {
      DecayDB[ 61].push_back(new G4QDecayChan(.020,211,-211,223));
      DecayDB[ 61].push_back(new G4QDecayChan(.040,111, 111,223));
      DecayDB[ 61].push_back(new G4QDecayChan(.060, 211,-213));
      DecayDB[ 61].push_back(new G4QDecayChan(.080,-211, 213));
      DecayDB[ 61].push_back(new G4QDecayChan(1.00, 111, 113));
    }
    if(limit<= 62 && nQ>= 62)    // K_3 0
    {
      DecayDB[ 62].push_back(new G4QDecayChan(.030, 111, 315));
      DecayDB[ 62].push_back(new G4QDecayChan(.060,-211, 325));
      DecayDB[ 62].push_back(new G4QDecayChan(.340, 311, 331));
      DecayDB[ 62].push_back(new G4QDecayChan(.400, 111, 313));
      DecayDB[ 62].push_back(new G4QDecayChan(.520,-211, 323));
      DecayDB[ 62].push_back(new G4QDecayChan(.620, 311, 113));
      DecayDB[ 62].push_back(new G4QDecayChan(.820, 321,-213));
      DecayDB[ 62].push_back(new G4QDecayChan(.940,-211, 321));
      DecayDB[ 62].push_back(new G4QDecayChan(1.00, 111, 311));
    }
    if(limit<= 63 && nQ>= 63)    // K_3 +
    {
      DecayDB[ 63].push_back(new G4QDecayChan(.030, 211, 315));
      DecayDB[ 63].push_back(new G4QDecayChan(.060, 111, 325));
      DecayDB[ 63].push_back(new G4QDecayChan(.340, 321, 331));
      DecayDB[ 63].push_back(new G4QDecayChan(.400, 211, 313));
      DecayDB[ 63].push_back(new G4QDecayChan(.520, 111, 323));
      DecayDB[ 63].push_back(new G4QDecayChan(.620, 311, 213));
      DecayDB[ 63].push_back(new G4QDecayChan(.820, 321, 113));
      DecayDB[ 63].push_back(new G4QDecayChan(.940, 211, 311));
      DecayDB[ 63].push_back(new G4QDecayChan(1.00, 111, 321));
    }
    if(limit<= 64 && nQ>= 64)    // phi_3
    {
      DecayDB[ 64].push_back(new G4QDecayChan(.250, 321,-321));
      DecayDB[ 64].push_back(new G4QDecayChan(.500, 311,-311));
      DecayDB[ 64].push_back(new G4QDecayChan(.625, 321,-323));
      DecayDB[ 64].push_back(new G4QDecayChan(.750,-321, 323));
      DecayDB[ 64].push_back(new G4QDecayChan(.875, 311,-313));
      DecayDB[ 64].push_back(new G4QDecayChan(1.00,-311, 313));
    }
    if(limit<= 65 && nQ>= 65)    // DELTA_7/2 -
    {
      DecayDB[ 65].push_back(new G4QDecayChan(.200, 2112,-213 ));
      DecayDB[ 65].push_back(new G4QDecayChan(.320,-211 , 2114));
      DecayDB[ 65].push_back(new G4QDecayChan(.500, 111 , 1114));
      DecayDB[ 65].push_back(new G4QDecayChan(1.00, 2112,-211 ));
    }
    if(limit<= 66 && nQ>= 66)    // DELTA_7/2 0
    {
      DecayDB[ 66].push_back(new G4QDecayChan(.133, 2112, 113 ));
      DecayDB[ 66].push_back(new G4QDecayChan(.200, 2212,-213 ));
      DecayDB[ 66].push_back(new G4QDecayChan(.360,-211 , 2214));
      DecayDB[ 66].push_back(new G4QDecayChan(.480, 211 , 1114));
      DecayDB[ 66].push_back(new G4QDecayChan(.500, 111 , 2114));
      DecayDB[ 66].push_back(new G4QDecayChan(.666, 2212,-211 ));
      DecayDB[ 66].push_back(new G4QDecayChan(1.00, 2112, 111 ));
    }
    if(limit<= 67 && nQ>= 67)    // DELTA_7/2 +
    {
      DecayDB[ 67].push_back(new G4QDecayChan(.133, 2112, 213 ));
      DecayDB[ 67].push_back(new G4QDecayChan(.200, 2212, 113 ));
      DecayDB[ 67].push_back(new G4QDecayChan(.360,-211 , 2224));
      DecayDB[ 67].push_back(new G4QDecayChan(.480, 211 , 2114));
      DecayDB[ 67].push_back(new G4QDecayChan(.500, 111 , 2214));
      DecayDB[ 67].push_back(new G4QDecayChan(.666, 2112, 211 ));
      DecayDB[ 67].push_back(new G4QDecayChan(1.00, 2212, 111 ));
    }
    if(limit<= 68 && nQ>= 68)    // DELTA_7/2 ++
    {
      DecayDB[ 68].push_back(new G4QDecayChan(.200, 2212, 213 ));
      DecayDB[ 68].push_back(new G4QDecayChan(.320, 211 , 2214));
      DecayDB[ 68].push_back(new G4QDecayChan(.500, 111 , 2224));
      DecayDB[ 68].push_back(new G4QDecayChan(1.00, 2212, 211 ));
    }
    if(limit<= 69 && nQ>= 69)     // LAMBDA_7/2
    {
      DecayDB[ 69].push_back(new G4QDecayChan(.160, 3122, 223 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.260, 2112,-313 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.360, 2212,-323 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.400, 3312, 321 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.440, 3322, 311 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.480, 3122, 221 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.520, 2112,-311 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.560, 2212,-321 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.600, 3112, 211 ));
      DecayDB[ 69].push_back(new G4QDecayChan(.800, 3222,-211 ));
      DecayDB[ 69].push_back(new G4QDecayChan(1.00, 3212, 111 ));
    }
    if(limit<= 70 && nQ>= 70)    // SIGMA_7/2 -
    {
      DecayDB[ 70].push_back(new G4QDecayChan(.030, 2112,-323 ));
      DecayDB[ 70].push_back(new G4QDecayChan(.165,-311 , 1114));
      DecayDB[ 70].push_back(new G4QDecayChan(.210,-321 , 2114));
      DecayDB[ 70].push_back(new G4QDecayChan(.390,-211 , 3124));
      DecayDB[ 70].push_back(new G4QDecayChan(.450,-211 , 3214));
      DecayDB[ 70].push_back(new G4QDecayChan(.510, 111 , 3114));
      DecayDB[ 70].push_back(new G4QDecayChan(.540, 311 , 3314));
      DecayDB[ 70].push_back(new G4QDecayChan(.570,-211 , 3212));
      DecayDB[ 70].push_back(new G4QDecayChan(.600, 111 , 3112));
      DecayDB[ 70].push_back(new G4QDecayChan(.780,-321 , 2112));
      DecayDB[ 70].push_back(new G4QDecayChan(1.00,-211 , 3122));
    }
    if(limit<= 71 && nQ>= 71)    // SIGMA_7/2 0
    {
      DecayDB[ 71].push_back(new G4QDecayChan(.015, 2112,-313 ));
      DecayDB[ 71].push_back(new G4QDecayChan(.030, 2212,-321 ));
      DecayDB[ 71].push_back(new G4QDecayChan(.120,-311 , 2114));
      DecayDB[ 71].push_back(new G4QDecayChan(.210,-321 , 2214));
      DecayDB[ 71].push_back(new G4QDecayChan(.390, 111 , 3124));
      DecayDB[ 71].push_back(new G4QDecayChan(.450,-211 , 3224));
      DecayDB[ 71].push_back(new G4QDecayChan(.510, 211 , 3114));
      DecayDB[ 71].push_back(new G4QDecayChan(.525, 311 , 3324));
      DecayDB[ 71].push_back(new G4QDecayChan(.540, 321 , 3314));
      DecayDB[ 71].push_back(new G4QDecayChan(.570,-211 , 3222));
      DecayDB[ 71].push_back(new G4QDecayChan(.600, 211 , 3112));
      DecayDB[ 71].push_back(new G4QDecayChan(.690,-311 , 2112));
      DecayDB[ 71].push_back(new G4QDecayChan(.780,-321 , 2212));
      DecayDB[ 71].push_back(new G4QDecayChan(1.00, 111 , 3122));
    }
    if(limit<= 72 && nQ>= 72)    // SIGMA_7/2 +
    {
      DecayDB[ 72].push_back(new G4QDecayChan(.030, 2212,-313 ));
      DecayDB[ 72].push_back(new G4QDecayChan(.165,-321 , 2224));
      DecayDB[ 72].push_back(new G4QDecayChan(.210,-311 , 2214));
      DecayDB[ 72].push_back(new G4QDecayChan(.390, 211 , 3124));
      DecayDB[ 72].push_back(new G4QDecayChan(.450, 211 , 3214));
      DecayDB[ 72].push_back(new G4QDecayChan(.510, 111 , 3224));
      DecayDB[ 72].push_back(new G4QDecayChan(.540, 321 , 3324));
      DecayDB[ 72].push_back(new G4QDecayChan(.570, 211 , 3212));
      DecayDB[ 72].push_back(new G4QDecayChan(.600, 111 , 3222));
      DecayDB[ 72].push_back(new G4QDecayChan(.780,-311 , 2212));
      DecayDB[ 72].push_back(new G4QDecayChan(1.00, 211 , 3122));
    }
    if(limit<= 73 && nQ>= 73)    // KSI_7/2 -
    {
      DecayDB[ 73].push_back(new G4QDecayChan(.400, 3112,-311));
      DecayDB[ 73].push_back(new G4QDecayChan(.800, 3212,-321));
      DecayDB[ 73].push_back(new G4QDecayChan(1.00, 3122,-321));
    }
    if(limit<= 74 && nQ>= 74)    // KSI_7/2 0
    {
      DecayDB[ 74].push_back(new G4QDecayChan(.400, 3212,-311));
      DecayDB[ 74].push_back(new G4QDecayChan(.800, 3222,-321));
      DecayDB[ 74].push_back(new G4QDecayChan(1.00, 3122,-311));
    }
    if(limit<= 75 && nQ>= 75)    // OMEGA_7/2 -
    {
      DecayDB[ 75].push_back(new G4QDecayChan(.250,-311 , 3314));
      DecayDB[ 75].push_back(new G4QDecayChan(.500,-321 , 3324));
      DecayDB[ 75].push_back(new G4QDecayChan(.750, 3312,-313 ));
      DecayDB[ 75].push_back(new G4QDecayChan(1.00, 3322,-323 ));
    }
    if(limit<= 76 && nQ>= 76)    // a_4 0
    {
      DecayDB[ 76].push_back(new G4QDecayChan(.200, 311,-311));
      DecayDB[ 76].push_back(new G4QDecayChan(.400, 321,-321));
      DecayDB[ 76].push_back(new G4QDecayChan(.600, 111, 221));
      DecayDB[ 76].push_back(new G4QDecayChan(1.00, 111, 211,-211));
    }
    if(limit<= 77 && nQ>= 77)    // a_4 +
    {
      DecayDB[ 77].push_back(new G4QDecayChan(.400, 321,-311));
      DecayDB[ 77].push_back(new G4QDecayChan(.600, 211, 221));
      DecayDB[ 77].push_back(new G4QDecayChan(.800, 211, 211,-211));
      DecayDB[ 77].push_back(new G4QDecayChan(1.00, 211, 111, 111));
    }
    if(limit<= 78 && nQ>= 78)    // f_4 0
    {
      DecayDB[ 78].push_back(new G4QDecayChan(.020, 333, 333));
      DecayDB[ 78].push_back(new G4QDecayChan(.340, 223, 223));
      DecayDB[ 78].push_back(new G4QDecayChan(.350, 221, 221));
      DecayDB[ 78].push_back(new G4QDecayChan(.360, 311,-311));
      DecayDB[ 78].push_back(new G4QDecayChan(.370, 321,-321));
      DecayDB[ 78].push_back(new G4QDecayChan(.670, 213,-213));
      DecayDB[ 78].push_back(new G4QDecayChan(.820, 113, 113));
      DecayDB[ 78].push_back(new G4QDecayChan(.940, 211,-211));
      DecayDB[ 78].push_back(new G4QDecayChan(1.00, 111, 111));
    }
    if(limit<= 79 && nQ>= 79)    // K_4 0
    {
      DecayDB[ 79].push_back(new G4QDecayChan(.060, 333, 313));
      DecayDB[ 79].push_back(new G4QDecayChan(.260, 223, 313));
      DecayDB[ 79].push_back(new G4QDecayChan(.380, 313, 113));
      DecayDB[ 79].push_back(new G4QDecayChan(.400, 323,-213));
      DecayDB[ 79].push_back(new G4QDecayChan(.500, 223, 311));
      DecayDB[ 79].push_back(new G4QDecayChan(.550, 311, 113));
      DecayDB[ 79].push_back(new G4QDecayChan(.600, 321,-213));
      DecayDB[ 79].push_back(new G4QDecayChan(.700, 311, 113));
      DecayDB[ 79].push_back(new G4QDecayChan(.800, 321,-213));
      DecayDB[ 79].push_back(new G4QDecayChan(.900, 311, 111));
      DecayDB[ 79].push_back(new G4QDecayChan(1.00, 321,-211));
    }
    if(limit<= 80 && nQ>= 80)    // K_4 +
    {
      DecayDB[ 80].push_back(new G4QDecayChan(.060, 333, 323));
      DecayDB[ 80].push_back(new G4QDecayChan(.260, 223, 323));
      DecayDB[ 80].push_back(new G4QDecayChan(.380, 313, 213));
      DecayDB[ 80].push_back(new G4QDecayChan(.400, 323, 113));
      DecayDB[ 80].push_back(new G4QDecayChan(.500, 223, 321));
      DecayDB[ 80].push_back(new G4QDecayChan(.550, 321, 113));
      DecayDB[ 80].push_back(new G4QDecayChan(.600, 311, 213));
      DecayDB[ 80].push_back(new G4QDecayChan(.700, 311, 211));
      DecayDB[ 80].push_back(new G4QDecayChan(.800, 321, 111));
      DecayDB[ 80].push_back(new G4QDecayChan(.900, 311, 211));
      DecayDB[ 80].push_back(new G4QDecayChan(1.00, 321, 111));
    }
    if(limit<=81&&nQ>=81)DecayDB[81].push_back(new G4QDecayChan(1., 333,333));//phi_4(2300)
    if(limit<=82&&nQ>=82)DecayDB[82].push_back(new G4QDecayChan(1.,2212, 2224));//pDelta++
    if(limit<=83&&nQ>=83)DecayDB[83].push_back(new G4QDecayChan(1.,2112, 1114));//nDelta-
    if(limit<=84&&nQ>=84)DecayDB[84].push_back(new G4QDecayChan(1.,2224, 2224));//D++D++
    if(limit<=85&&nQ>=85)DecayDB[85].push_back(new G4QDecayChan(1.,1114, 1114));//Del-Del-
    if(limit<=86&&nQ>=86)DecayDB[86].push_back(new G4QDecayChan(1.,2212,2212,2224));//ppD++
    if(limit<=87&&nQ>=87)DecayDB[87].push_back(new G4QDecayChan(1.,2112,2112,1114));//nnD-
    if(limit<=88&&nQ>=88)DecayDB[88].push_back(new G4QDecayChan(1.,2212,2224,2224));//p2D++
    if(limit<=89&&nQ>=89)DecayDB[89].push_back(new G4QDecayChan(1.,2112,1114,1114));//nD-D-
    //if(limit<=90&&nQ>=90)DecayDB[90] = 0; // neutron (n) as a quark exchange fragment
    //if(limit<=91&&nQ>=91)DecayDB[91] = 0; // proton (p)  as a quark exchange fragment
    //if(limit<=92&&nQ>=92)DecayDB[92] = 0; // neutron (L/Sigma0) as aQuarkExchangeFragment
    //if(limit<=93&&nQ>=93)DecayDB[93] = 0; // neutron (Sigma-) as a quarkExchangeFragment
    //if(limit<=94&&nQ>=94)DecayDB[94] = 0; // neutron (Sigma+) as a quarkExchangeFragment
    //if(limit<=95&&nQ>=95)DecayDB[95] = 0; // neutron (Xi-) as a quark exchange fragment
    //if(limit<=96&&nQ>=96)DecayDB[96] = 0; // neutron (Xi0) as a quark exchange fragment
    //if(limit<=97&&nQ>=97)DecayDB[97] = 0; // neutron (Omega-) as a quarkExchangeFragment
    if(limit<=98&&nQ>=98)DecayDB[98].push_back(new G4QDecayChan(1.,2112, 2112)); //nn (src)
    if(limit<=99&&nQ>=99)DecayDB[99].push_back(new G4QDecayChan(1.,2212, 2112)); //d/pn(?)
    if(limit<=100&&nQ>=100)DecayDB[100].push_back(new G4QDecayChan(1.,2212,2212));//pp(src)
    if(limit<=101&&nQ>=101)DecayDB[101].push_back(new G4QDecayChan(1.,3122,2112));//Ln
    if(limit<=102&&nQ>=102)DecayDB[102].push_back(new G4QDecayChan(1.,3122,2212));//Lp
    if(limit<=104&&nQ>=104)DecayDB[104].push_back(new G4QDecayChan(1.,3112,2112));//nSig-
    if(limit<=103&&nQ>=103)DecayDB[103].push_back(new G4QDecayChan(1.,3122,3122));//LL
    if(limit<=105&&nQ>=105)DecayDB[105].push_back(new G4QDecayChan(1.,3222,2212));//pSig+
    //if(limit<=106&&nQ>=106)DecayDB[106] = 0; // t
    //if(limit<=107&&nQ>=107)DecayDB[107] = 0; // He3
    //Lnn=>Lambda+n+n decay to avoid the final state Hypernucleus
    if(limit<=108&&nQ>=108)DecayDB[108].push_back(new G4QDecayChan(1.,3122,2112,2112));
    if(limit<=109&&nQ>=109)DecayDB[109].push_back(new G4QDecayChan(1.,3122,90001001));// Ld
    //Lpp=>Lambda+p+p decay to avoid the final state Hypernucleus
    if(limit<=110&&nQ>=110)DecayDB[110].push_back(new G4QDecayChan(1.,3122,2212,2212));
    //LLn=>Lambda+Lambda+n decay to avoid the final state Hypernucleus
    if(limit<=111&&nQ>=111)DecayDB[111].push_back(new G4QDecayChan(1.,3122,3122,2112));
    //LLp=>Lambda+Lambda+p decay to avoid the final state Hypernucleus
    if(limit<=112&&nQ>=112)DecayDB[112].push_back(new G4QDecayChan(1.,3122,3122,2212));
    // nnSigma-=>n+n+Sigma- decay to avoid the final state Hypernucleus
    if(limit<=113&&nQ>=113)DecayDB[113].push_back(new G4QDecayChan(1.,2112,2112,3112));
    // ------- Nuclear fragments
    //if(limit<=114 && nQ>=114)
    //{
    //  if(limit<114) limit=101;
    //  for (int i=limit; i<nQ; i++) DecayDB[i] = 0;
    //}
    //Update the limit
    limit=nQ+1;
#ifdef debug
    G4cout<<"G4QParticle::InitDecayVector: limit is set to "<<limit<<G4endl;
#endif
  }
  //if(!nQ)G4cout<<"G4QParticle::InitDecayVector:Q=0,nD="<<DecayDB[abs(nQ)].size()<<G4endl;
  return DecayDB[std::abs(nQ)];
}

// Initialize the Particle by a Q Code
void G4QParticle::InitQParticle(G4int theQCode)
//    =============================================
{
  aQPDG.InitByQCode(theQCode);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(theQCode);
  //if(!theQCode)G4cout<<"G4QPar::InitQP:PDG="<<GetPDGCode()<<",n="<<aDecay.size()<<G4endl;
}

// Initialize the Particle by a PDG Code
void G4QParticle::InitPDGParticle(G4int thePDGCode)
//    =============================================
{
  aQPDG      = G4QPDGCode(thePDGCode);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(aQPDG.GetQCode());
}
