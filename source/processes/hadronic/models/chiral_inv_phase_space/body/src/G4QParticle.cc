// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QParticle.cc,v 1.1 2000-08-17 13:55:49 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QParticle ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Particles in the CHIPS Model
// -------------------------------------------------------------------
 
//#define debug
//#define pdebug

#include "G4QParticleVector.hh"

G4QParticle::G4QParticle() {};

G4QParticle::G4QParticle(G4int thePDG)
{
  aQPDG      = G4QPDGCode(thePDG);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(aQPDG.GetQCode());
}

G4QParticle::~G4QParticle() {aDecay.clearAndDestroy();};

// Standard output for QParticle
ostream& operator<<(ostream& lhs, G4QParticle& rhs)
//       =========================================
{
  G4QPDGCode rhsQPDG = rhs.GetQPDG();
  lhs << endl << "Particle with PDG=" << rhsQPDG << ", Spin=" << rhs.GetSpin()
      << ", mass=" << rhs.GetMass() << ", width=" << rhs.GetWidth() << endl;
  lhs << " Quark Content of the Particle=" << rhs.GetQContent() << ", Decay Channels:" << endl;
  G4QDecayChanVector DCV = rhs.GetDecayVector();
  G4int n = DCV.entries();
  for (int i=0; i<n; i++)
  {
    lhs << DCV[i]->GetDecayChanLimit() << "PDG codes";
    G4QPDGCodeVector PCV=DCV[i]->GetVecOfSecHadrons();
    G4int m = PCV.entries();
    for (int j=0; j<n; j++)
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
  static G4int nP = 486;                  // Up to A=80
  static G4QDecayChanVector* DecayDB = new G4QDecayChanVector[nP];
  static int limit= 0;
  if(nQ>=limit && nQ<nP)
  {
    // *** Secondary PDG-particles should be ordered in a Channel by increasing width ***!!!***
    // *** Channels should be ordered by increasing minimum mass of the secondary particles ***
    //if(limit<=  0 && nQ>=  0)DecayDB[  0] = NULL;    // gamma
    if(limit<=  1 && nQ>=  1)    // Low sigma=pi,pi S-wave : f_0 (800)
	{
      DecayDB[  1].insert(new G4QDecayChan(.333,211,-211));
      DecayDB[  1].insert(new G4QDecayChan(1.00,111, 111));
	}
    if(limit<=  2 && nQ>=  2)    // Midle Regeon-Pomeron   : f_0 (980)
	{
      DecayDB[  2].insert(new G4QDecayChan(.333,211,-211));
      DecayDB[  2].insert(new G4QDecayChan(1.00,111, 111));
	}
    if(limit<=  3 && nQ>=  3)    // High Regeon-Pomeron    : f_0 (1500)
	{
      DecayDB[  3].insert(new G4QDecayChan(.300,221, 331));
      DecayDB[  3].insert(new G4QDecayChan(.600,221, 221));
      DecayDB[  3].insert(new G4QDecayChan(.650,311,-311));
      DecayDB[  3].insert(new G4QDecayChan(.700,321,-321));
      DecayDB[  3].insert(new G4QDecayChan(.800,211,-211));
      DecayDB[  3].insert(new G4QDecayChan(1.00,111, 111));
	}
    //if(limit<=  4 && nQ>=  4)DecayDB[  4].insert(new G4QDecayChan(1.00, 22,  22));//Pi 0 @@!@@
    //if(limit<=  5 && nQ>=  5)DecayDB[  5] = NULL;    // Pi +
    if(limit<=  6 && nQ>=  6)    // eta
	{
      DecayDB[  6].insert(new G4QDecayChan(.231,211,-211,111));
      DecayDB[  6].insert(new G4QDecayChan(.553,111, 111,111));
      DecayDB[  6].insert(new G4QDecayChan(.603,211,-211, 22));
      DecayDB[  6].insert(new G4QDecayChan(1.00, 22,  22));
	}
    //if(limit<=  7 && nQ>=  7)    // K 0 (K_short - probab 1/2)
	//{
    //  DecayDB[  7].insert(new G4QDecayChan(.6861,211,-211));
    //  DecayDB[  7].insert(new G4QDecayChan(1.00, 111, 111));
	//}
    //if(limit<=  8 && nQ>=  8)DecayDB[  8] = NULL;    // K +
    if(limit<=  9 && nQ>=  9)    // eta'
	{
      DecayDB[  9].insert(new G4QDecayChan(.438,211,-211,221));
      DecayDB[  9].insert(new G4QDecayChan(.645,111, 111,221));
      DecayDB[  9].insert(new G4QDecayChan(.675, 22, 223));
      DecayDB[  9].insert(new G4QDecayChan(.677,111, 111,111));
      DecayDB[  9].insert(new G4QDecayChan(.979, 22, 113));
      DecayDB[  9].insert(new G4QDecayChan(1.00, 22,  22));
	}
    //if(limit<= 10 && nQ>= 10)DecayDB[ 10] = NULL;    // n
    //if(limit<= 11 && nQ>= 11)DecayDB[ 11] = NULL;    // p
    //if(limit<= 12 && nQ>= 12)    // Lambda ===>>> all week decays are closed at this time
	//{
    //  DecayDB[ 12].insert(new G4QDecayChan(.640,2212,-211));
    //  DecayDB[ 12].insert(new G4QDecayChan(1.00,2112, 111));
	//}
    //if(limit<= 13 && nQ>= 13)DecayDB[ 13].insert(new G4QDecayChan(1.00,2112,-211)); // Sigma -
    if(limit<= 14 && nQ>= 14)DecayDB[ 14].insert(new G4QDecayChan(1.00,3122,  22));//Sigma 0(EM)
    //if(limit<= 15 && nQ>= 15)    // Sigma +
	//{
    //  DecayDB[ 15].insert(new G4QDecayChan(.484,2112, 211));
    //  DecayDB[ 15].insert(new G4QDecayChan(1.00,2212, 111));
	//}
    //if(limit<= 16 && nQ>= 16)DecayDB[ 16].insert(new G4QDecayChan(1.00,3122,-211));   // Ksi -
    //if(limit<= 17 && nQ>= 17)DecayDB[ 17].insert(new G4QDecayChan(1.00,3122, 111));   // Ksi 0
    if(limit<= 18 && nQ>= 18)DecayDB[ 18].insert(new G4QDecayChan(1.00, 211,-211));   // rho 0
    if(limit<= 19 && nQ>= 19)DecayDB[ 19].insert(new G4QDecayChan(1.00, 211, 111));   // rho +
    if(limit<= 20 && nQ>= 20)    // omega
	{
      DecayDB[ 20].insert(new G4QDecayChan(.892, 211,-211,111));
      DecayDB[ 20].insert(new G4QDecayChan(.914, 211,-211));
      DecayDB[ 20].insert(new G4QDecayChan(1.00,  22, 111));
	}
    if(limit<= 21 && nQ>= 21)    // K* 0
	{
      DecayDB[ 21].insert(new G4QDecayChan(.667,-211, 321));
      DecayDB[ 21].insert(new G4QDecayChan(1.00, 111, 311));
	}
    if(limit<= 22 && nQ>= 22)    // K* +
	{
      DecayDB[ 22].insert(new G4QDecayChan(.667, 211, 311));
      DecayDB[ 22].insert(new G4QDecayChan(1.00, 111, 321));
	}
    if(limit<= 23 && nQ>= 23)    // phi
	{
      DecayDB[ 23].insert(new G4QDecayChan(.341, 311,-311));
      DecayDB[ 23].insert(new G4QDecayChan(.832, 321,-321));
      DecayDB[ 23].insert(new G4QDecayChan(.850,  22, 221));
      DecayDB[ 23].insert(new G4QDecayChan(.900, 211,-213));
      DecayDB[ 23].insert(new G4QDecayChan(.950,-211, 213));
      DecayDB[ 23].insert(new G4QDecayChan(1.00, 111, 113));
	}
    if(limit<= 24 && nQ>= 24)DecayDB[ 24].insert(new G4QDecayChan(1.00,2112,-211)); // Delta -
    if(limit<= 25 && nQ>= 25)    // Delta 0
	{
      DecayDB[ 25].insert(new G4QDecayChan(.333,2212,-211));
      DecayDB[ 25].insert(new G4QDecayChan(1.00,2112, 111));
	}
    if(limit<= 26 && nQ>= 26)    // Delta +
	{
      DecayDB[ 26].insert(new G4QDecayChan(.333,2112, 211));
      DecayDB[ 26].insert(new G4QDecayChan(1.00,2212, 111));
	}
    if(limit<= 27 && nQ>= 27)DecayDB[ 27].insert(new G4QDecayChan(1.00,2212, 211)); // Delta ++
    if(limit<= 28 && nQ>= 28)    // Lambda* (1520)
	{
      DecayDB[ 28].insert(new G4QDecayChan(.230,3112,-311));
      DecayDB[ 28].insert(new G4QDecayChan(.460,3222,-321));
      DecayDB[ 28].insert(new G4QDecayChan(.463,3112,211,111));
      DecayDB[ 28].insert(new G4QDecayChan(.466,3212,211,-211));
      DecayDB[ 28].insert(new G4QDecayChan(.467,3212,111,111));
      DecayDB[ 28].insert(new G4QDecayChan(.470,3222,-211,111));
      DecayDB[ 28].insert(new G4QDecayChan(.540,3122,211,-211));
      DecayDB[ 28].insert(new G4QDecayChan(.570,3122,111,111));
      DecayDB[ 28].insert(new G4QDecayChan(.710,3222,-211));
      DecayDB[ 28].insert(new G4QDecayChan(.850,3212, 111));
      DecayDB[ 28].insert(new G4QDecayChan(.990,3112, 211));
      DecayDB[ 28].insert(new G4QDecayChan(1.00,3122,  22));
	}
    if(limit<= 29 && nQ>= 29)    // Sigma* -
	{
      DecayDB[ 29].insert(new G4QDecayChan(.060,3112, 111));
      DecayDB[ 29].insert(new G4QDecayChan(.120,3212,-211));
      DecayDB[ 29].insert(new G4QDecayChan(1.00,3122,-211));
	}
    if(limit<= 30 && nQ>= 30)    // Sigma* 0
	{
      DecayDB[ 30].insert(new G4QDecayChan(.040,3112, 211));
      DecayDB[ 30].insert(new G4QDecayChan(.080,3222,-211));
      DecayDB[ 30].insert(new G4QDecayChan(.120,3212, 111));
      DecayDB[ 30].insert(new G4QDecayChan(1.00,3122, 111));
	}
    if(limit<= 31 && nQ>= 31)    // Sigma* +
	{
      DecayDB[ 31].insert(new G4QDecayChan(.060,3212, 211));
      DecayDB[ 31].insert(new G4QDecayChan(.120,3222, 111));
      DecayDB[ 31].insert(new G4QDecayChan(1.00,3122, 211));
	}
    if(limit<= 32 && nQ>= 32)    // Ksi* -
	{
      DecayDB[ 32].insert(new G4QDecayChan(.667,3322,-211));
      DecayDB[ 32].insert(new G4QDecayChan(1.00,3312, 111));
	}
    if(limit<= 33 && nQ>= 33)    // Ksi* 0
	{
      DecayDB[ 33].insert(new G4QDecayChan(.667,3312, 211));
      DecayDB[ 33].insert(new G4QDecayChan(1.00,3322, 111));
	}
    //if(limit<= 34 && nQ>= 34)    // OMEGA - (Weak)
	//{
    //  DecayDB[ 34].insert(new G4QDecayChan(.678,3122, 321));
    //  DecayDB[ 34].insert(new G4QDecayChan(.914,3322,-211));
    //  DecayDB[ 34].insert(new G4QDecayChan(1.00,3312, 111));
	//}
    if(limit<= 35 && nQ>= 35)    // a_2 0
	{
      DecayDB[ 35].insert(new G4QDecayChan(.068, 211,-211,223));
      DecayDB[ 35].insert(new G4QDecayChan(.102, 111, 111,223));
      DecayDB[ 35].insert(new G4QDecayChan(.126, 321,-321));
      DecayDB[ 35].insert(new G4QDecayChan(.150, 311,-311));
      DecayDB[ 35].insert(new G4QDecayChan(.298, 111, 221));
      DecayDB[ 35].insert(new G4QDecayChan(.532,-211, 213));
      DecayDB[ 35].insert(new G4QDecayChan(.766, 211,-213));
      DecayDB[ 35].insert(new G4QDecayChan(1.00, 111, 113));
	}
    if(limit<= 36 && nQ>= 36)    // a_2 +
	{
      DecayDB[ 36].insert(new G4QDecayChan(.102,111,211,223));
      DecayDB[ 36].insert(new G4QDecayChan(.150, 321,-311));
      DecayDB[ 36].insert(new G4QDecayChan(.298, 211, 221));
      DecayDB[ 36].insert(new G4QDecayChan(.649, 211, 113));
      DecayDB[ 36].insert(new G4QDecayChan(1.00, 111, 213));
	}
    if(limit<= 37 && nQ>= 37)    // f_2 0
	{
      DecayDB[ 37].insert(new G4QDecayChan(.004, 221, 221));
      DecayDB[ 37].insert(new G4QDecayChan(.027, 311,-311));
      DecayDB[ 37].insert(new G4QDecayChan(.050, 321,-321));
      DecayDB[ 37].insert(new G4QDecayChan(.122, 111, 113));
      DecayDB[ 37].insert(new G4QDecayChan(.125, 111, 221));
      DecayDB[ 37].insert(new G4QDecayChan(.153, 211,-211,113));
      DecayDB[ 37].insert(new G4QDecayChan(.717, 211,-211));
      DecayDB[ 37].insert(new G4QDecayChan(.999, 111, 111));
      DecayDB[ 37].insert(new G4QDecayChan(1.00,  22,  22));
	}
    if(limit<= 38 && nQ>= 38)    // K_2 0
	{
      DecayDB[ 38].insert(new G4QDecayChan(.028, 311, 223));
      DecayDB[ 38].insert(new G4QDecayChan(.074, 211,-211,313));
      DecayDB[ 38].insert(new G4QDecayChan(.143,111,-211,323));
      DecayDB[ 38].insert(new G4QDecayChan(.166,111, 111,313));
      DecayDB[ 38].insert(new G4QDecayChan(.190,-211, 323));
      DecayDB[ 38].insert(new G4QDecayChan(.314, 111, 313));
      DecayDB[ 38].insert(new G4QDecayChan(.357, 311, 113));
      DecayDB[ 38].insert(new G4QDecayChan(.500, 321,-213));
      DecayDB[ 38].insert(new G4QDecayChan(.750,-211, 321));
      DecayDB[ 38].insert(new G4QDecayChan(1.00, 111, 311));
	}
    if(limit<= 39 && nQ>= 39)    // K_2 +
	{
      DecayDB[ 39].insert(new G4QDecayChan(.028, 321, 223));
      DecayDB[ 39].insert(new G4QDecayChan(.074,211,-211,323));
      DecayDB[ 39].insert(new G4QDecayChan(.143,111, 211,313));
      DecayDB[ 39].insert(new G4QDecayChan(.166,111, 111,323));
      DecayDB[ 39].insert(new G4QDecayChan(.190, 211, 313));
      DecayDB[ 39].insert(new G4QDecayChan(.314, 111, 323));
      DecayDB[ 39].insert(new G4QDecayChan(.357, 311, 213));
      DecayDB[ 39].insert(new G4QDecayChan(.500, 321, 113));
      DecayDB[ 39].insert(new G4QDecayChan(.750, 211, 311));
      DecayDB[ 39].insert(new G4QDecayChan(1.00, 111, 321));
	}
    if(limit<= 40 && nQ>= 40)    // f_2' 0
	{
      DecayDB[ 40].insert(new G4QDecayChan(.103, 221, 221));
      DecayDB[ 40].insert(new G4QDecayChan(.547, 311,-311));
      DecayDB[ 40].insert(new G4QDecayChan(.991, 321,-321));
      DecayDB[ 40].insert(new G4QDecayChan(.997, 211,-211));
      DecayDB[ 40].insert(new G4QDecayChan(1.00, 111, 111));
	}
    if(limit<= 41 && nQ>= 41)    // N_5/2 0
	{
      DecayDB[ 41].insert(new G4QDecayChan(.040, 211, 1114));
      DecayDB[ 41].insert(new G4QDecayChan(.080, 111, 2114));
      DecayDB[ 41].insert(new G4QDecayChan(.120,-211, 2214));
      DecayDB[ 41].insert(new G4QDecayChan(.180, 2112, 113));
      DecayDB[ 41].insert(new G4QDecayChan(.210, 2212,-213));
      DecayDB[ 41].insert(new G4QDecayChan(.340, 2112, 110));
      DecayDB[ 41].insert(new G4QDecayChan(.780, 2212,-211));
      DecayDB[ 41].insert(new G4QDecayChan(1.00, 2112, 111));
	}
    if(limit<= 42 && nQ>= 42)    // N_5/2 +
	{
      DecayDB[ 42].insert(new G4QDecayChan(.040,-211, 2224));
      DecayDB[ 42].insert(new G4QDecayChan(.080, 211, 2114));
      DecayDB[ 42].insert(new G4QDecayChan(.120, 111, 2214));
      DecayDB[ 42].insert(new G4QDecayChan(.180, 2112, 213));
      DecayDB[ 42].insert(new G4QDecayChan(.210, 2212, 113));
      DecayDB[ 42].insert(new G4QDecayChan(.340, 2212, 229));
      DecayDB[ 42].insert(new G4QDecayChan(.780, 2112, 211));
      DecayDB[ 42].insert(new G4QDecayChan(1.00, 2212, 111));
	}
    if(limit<= 43 && nQ>= 43)    // LAMBDA_5/2
	{
      DecayDB[ 43].insert(new G4QDecayChan(.350, 2112,-311));
      DecayDB[ 43].insert(new G4QDecayChan(.700, 2212,-321));
      DecayDB[ 43].insert(new G4QDecayChan(.740, 211, 3114));
      DecayDB[ 43].insert(new G4QDecayChan(.780,-211, 3224));
      DecayDB[ 43].insert(new G4QDecayChan(.820, 111, 3214));
      DecayDB[ 43].insert(new G4QDecayChan(.880, 3112, 211));
      DecayDB[ 43].insert(new G4QDecayChan(.940, 3222,-211));
      DecayDB[ 43].insert(new G4QDecayChan(1.00, 3212, 111));
	}
    if(limit<= 44 && nQ>= 44)    // SIGMA_5/2 -
	{
      DecayDB[ 44].insert(new G4QDecayChan(.600, 2112,-321));
      DecayDB[ 44].insert(new G4QDecayChan(.660,-211, 3214));
      DecayDB[ 44].insert(new G4QDecayChan(.720, 111, 3114));
      DecayDB[ 44].insert(new G4QDecayChan(.810, 3212,-211));
      DecayDB[ 44].insert(new G4QDecayChan(.900, 3112, 111));
      DecayDB[ 44].insert(new G4QDecayChan(1.00, 3122,-211));
	}
    if(limit<= 45 && nQ>= 45)    // SIGMA_5/2 0
	{
      DecayDB[ 45].insert(new G4QDecayChan(.300, 2112,-311));
      DecayDB[ 45].insert(new G4QDecayChan(.600, 2212,-321));
      DecayDB[ 45].insert(new G4QDecayChan(.640, 211, 3114));
      DecayDB[ 45].insert(new G4QDecayChan(.680,-211, 3224));
      DecayDB[ 45].insert(new G4QDecayChan(.720, 111, 3214));
      DecayDB[ 45].insert(new G4QDecayChan(.780, 3112, 211));
      DecayDB[ 45].insert(new G4QDecayChan(.840, 3222,-211));
      DecayDB[ 45].insert(new G4QDecayChan(.900, 3212, 111));
      DecayDB[ 45].insert(new G4QDecayChan(1.00, 3122, 111));
	}
    if(limit<= 46 && nQ>= 46)    // SIGMA_5/2 +
	{
      DecayDB[ 46].insert(new G4QDecayChan(.600, 2212,-311));
      DecayDB[ 46].insert(new G4QDecayChan(.660, 211, 3214));
      DecayDB[ 46].insert(new G4QDecayChan(.720, 111, 3224));
      DecayDB[ 46].insert(new G4QDecayChan(.810, 3212, 211));
      DecayDB[ 46].insert(new G4QDecayChan(.900, 3222, 111));
      DecayDB[ 46].insert(new G4QDecayChan(1.00, 3122, 211));
	}
    if(limit<= 47 && nQ>= 47)    // KSI_5/2 -
	{
      DecayDB[ 47].insert(new G4QDecayChan(.400, 3112,-311));
      DecayDB[ 47].insert(new G4QDecayChan(.800, 3212,-321));
      DecayDB[ 47].insert(new G4QDecayChan(1.00, 3122,-321));
	}
    if(limit<= 48 && nQ>= 48)    // KSI_5/2 0
	{
      DecayDB[ 48].insert(new G4QDecayChan(.400, 3212,-311));
      DecayDB[ 48].insert(new G4QDecayChan(.800, 3222,-321));
      DecayDB[ 48].insert(new G4QDecayChan(1.00, 3122,-311));
	}
    if(limit<= 49 && nQ>= 49)    // rho_3 0
	{
      DecayDB[ 49].insert(new G4QDecayChan(.019,311,-313));
      DecayDB[ 49].insert(new G4QDecayChan(.038,321,-323));
      DecayDB[ 49].insert(new G4QDecayChan(.046,311,-311));
      DecayDB[ 49].insert(new G4QDecayChan(.054,321,-321));
      DecayDB[ 49].insert(new G4QDecayChan(.224,111, 223));
      DecayDB[ 49].insert(new G4QDecayChan(.404,111,-211,213));
      DecayDB[ 49].insert(new G4QDecayChan(.584,111, 211,-213));
      DecayDB[ 49].insert(new G4QDecayChan(.764,111, 111,113));
      DecayDB[ 49].insert(new G4QDecayChan(1.00,211,-211));
	}
    if(limit<= 50 && nQ>= 50)    // rho_3 +
	{
      DecayDB[ 50].insert(new G4QDecayChan(.019, 321,-313));
      DecayDB[ 50].insert(new G4QDecayChan(.038,-311, 323));
      DecayDB[ 50].insert(new G4QDecayChan(.054, 321,-311));
      DecayDB[ 50].insert(new G4QDecayChan(.224, 211, 223));
      DecayDB[ 50].insert(new G4QDecayChan(.404,211,-211,213));
      DecayDB[ 50].insert(new G4QDecayChan(.584,211,211,-213));
      DecayDB[ 50].insert(new G4QDecayChan(.764,211,111,113));
      DecayDB[ 50].insert(new G4QDecayChan(1.00, 211, 111));
	}
    if(limit<= 51 && nQ>= 51)    // omega_3
	{
      DecayDB[ 51].insert(new G4QDecayChan(.020,211,-211,223));
      DecayDB[ 51].insert(new G4QDecayChan(.040,111, 111,223));
      DecayDB[ 51].insert(new G4QDecayChan(.060, 211,-213));
      DecayDB[ 51].insert(new G4QDecayChan(.080,-211, 213));
      DecayDB[ 51].insert(new G4QDecayChan(1.00, 111, 113));
	}
    if(limit<= 52 && nQ>= 52)    // K_3 0
	{
      DecayDB[ 52].insert(new G4QDecayChan(.030, 111, 315));
      DecayDB[ 52].insert(new G4QDecayChan(.060,-211, 325));
      DecayDB[ 52].insert(new G4QDecayChan(.340, 311, 331));
      DecayDB[ 52].insert(new G4QDecayChan(.400, 111, 313));
      DecayDB[ 52].insert(new G4QDecayChan(.520,-211, 323));
      DecayDB[ 52].insert(new G4QDecayChan(.620, 311, 113));
      DecayDB[ 52].insert(new G4QDecayChan(.820, 321,-213));
      DecayDB[ 52].insert(new G4QDecayChan(.940,-211, 321));
      DecayDB[ 52].insert(new G4QDecayChan(1.00, 111, 311));
	}
    if(limit<= 53 && nQ>= 53)    // K_3 +
	{
      DecayDB[ 53].insert(new G4QDecayChan(.030, 211, 315));
      DecayDB[ 53].insert(new G4QDecayChan(.060, 111, 325));
      DecayDB[ 53].insert(new G4QDecayChan(.340, 321, 331));
      DecayDB[ 53].insert(new G4QDecayChan(.400, 211, 313));
      DecayDB[ 53].insert(new G4QDecayChan(.520, 111, 323));
      DecayDB[ 53].insert(new G4QDecayChan(.620, 311, 213));
      DecayDB[ 53].insert(new G4QDecayChan(.820, 321, 113));
      DecayDB[ 53].insert(new G4QDecayChan(.940, 211, 311));
      DecayDB[ 53].insert(new G4QDecayChan(1.00, 111, 321));
	}
    if(limit<= 54 && nQ>= 54)    // phi_3
	{
      DecayDB[ 54].insert(new G4QDecayChan(.250, 321,-321));
      DecayDB[ 54].insert(new G4QDecayChan(.500, 311,-311));
      DecayDB[ 54].insert(new G4QDecayChan(.625, 321,-323));
      DecayDB[ 54].insert(new G4QDecayChan(.750,-321, 323));
      DecayDB[ 54].insert(new G4QDecayChan(.875, 311,-313));
      DecayDB[ 54].insert(new G4QDecayChan(1.00,-311, 313));
	}
    if(limit<= 55 && nQ>= 55)    // DELTA_7/2 -
	{
      DecayDB[ 55].insert(new G4QDecayChan(.200, 2112,-213 ));
      DecayDB[ 55].insert(new G4QDecayChan(.320,-211 , 2114));
      DecayDB[ 55].insert(new G4QDecayChan(.500, 111 , 1114));
      DecayDB[ 55].insert(new G4QDecayChan(1.00, 2112,-211 ));
	}
    if(limit<= 56 && nQ>= 56)    // DELTA_7/2 0
	{
      DecayDB[ 56].insert(new G4QDecayChan(.133, 2112, 113 ));
      DecayDB[ 56].insert(new G4QDecayChan(.200, 2212,-213 ));
      DecayDB[ 56].insert(new G4QDecayChan(.360,-211 , 2214));
      DecayDB[ 56].insert(new G4QDecayChan(.480, 211 , 1114));
      DecayDB[ 56].insert(new G4QDecayChan(.500, 111 , 2114));
      DecayDB[ 56].insert(new G4QDecayChan(.666, 2212,-211 ));
      DecayDB[ 56].insert(new G4QDecayChan(1.00, 2112, 111 ));
	}
    if(limit<= 57 && nQ>= 57)    // DELTA_7/2 +
	{
      DecayDB[ 57].insert(new G4QDecayChan(.133, 2112, 213 ));
      DecayDB[ 57].insert(new G4QDecayChan(.200, 2212, 113 ));
      DecayDB[ 57].insert(new G4QDecayChan(.360,-211 , 2224));
      DecayDB[ 57].insert(new G4QDecayChan(.480, 211 , 2114));
      DecayDB[ 57].insert(new G4QDecayChan(.500, 111 , 2214));
      DecayDB[ 57].insert(new G4QDecayChan(.666, 2112, 211 ));
      DecayDB[ 57].insert(new G4QDecayChan(1.00, 2212, 111 ));
	}
    if(limit<= 58 && nQ>= 58)    // DELTA_7/2 ++
	{
      DecayDB[ 58].insert(new G4QDecayChan(.200, 2212, 213 ));
      DecayDB[ 58].insert(new G4QDecayChan(.320, 211 , 2214));
      DecayDB[ 58].insert(new G4QDecayChan(.500, 111 , 2224));
      DecayDB[ 58].insert(new G4QDecayChan(1.00, 2112, 211 ));
	}
    if(limit<= 59 && nQ>= 59)     // LAMBDA_7/2
	{
      DecayDB[ 59].insert(new G4QDecayChan(.160, 3122, 223 ));
      DecayDB[ 59].insert(new G4QDecayChan(.260, 2112, 313 ));
      DecayDB[ 59].insert(new G4QDecayChan(.360, 2212, 323 ));
      DecayDB[ 59].insert(new G4QDecayChan(.400, 3312, 321 ));
      DecayDB[ 59].insert(new G4QDecayChan(.440, 3322, 311 ));
      DecayDB[ 59].insert(new G4QDecayChan(.480, 3122, 221 ));
      DecayDB[ 59].insert(new G4QDecayChan(.520, 2112, 311 ));
      DecayDB[ 59].insert(new G4QDecayChan(.560, 2212, 321 ));
      DecayDB[ 59].insert(new G4QDecayChan(.600, 3112, 211 ));
      DecayDB[ 59].insert(new G4QDecayChan(.800, 3222,-211 ));
      DecayDB[ 59].insert(new G4QDecayChan(1.00, 3212, 111 ));
	}
    if(limit<= 60 && nQ>= 60)    // SIGMA_7/2 -
	{
      DecayDB[ 60].insert(new G4QDecayChan(.030, 2112,-323 ));
      DecayDB[ 60].insert(new G4QDecayChan(.165,-311 , 1114));
      DecayDB[ 60].insert(new G4QDecayChan(.210,-321 , 2114));
      DecayDB[ 60].insert(new G4QDecayChan(.390,-211 , 3124));
      DecayDB[ 60].insert(new G4QDecayChan(.450,-211 , 3214));
      DecayDB[ 60].insert(new G4QDecayChan(.510, 111 , 3114));
      DecayDB[ 60].insert(new G4QDecayChan(.540, 311 , 3314));
      DecayDB[ 60].insert(new G4QDecayChan(.570,-211 , 3212));
      DecayDB[ 60].insert(new G4QDecayChan(.600, 111 , 3112));
      DecayDB[ 60].insert(new G4QDecayChan(.780,-321 , 2112));
      DecayDB[ 60].insert(new G4QDecayChan(1.00,-211 , 3122));
	}
    if(limit<= 61 && nQ>= 61)    // SIGMA_7/2 0
	{
      DecayDB[ 61].insert(new G4QDecayChan(.015, 2112,-313 ));
      DecayDB[ 61].insert(new G4QDecayChan(.030, 2212,-321 ));
      DecayDB[ 61].insert(new G4QDecayChan(.120,-311 , 2114));
      DecayDB[ 61].insert(new G4QDecayChan(.210,-321 , 2214));
      DecayDB[ 61].insert(new G4QDecayChan(.390, 111 , 3124));
      DecayDB[ 61].insert(new G4QDecayChan(.450,-211 , 3224));
      DecayDB[ 61].insert(new G4QDecayChan(.510, 211 , 3114));
      DecayDB[ 61].insert(new G4QDecayChan(.525, 311 , 3324));
      DecayDB[ 61].insert(new G4QDecayChan(.540, 321 , 3314));
      DecayDB[ 61].insert(new G4QDecayChan(.570,-211 , 3222));
      DecayDB[ 61].insert(new G4QDecayChan(.600, 211 , 3112));
      DecayDB[ 61].insert(new G4QDecayChan(.690,-311 , 2112));
      DecayDB[ 61].insert(new G4QDecayChan(.780,-321 , 2212));
      DecayDB[ 61].insert(new G4QDecayChan(1.00, 111 , 3122));
	}
    if(limit<= 62 && nQ>= 62)    // SIGMA_7/2 +
	{
      DecayDB[ 62].insert(new G4QDecayChan(.030, 2212,-313 ));
      DecayDB[ 62].insert(new G4QDecayChan(.165,-321 , 2224));
      DecayDB[ 62].insert(new G4QDecayChan(.210,-311 , 2214));
      DecayDB[ 62].insert(new G4QDecayChan(.390, 211 , 3124));
      DecayDB[ 62].insert(new G4QDecayChan(.450, 211 , 3214));
      DecayDB[ 62].insert(new G4QDecayChan(.510, 111 , 3224));
      DecayDB[ 62].insert(new G4QDecayChan(.540, 321 , 3324));
      DecayDB[ 62].insert(new G4QDecayChan(.570, 211 , 3212));
      DecayDB[ 62].insert(new G4QDecayChan(.600, 111 , 3222));
      DecayDB[ 62].insert(new G4QDecayChan(.780,-311 , 2212));
      DecayDB[ 62].insert(new G4QDecayChan(1.00, 211 , 3122));
	}
    if(limit<= 63 && nQ>= 63)    // KSI_7/2 -
	{
      DecayDB[ 63].insert(new G4QDecayChan(.400, 3112,-311));
      DecayDB[ 63].insert(new G4QDecayChan(.800, 3212,-321));
      DecayDB[ 63].insert(new G4QDecayChan(1.00, 3122,-321));
	}
    if(limit<= 64 && nQ>= 64)    // KSI_7/2 0
	{
      DecayDB[ 64].insert(new G4QDecayChan(.400, 3212,-311));
      DecayDB[ 64].insert(new G4QDecayChan(.800, 3222,-321));
      DecayDB[ 64].insert(new G4QDecayChan(1.00, 3122,-311));
	}
    if(limit<= 65 && nQ>= 65)    // OMEGA_7/2 -
	{
      DecayDB[ 65].insert(new G4QDecayChan(.250, 311 ,-3314));
      DecayDB[ 65].insert(new G4QDecayChan(.500, 321 ,-3324));
      DecayDB[ 65].insert(new G4QDecayChan(.750,-3312, 313 ));
      DecayDB[ 65].insert(new G4QDecayChan(1.00,-3322, 323 ));
	}
    if(limit<= 66 && nQ>= 66)    // a_4 0
	{
      DecayDB[ 66].insert(new G4QDecayChan(.200, 311,-311));
      DecayDB[ 66].insert(new G4QDecayChan(.400, 321,-321));
      DecayDB[ 66].insert(new G4QDecayChan(.600, 111, 221));
      DecayDB[ 66].insert(new G4QDecayChan(1.00, 111, 211,-211));
	}
    if(limit<= 67 && nQ>= 67)    // a_4 +
	{
      DecayDB[ 67].insert(new G4QDecayChan(.400, 321,-311));
      DecayDB[ 67].insert(new G4QDecayChan(.600, 211, 221));
      DecayDB[ 67].insert(new G4QDecayChan(.800, 211, 211,-211));
      DecayDB[ 67].insert(new G4QDecayChan(1.00, 211, 111, 111));
	}
    if(limit<= 68 && nQ>= 68)    // f_4 0
	{
      DecayDB[ 68].insert(new G4QDecayChan(.020, 333, 333));
      DecayDB[ 68].insert(new G4QDecayChan(.340, 223, 223));
      DecayDB[ 68].insert(new G4QDecayChan(.350, 221, 221));
      DecayDB[ 68].insert(new G4QDecayChan(.360, 311,-311));
      DecayDB[ 68].insert(new G4QDecayChan(.370, 321,-321));
      DecayDB[ 68].insert(new G4QDecayChan(.670, 213,-213));
      DecayDB[ 68].insert(new G4QDecayChan(.820, 113, 113));
      DecayDB[ 68].insert(new G4QDecayChan(.940, 211,-211));
      DecayDB[ 68].insert(new G4QDecayChan(1.00, 111, 111));
	}
    if(limit<= 69 && nQ>= 69)    // K_4 0
	{
      DecayDB[ 69].insert(new G4QDecayChan(.060, 333, 313));
      DecayDB[ 69].insert(new G4QDecayChan(.260, 223, 313));
      DecayDB[ 69].insert(new G4QDecayChan(.380, 313, 113));
      DecayDB[ 69].insert(new G4QDecayChan(.400, 323,-213));
      DecayDB[ 69].insert(new G4QDecayChan(.500, 223, 311));
      DecayDB[ 69].insert(new G4QDecayChan(.550, 311, 113));
      DecayDB[ 69].insert(new G4QDecayChan(.600, 321,-213));
      DecayDB[ 69].insert(new G4QDecayChan(.700, 311, 113));
      DecayDB[ 69].insert(new G4QDecayChan(.800, 321,-213));
      DecayDB[ 69].insert(new G4QDecayChan(.900, 311, 111));
      DecayDB[ 69].insert(new G4QDecayChan(1.00, 321,-211));
	}
    if(limit<= 70 && nQ>= 70)    // K_4 +
	{
      DecayDB[ 70].insert(new G4QDecayChan(.060, 333, 323));
      DecayDB[ 70].insert(new G4QDecayChan(.260, 223, 323));
      DecayDB[ 70].insert(new G4QDecayChan(.380, 313, 213));
      DecayDB[ 70].insert(new G4QDecayChan(.400, 323, 113));
      DecayDB[ 70].insert(new G4QDecayChan(.500, 223, 321));
      DecayDB[ 70].insert(new G4QDecayChan(.550, 321, 113));
      DecayDB[ 70].insert(new G4QDecayChan(.600, 311, 213));
      DecayDB[ 70].insert(new G4QDecayChan(.700, 311, 211));
      DecayDB[ 70].insert(new G4QDecayChan(.800, 321, 111));
      DecayDB[ 70].insert(new G4QDecayChan(.900, 311, 211));
      DecayDB[ 70].insert(new G4QDecayChan(1.00, 321, 111));
	}
    if(limit<= 71 && nQ>= 71)DecayDB[ 71].insert(new G4QDecayChan(1.00, 333, 333));//phi_4(2300)
    if(limit<= 75 && nQ>= 75)DecayDB[ 75].insert(new G4QDecayChan(1.00, 2112, 2112)); //nn
    if(limit<= 76 && nQ>= 76)DecayDB[ 76].insert(new G4QDecayChan(1.00, 2212, 2112)); //pn
    if(limit<= 77 && nQ>= 77)DecayDB[ 77].insert(new G4QDecayChan(1.00, 2212, 2212)); //pp
    // ------- Nuclear fragments
    //if(limit<= 72 && nQ>=72)
	//{
    //  if(limit<72) limit=72;
    //  for (int i=limit; i<nQ; i++) DecayDB[i] = NULL;
    //}

	//Update the limit
    limit=nQ+1;
#ifdef debug
	cout<<"G4QParticle::InitDecayVector: limit is set to "<<limit<<endl;
#endif
  }
  return DecayDB[abs(nQ)];
}

// Initialize the Particle by a Q Code
void G4QParticle::InitQParticle(G4int theQCode)
//    =============================================
{
  aQPDG.InitByQCode(theQCode);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(aQPDG.GetQCode());
}

// Initialize the Particle by a PDG Code
void G4QParticle::InitPDGParticle(G4int thePDGCode)
//    =============================================
{
  aQPDG      = G4QPDGCode(thePDGCode);
  aQuarkCont = aQPDG.GetQuarkContent();
  aDecay     = InitDecayVector(aQPDG.GetQCode());
}



