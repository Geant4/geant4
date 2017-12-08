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
#include "G4QGSMSplitableHadron.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh" 
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

#include "G4Log.hh"
#include "G4Pow.hh"


// based on prototype by Maxim Komogorov
// Splitting into methods, and centralizing of model parameters HPW Feb 1999
// restructuring HPW Feb 1999
// fixing bug in the sampling of 'x', HPW Feb 1999
// fixing bug in sampling pz, HPW Feb 1999. 
// Code now also good for p-nucleus scattering (before only p-p), HPW Feb 1999.
// Using Parton more directly, HPW Feb 1999.
// Shortening the algorithm for sampling x, HPW Feb 1999.
// sampling of x replaced by formula, taking X_min into account in the correlated sampling. HPW, Feb 1999.
// logic much clearer now. HPW Feb 1999
// Removed the ordering problem. No Direction needed in selection of valence quark types. HPW Mar'99.
// Fixing p-t distributions for scattering of nuclei.
// Separating out parameters.

void G4QGSMSplitableHadron::InitParameters()
{
  // changing rapidity distribution for all
  alpha = -0.5; // Note that this number is still assumed in the algorithm
  // needs to be generalized.
  // changing rapidity distribution for projectile like
  beta = 2.5;// Note that this number is still assumed in the algorithm
  // needs to be generalized.
  theMinPz = 0.5*G4PionMinus::PionMinus()->GetPDGMass();
  //theMinPz = 0.1*G4PionMinus::PionMinus()->GetPDGMass();
  //theMinPz = G4PionMinus::PionMinus()->GetPDGMass();
  // as low as possible, otherwise, we have unphysical boundary conditions in the sampling.
  StrangeSuppress = 0.48;
  sigmaPt = 0.*GeV; // widens eta slightly, if increased to 1.7,
           // but Maxim's algorithm breaks energy conservation to be revised.
  widthOfPtSquare = 0.01*GeV*GeV;
  Direction = FALSE;
  minTransverseMass = 1*keV;
} 

G4QGSMSplitableHadron::G4QGSMSplitableHadron() 
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary, G4bool aDirection)
: G4VSplitableHadron(aPrimary)
{
  InitParameters();
  Direction = aDirection;
}


G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary)
: G4VSplitableHadron(aPrimary)
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4Nucleon & aNucleon)
: G4VSplitableHadron(aNucleon)
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4Nucleon & aNucleon, G4bool aDirection)
: G4VSplitableHadron(aNucleon)
{
  InitParameters();
  Direction = aDirection;
}

G4QGSMSplitableHadron::~G4QGSMSplitableHadron(){}

//**************************************************************************************************************************

void G4QGSMSplitableHadron::SplitUp()
{  
  if (IsSplit()) return;
  Splitting();
  if (Color.size()!=0) return;
  if (GetSoftCollisionCount() == 0)
  {
    DiffractiveSplitUp();
  } else {
    SoftSplitUp();
  }
}

void G4QGSMSplitableHadron::DiffractiveSplitUp()
{   
  // take the particle definitions and get the partons HPW
  G4Parton * Left = NULL;
  G4Parton * Right = NULL;
  GetValenceQuarkFlavors(GetDefinition(), Left, Right);
  Left->SetPosition(GetPosition());
  Right->SetPosition(GetPosition());

  G4LorentzVector HadronMom = Get4Momentum();
  //std::cout << "DSU 1 - "<<HadronMom<<std::endl;

  // momenta of string ends
  G4double pt2 = HadronMom.perp2();
  G4double transverseMass2 = HadronMom.plus()*HadronMom.minus();
  G4double maxAvailMomentum2 = sqr(std::sqrt(transverseMass2) - std::sqrt(pt2));
  G4ThreeVector pt(minTransverseMass, minTransverseMass, 0);
  if(maxAvailMomentum2/widthOfPtSquare>0.01) pt = GaussianPt(widthOfPtSquare, maxAvailMomentum2);
  //std::cout << "DSU 1.1 - "<< maxAvailMomentum2<< pt <<std::endl;

  G4LorentzVector LeftMom(pt, 0.);
  G4LorentzVector RightMom;
  RightMom.setPx(HadronMom.px() - pt.x());
  RightMom.setPy(HadronMom.py() - pt.y());
  //std::cout << "DSU 2 - "<<RightMom<<" "<< LeftMom <<std::endl;

  G4double Local1 = HadronMom.minus() + (RightMom.perp2() - LeftMom.perp2())/HadronMom.plus();
  G4double Local2 = std::sqrt(std::max(0., sqr(Local1) - 4.*RightMom.perp2()*HadronMom.minus()/HadronMom.plus()));
  //std::cout << "DSU 3 - "<< Local1 <<" "<< Local2 <<std::endl;
  if (Direction) Local2 = -Local2;
  G4double RightMinus   = 0.5*(Local1 + Local2);
  G4double LeftMinus = HadronMom.minus() - RightMinus;
  //std::cout << "DSU 4 - "<< RightMinus <<" "<< LeftMinus << " "<<HadronMom.minus() <<std::endl;

  G4double LeftPlus  = LeftMom.perp2()/LeftMinus;
  G4double RightPlus = HadronMom.plus() - LeftPlus;
  //std::cout << "DSU 5 - "<< RightPlus <<" "<< LeftPlus <<std::endl;
  LeftMom.setPz(0.5*(LeftPlus - LeftMinus));
  LeftMom.setE (0.5*(LeftPlus + LeftMinus));
  RightMom.setPz(0.5*(RightPlus - RightMinus));
  RightMom.setE (0.5*(RightPlus + RightMinus));
  //std::cout << "DSU 6 - "<< LeftMom <<" "<< RightMom <<std::endl;
  Left->Set4Momentum(LeftMom);
  Right->Set4Momentum(RightMom);
  Color.push_back(Left);
  AntiColor.push_back(Right);
}


void G4QGSMSplitableHadron::SoftSplitUp()
{
  //... sample transversal momenta for sea and valence quarks
  G4double phi, pts;
  G4double SumPy = 0.;
  G4double SumPx = 0.;
  G4ThreeVector Pos    = GetPosition();
  G4int nSeaPair = GetSoftCollisionCount()-1;

  // here the condition,to ensure viability of splitting, also in cases
  // where difractive excitation occured together with soft scattering.
  //   G4double LightConeMomentum = (Direction)? Get4Momentum().plus() : Get4Momentum().minus();
  //   G4double Xmin = theMinPz/LightConeMomentum;
  G4double Xmin = theMinPz/( Get4Momentum().e() - GetDefinition()->GetPDGMass() );
  while(Xmin>=1-(2*nSeaPair+1)*Xmin) Xmin*=0.95;  /* Loop checking, 26.10.2015, A.Ribon */

  G4int aSeaPair;
  for (aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
  {
    //  choose quark flavour, d:u:s = 1:1:(1/StrangeSuppress-2)

    G4int aPDGCode = 1 + (G4int)(G4UniformRand()/StrangeSuppress);

    //  BuildSeaQuark() determines quark spin, isospin and colour
    //  via parton-constructor G4Parton(aPDGCode)

    G4Parton * aParton = BuildSeaQuark(false, aPDGCode, nSeaPair);

    //G4cerr << "G4QGSMSplitableHadron::SoftSplitUp()" << G4endl;

    //G4cerr << "Parton 1: "
    //       << " PDGcode: "  << aPDGCode
    //       << " - Name: "   << aParton->GetDefinition()->GetParticleName()
    //       << " - Type: "   << aParton->GetDefinition()->GetParticleType()
    //       << " - Spin-3: " << aParton->GetSpinZ()
    //       << " - Colour: " << aParton->GetColour() << G4endl;

    // save colour a spin-3 for anti-quark

    G4int firstPartonColour = aParton->GetColour();
    G4double firstPartonSpinZ = aParton->GetSpinZ();

    SumPx += aParton->Get4Momentum().px();
    SumPy += aParton->Get4Momentum().py();
    Color.push_back(aParton);

    // create anti-quark

    aParton = BuildSeaQuark(true, aPDGCode, nSeaPair);
    aParton->SetSpinZ(-firstPartonSpinZ);
    aParton->SetColour(-firstPartonColour);

    //G4cerr << "Parton 2: "
    //       << " PDGcode: "  << -aPDGCode
    //	     << " - Name: "   << aParton->GetDefinition()->GetParticleName()
    //	     << " - Type: "   << aParton->GetDefinition()->GetParticleType()
    //	     << " - Spin-3: " << aParton->GetSpinZ()
    //	     << " - Colour: " << aParton->GetColour() << G4endl;
    //G4cerr << "------------" << G4endl;

    SumPx += aParton->Get4Momentum().px();
    SumPy += aParton->Get4Momentum().py();
    AntiColor.push_back(aParton);
  }

  // Valence quark
  G4Parton* pColorParton = NULL;
  G4Parton* pAntiColorParton = NULL;
  GetValenceQuarkFlavors(GetDefinition(), pColorParton, pAntiColorParton);
  G4int ColorEncoding = pColorParton->GetPDGcode();

  pts   =  sigmaPt*std::sqrt(-G4Log(G4UniformRand()));
  phi   = 2.*pi*G4UniformRand();
  G4double Px = pts*std::cos(phi);
  G4double Py = pts*std::sin(phi);
  SumPx += Px;
  SumPy += Py;

  if (ColorEncoding < 0) // use particle definition
  {
    G4LorentzVector ColorMom(-SumPx, -SumPy, 0, 0);
    pColorParton->Set4Momentum(ColorMom);
    G4LorentzVector AntiColorMom(Px, Py, 0, 0);
    pAntiColorParton->Set4Momentum(AntiColorMom);
  } else {
    G4LorentzVector ColorMom(Px, Py, 0, 0);
    pColorParton->Set4Momentum(ColorMom);
    G4LorentzVector AntiColorMom(-SumPx, -SumPy, 0, 0);
    pAntiColorParton->Set4Momentum(AntiColorMom);
  }
  Color.push_back(pColorParton);
  AntiColor.push_back(pAntiColorParton);

  // Sample X
  G4int nAttempt = 0;
  G4double SumX = 0;
  G4double aBeta = beta;
  G4double ColorX, AntiColorX;
  if (GetDefinition() == G4PionMinus::PionMinusDefinition()) aBeta = 1.;
  if (GetDefinition() == G4Gamma::GammaDefinition()) aBeta = 1.;
  if (GetDefinition() == G4PionPlus::PionPlusDefinition()) aBeta = 1.;
  if (GetDefinition() == G4PionZero::PionZeroDefinition()) aBeta = 1.;
  if (GetDefinition() == G4KaonPlus::KaonPlusDefinition()) aBeta = 0.;
  if (GetDefinition() == G4KaonMinus::KaonMinusDefinition()) aBeta = 0.;
  const G4int maxNumberOfAttempts = 1000;
  do
  {
    SumX = 0;
    nAttempt++;
    G4int NumberOfUnsampledSeaQuarks = 2*nSeaPair;
    ColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
    Color.back()->SetX(SumX = ColorX);// this is the valenz quark.
    for(G4int aPair = 0; aPair < nSeaPair; aPair++)
    {
      NumberOfUnsampledSeaQuarks--;
      ColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
      Color[aPair]->SetX(ColorX);
      SumX += ColorX;
      NumberOfUnsampledSeaQuarks--;
      AntiColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
      AntiColor[aPair]->SetX(AntiColorX); // the 'sea' partons
      SumX += AntiColorX;
      if (1. - SumX <= Xmin)  break;
    }
  } while ( (1. - SumX <= Xmin) && nAttempt < maxNumberOfAttempts );  /* Loop checking, 26.10.2015, A.Ribon */   
  if ( nAttempt >= maxNumberOfAttempts ) return;

  (*(AntiColor.end()-1))->SetX(1. - SumX); // the di-quark takes the rest, then go to momentum
  G4double lightCone  = ((!Direction) ? Get4Momentum().minus() : Get4Momentum().plus());
  G4double lightCone2 = ((!Direction) ? Get4Momentum().plus() : Get4Momentum().minus());
  for(aSeaPair = 0; aSeaPair < nSeaPair+1; aSeaPair++)
  {
    G4Parton* aParton = Color[aSeaPair];
    aParton->DefineMomentumInZ(lightCone, lightCone2, Direction);

    aParton = AntiColor[aSeaPair];
    aParton->DefineMomentumInZ(lightCone, lightCone2, Direction);
  }
  return;
} 

void G4QGSMSplitableHadron::
GetValenceQuarkFlavors(const G4ParticleDefinition * aPart, G4Parton *& Parton1, G4Parton *& Parton2)
{
  // Note! convention aEnd = q or (qq)bar and bEnd = qbar or qq.
  G4int aEnd;
  G4int bEnd;
  G4int HadronEncoding = aPart->GetPDGEncoding();
  if (aPart->GetBaryonNumber() == 0)
  {
    theMesonSplitter.SplitMeson(HadronEncoding, &aEnd, &bEnd);
  } else {
    theBaryonSplitter.SplitBarion(HadronEncoding, &aEnd, &bEnd);
  }

  Parton1 = new G4Parton(aEnd);
  Parton1->SetPosition(GetPosition());

  //G4cerr << "G4QGSMSplitableHadron::GetValenceQuarkFlavors()" << G4endl;
  //G4cerr << "Parton 1: "
  //	   << " PDGcode: "  << aEnd
  //	   << " - Name: "   << Parton1->GetDefinition()->GetParticleName()
  //	   << " - Type: "   << Parton1->GetDefinition()->GetParticleType()
  //	   << " - Spin-3: " << Parton1->GetSpinZ()
  //	   << " - Colour: " << Parton1->GetColour() << G4endl;

  Parton2 = new G4Parton(bEnd);
  Parton2->SetPosition(GetPosition());

  //G4cerr << "Parton 2: "
  //	   << " PDGcode: "  << bEnd
  //	   << " - Name: "   << Parton2->GetDefinition()->GetParticleName()
  //	   << " - Type: "   << Parton2->GetDefinition()->GetParticleType()
  //	   << " - Spin-3: " << Parton2->GetSpinZ()
  //	   << " - Colour: " << Parton2->GetColour() << G4endl;
  //G4cerr << "... now checking for color and spin conservation - yielding: " << G4endl;

  // colour of parton 1 choosen at random by G4Parton(aEnd)
  // colour of parton 2 is the opposite:

  Parton2->SetColour(-(Parton1->GetColour()));

  // isospin-3 of both partons is handled by G4Parton(PDGCode)

  // spin-3 of parton 1 and 2 choosen at random by G4Parton(aEnd)
  // spin-3 of parton 2 may be constrained by spin of original particle:

  if ( std::abs(Parton1->GetSpinZ() + Parton2->GetSpinZ()) > aPart->GetPDGSpin())
  {
    Parton2->SetSpinZ(-(Parton2->GetSpinZ()));    
  }

  //G4cerr << "Parton 2: "
  //       << " PDGcode: "  << bEnd
  //       << " - Name: "   << Parton2->GetDefinition()->GetParticleName()
  //	   << " - Type: "   << Parton2->GetDefinition()->GetParticleType()
  //	   << " - Spin-3: " << Parton2->GetSpinZ()
  //	   << " - Colour: " << Parton2->GetColour() << G4endl;
  //G4cerr << "------------" << G4endl;
}


G4ThreeVector G4QGSMSplitableHadron::GaussianPt(G4double widthSquare, G4double maxPtSquare)
{
  G4double R;
  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = -1;
  while ( ((R = -widthSquare*G4Log(G4UniformRand())) > maxPtSquare) &&  /* Loop checking, 26.10.2015, A.Ribon */
          ++loopCounter < maxNumberOfLoops ) {;}
  if ( loopCounter >= maxNumberOfLoops ) R = 0.0;
  R = std::sqrt(R);
  G4double phi = twopi*G4UniformRand();
  return G4ThreeVector (R*std::cos(phi), R*std::sin(phi), 0.);
}

G4Parton * G4QGSMSplitableHadron::
BuildSeaQuark(G4bool isAntiQuark, G4int aPDGCode, G4int /* nSeaPair*/)
{
  if (isAntiQuark) aPDGCode*=-1;
  G4Parton* result = new G4Parton(aPDGCode);
  result->SetPosition(GetPosition());
  G4ThreeVector aPtVector = GaussianPt(sigmaPt, DBL_MAX);
  G4LorentzVector a4Momentum(aPtVector, 0);
  result->Set4Momentum(a4Momentum);
  return result;
}

G4double G4QGSMSplitableHadron::
SampleX(G4double anXmin, G4int nSea, G4int totalSea, G4double aBeta)
{
  G4double result;
  G4double x1, x2;
  G4double ymax = 0;
  for(G4int ii=1; ii<100; ii++)
  {
    G4double y = G4Pow::GetInstance()->powA(1./G4double(ii), alpha);
    y *= G4Pow::GetInstance()->powN( G4Pow::GetInstance()->powA(1-anXmin-totalSea*anXmin, alpha+1) - G4Pow::GetInstance()->powA(anXmin, alpha+1), nSea);
    y *= G4Pow::GetInstance()->powA(1-anXmin-totalSea*anXmin, aBeta+1) - G4Pow::GetInstance()->powA(anXmin, aBeta+1);
    if(y>ymax) ymax = y;
  }
  G4double y;
  G4double xMax=1-(totalSea+1)*anXmin;
  if(anXmin > xMax)
  {
    G4cout << "anXmin = "<<anXmin<<" nSea = "<<nSea<<" totalSea = "<< totalSea<<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "G4QGSMSplitableHadron - Fatal: Cannot sample parton densities under these constraints.");
  }
  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = -1;
  do
  {
    x1 = G4RandFlat::shoot(anXmin, xMax);
    y = G4Pow::GetInstance()->powA(x1, alpha);
    y *= G4Pow::GetInstance()->powN( G4Pow::GetInstance()->powA(1-x1-totalSea*anXmin, alpha+1) - G4Pow::GetInstance()->powA(anXmin, alpha+1), nSea);
    y *= G4Pow::GetInstance()->powA(1-x1-totalSea*anXmin, aBeta+1) - G4Pow::GetInstance()->powA(anXmin, aBeta+1);
    x2 = ymax*G4UniformRand();
  } while ( (x2>y) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 26.10.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    G4ExceptionDescription ed;
    ed << " Failed sampling after maxNumberOfLoops attempts : forced exit! " << G4endl;
    G4Exception( "G4QGSMSplitableHadron::SampleX ", "HAD_QGS_002", JustWarning, ed );
  }
  result = x1;
  return result;
}

