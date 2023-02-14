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
  //  theMinPz = 0.1*G4PionMinus::PionMinus()->GetPDGMass();
  //  theMinPz = G4PionMinus::PionMinus()->GetPDGMass();
  // as low as possible, otherwise, we have unphysical boundary conditions in the sampling.
  StrangeSuppress = 0.48;
  sigmaPt = 0.*GeV; // widens eta slightly, if increased to 1.7,
  // but Maxim's algorithm breaks energy conservation to be revised.
  widthOfPtSquare = 0.5*sqr(GeV);
  Direction = FALSE;
  minTransverseMass = 1*keV;
  iP  =0;//     Color.begin();
  iAP =0;// AntiColor.begin();
} 

G4QGSMSplitableHadron::G4QGSMSplitableHadron() 
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary, G4bool aDirection)
:G4VSplitableHadron(aPrimary)
{
  InitParameters();
  Direction = aDirection;
}


G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4ReactionProduct & aPrimary)
:  G4VSplitableHadron(aPrimary)
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4Nucleon & aNucleon)
:  G4VSplitableHadron(aNucleon)
{
  InitParameters();
}

G4QGSMSplitableHadron::G4QGSMSplitableHadron(const G4Nucleon & aNucleon, G4bool aDirection)
:  G4VSplitableHadron(aNucleon)
{
  InitParameters();
  Direction = aDirection;
}

G4QGSMSplitableHadron::~G4QGSMSplitableHadron() {}


//**************************************************************************************************************************

void G4QGSMSplitableHadron::SplitUp()
{  
  if (IsSplit()) return;
  Splitting();  // To mark that a hadron is split
  if (Color.size()!=0) return;
  if (GetSoftCollisionCount() == 0)  // GetSoftCollisionCount() from G4VSplitableHadron
  {
    DiffractiveSplitUp();
  }
  else
  {
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

  G4double maxAvailMomentum2 = sqr(HadronMom.mag()/2.);

  G4ThreeVector pt(minTransverseMass, minTransverseMass, 0);
  if (maxAvailMomentum2/widthOfPtSquare>0.01) pt = GaussianPt(widthOfPtSquare, maxAvailMomentum2);

  G4LorentzVector LeftMom(pt, 0.);
  G4LorentzVector RightMom;
  RightMom.setPx(HadronMom.px() - pt.x());
  RightMom.setPy(HadronMom.py() - pt.y());

  G4double Local1 = HadronMom.minus() + (RightMom.perp2() - LeftMom.perp2())/HadronMom.plus();
  G4double Local2 = std::sqrt(std::max(0., sqr(Local1) - 4.*RightMom.perp2()*HadronMom.minus()/HadronMom.plus()));

  if (Direction) Local2 = -Local2;   
  G4double RightMinus   = 0.5*(Local1 + Local2);
  G4double LeftMinus = HadronMom.minus() - RightMinus;

  if (LeftMinus <= 0.) {
    RightMinus   = 0.5*(Local1 - Local2);
    LeftMinus = HadronMom.minus() - RightMinus;
  }

  G4double LeftPlus  = LeftMom.perp2()/LeftMinus;
  G4double RightPlus = HadronMom.plus() - LeftPlus;

  LeftMom.setPz(0.5*(LeftPlus - LeftMinus));
  LeftMom.setE (0.5*(LeftPlus + LeftMinus));
  RightMom.setPz(0.5*(RightPlus - RightMinus));
  RightMom.setE (0.5*(RightPlus + RightMinus));

  Left->Set4Momentum(LeftMom);
  Right->Set4Momentum(RightMom);

  Color.push_back(Left);
  AntiColor.push_back(Right);
  iP=0; iAP=0;
}


void G4QGSMSplitableHadron::SoftSplitUp()
{
  G4int nSeaPair = GetSoftCollisionCount()-1;

  G4LorentzVector tmp(0., 0., 0., 0.);

  G4int aSeaPair;
  for (aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
  {
    //  choose quark flavour, d:u:s = 1:1:(1/StrangeSuppress-2)
    G4int aPDGCode = 1 + (G4int)(G4UniformRand()/StrangeSuppress);

    //  BuildSeaQuark() determines quark spin, isospin and colour
    //  via parton-constructor G4Parton(aPDGCode)
    G4Parton * aParton = BuildSeaQuark(false, aPDGCode, nSeaPair);

    G4int firstPartonColour = aParton->GetColour();
    G4double firstPartonSpinZ = aParton->GetSpinZ();

    aParton->Set4Momentum(tmp);
    Color.push_back(aParton);

    // create anti-quark
    aParton = BuildSeaQuark(true, aPDGCode, nSeaPair);
    aParton->SetSpinZ(-firstPartonSpinZ);
    aParton->SetColour(-firstPartonColour);
    AntiColor.push_back(aParton);
  }

  // Valence quark
  G4Parton* pColorParton = NULL;
  G4Parton* pAntiColorParton = NULL;
  GetValenceQuarkFlavors(GetDefinition(), pColorParton, pAntiColorParton);

  pColorParton->Set4Momentum(tmp);
  pAntiColorParton->Set4Momentum(tmp);

  Color.push_back(pColorParton);
  AntiColor.push_back(pAntiColorParton);

  iP=0; iAP=0;

  return;
} 


void G4QGSMSplitableHadron::GetValenceQuarkFlavors(const G4ParticleDefinition * aPart, 
                                                   G4Parton *& Parton1, G4Parton *& Parton2)
{
  // Note! convention aEnd = q or (qq)bar and bEnd = qbar or qq.
  G4int aEnd=0;
  G4int bEnd=0;
  G4int HadronEncoding = aPart->GetPDGEncoding();
  if (aPart->GetBaryonNumber() == 0)
  {
    theMesonSplitter.SplitMeson(HadronEncoding, &aEnd, &bEnd);
  }
  else
  {
    theBaryonSplitter.SplitBarion(HadronEncoding, &aEnd, &bEnd);
  }

  Parton1 = new G4Parton(aEnd);
  Parton1->SetPosition(GetPosition());

  Parton2 = new G4Parton(bEnd);
  Parton2->SetPosition(GetPosition());

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
}


G4ThreeVector G4QGSMSplitableHadron::GaussianPt(G4double widthSquare, G4double maxPtSquare)
{
  G4double R;
  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = -1;
  while( ((R = -widthSquare*G4Log(G4UniformRand())) > maxPtSquare) && 
         ++loopCounter < maxNumberOfLoops ) {;}  /* Loop checking, 07.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    R = 0.99*maxPtSquare;  // Just an acceptable value, without any physics consideration.
  }
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
    y *= G4Pow::GetInstance()->powN( G4Pow::GetInstance()->powA(1-anXmin-totalSea*anXmin, alpha+1) - 
                                     G4Pow::GetInstance()->powA(anXmin, alpha+1), nSea );
    y *= G4Pow::GetInstance()->powA(1-anXmin-totalSea*anXmin, aBeta+1) - 
         G4Pow::GetInstance()->powA(anXmin, aBeta+1);
   if (y>ymax) ymax = y;
  }
  G4double y;
  G4double xMax=1-(totalSea+1)*anXmin;
  if (anXmin > xMax)
  {
    throw G4HadronicException(__FILE__, __LINE__, 
            "G4QGSMSplitableHadron - Fatal: Cannot sample parton densities under these constraints.");
  }
  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = 0;
  do
  {
    x1 = G4RandFlat::shoot(anXmin, xMax);
    y = G4Pow::GetInstance()->powA(x1, alpha);
    y *= G4Pow::GetInstance()->powN( G4Pow::GetInstance()->powA(1-x1-totalSea*anXmin, alpha+1) - 
                                     G4Pow::GetInstance()->powA(anXmin, alpha+1), nSea );
    y *= G4Pow::GetInstance()->powA(1-x1-totalSea*anXmin, aBeta+1) -
         G4Pow::GetInstance()->powA(anXmin, aBeta+1);
    x2 = ymax*G4UniformRand();
  } while( (x2>y) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) {
    x1 = 0.5*( anXmin + xMax );  // Just an acceptable value, without any physics consideration.         
  }
  result = x1;
  return result;
}
