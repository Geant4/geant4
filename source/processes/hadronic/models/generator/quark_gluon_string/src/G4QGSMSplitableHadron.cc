#include "G4QGSMSplitableHadron.hh"
#include "G4ParticleTable.hh" 
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

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
  theMinPz = 0.8*G4PionMinus::PionMinus()->GetPDGMass(); 
//  theMinPz = 0.1*G4PionMinus::PionMinus()->GetPDGMass(); 
//  theMinPz = G4PionMinus::PionMinus()->GetPDGMass(); 
  // as low as possible, otherwise, we have unphysical boundary conditions in the sampling.
  StrangeSuppress = 0.48;
  sigmaPt = 0.*GeV; // widens eta slightly, if increased to 1.7, 
                    // but Maxim's algorithm breaks energy conservation
		    // to be revised.
  widthOfPtSquare = 0.01*GeV*GeV;
  Direction = FALSE;
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

G4QGSMSplitableHadron::~G4QGSMSplitableHadron(){}

const G4QGSMSplitableHadron & G4QGSMSplitableHadron::operator=(const G4QGSMSplitableHadron &right)
{
  G4Exception("G4QGSMSplitableHadron::operator= meant to not be accessable");
  return *this;
}


//**************************************************************************************************************************
    
void G4QGSMSplitableHadron::SplitUp()
{  
  if (!Color.isEmpty()) return;
  if (GetSoftCollisionCount() == 0)
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

  // momenta of string ends 
  G4double pt2 = HadronMom.perp2();
  G4double transverseMass2 = HadronMom.plus()*HadronMom.minus();
  G4double maxAvailMomentum2 = sqr(sqrt(transverseMass2) - sqrt(pt2));
  G4ThreeVector pt = GaussianPt(widthOfPtSquare, maxAvailMomentum2);

  G4LorentzVector LeftMom(pt, 0.);
  G4LorentzVector RightMom;
  RightMom.setPx(HadronMom.px() - pt.x());
  RightMom.setPy(HadronMom.py() - pt.y());

  G4double Local1 = HadronMom.minus() + (RightMom.perp2() - LeftMom.perp2())/HadronMom.plus();
  G4double Local2 = sqrt(G4std::max(0., sqr(Local1) - 4.*RightMom.perp2()*HadronMom.minus()/HadronMom.plus()));
  if (Direction) Local2 = -Local2;
  G4double RightMinus   = 0.5*(Local1 + Local2);
  G4double LeftMinus = HadronMom.minus() - RightMinus;

  G4double LeftPlus  = LeftMom.perp2()/LeftMinus;
  G4double RightPlus = HadronMom.plus() - LeftPlus;
  LeftMom.setPz(0.5*(LeftPlus - LeftMinus));
  LeftMom.setE (0.5*(LeftPlus + LeftMinus));
  RightMom.setPz(0.5*(RightPlus - RightMinus));
  RightMom.setE (0.5*(RightPlus + RightMinus));
  Left->Set4Momentum(LeftMom);
  Right->Set4Momentum(RightMom);
  Color.insert(Left);
  AntiColor.insert(Right);
}


void G4QGSMSplitableHadron::SoftSplitUp()
{
   //... sample transversal momenta for sea and valence quarks
   G4double phi, pts;
   G4double SumPy = 0.;
   G4double SumPx = 0.;
   G4ThreeVector Pos    = GetPosition();
   G4int nSeaPair = GetSoftCollisionCount()-1; 
   G4int aSeaPair;
   for (aSeaPair = 0; aSeaPair < nSeaPair; aSeaPair++)
   {
     G4int aPDGCode = 1 + (G4int)(G4UniformRand()/StrangeSuppress);

     G4Parton * aParton = BuildSeaQuark(false, aPDGCode, nSeaPair);
     SumPx += aParton->Get4Momentum().px();
     SumPy += aParton->Get4Momentum().py();
     Color.insert(aParton);

     aParton = BuildSeaQuark(true, aPDGCode, nSeaPair);
     SumPx += aParton->Get4Momentum().px();
     SumPy += aParton->Get4Momentum().py();
     AntiColor.insert(aParton);
   }
   // Valence quark    
   G4Parton* pColorParton = NULL;   
   G4Parton* pAntiColorParton = NULL;   
   GetValenceQuarkFlavors(GetDefinition(), pColorParton, pAntiColorParton);
   G4int ColorEncoding = pColorParton->GetPDGcode();
   G4int AntiColorEncoding = pAntiColorParton->GetPDGcode();
   
   pts   =  sigmaPt*sqrt(-log(G4UniformRand()));
   phi   = 2.*pi*G4UniformRand();
   G4double Px = pts*cos(phi);
   G4double Py = pts*sin(phi);
   SumPx += Px;
   SumPy += Py;

   if (ColorEncoding < 0) // use particle definition
   {
     G4LorentzVector ColorMom(-SumPx, -SumPy, 0, 0);
     pColorParton->Set4Momentum(ColorMom);
     G4LorentzVector AntiColorMom(Px, Py, 0, 0);
     pAntiColorParton->Set4Momentum(AntiColorMom);
   }
   else
   {
     G4LorentzVector ColorMom(Px, Py, 0, 0);
     pColorParton->Set4Momentum(ColorMom);
     G4LorentzVector AntiColorMom(-SumPx, -SumPy, 0, 0);
     pAntiColorParton->Set4Momentum(AntiColorMom);
   }
   Color.insert(pColorParton);
   AntiColor.insert(pAntiColorParton);

   // Sample X
   G4double LightConeMomentum = (Direction)? Get4Momentum().plus() : Get4Momentum().minus();
   G4double Xmin = theMinPz/LightConeMomentum;
   G4int nAttempt = 0;
   G4double SumX = 0;
   G4double aBeta = beta;
   G4double ColorX, AntiColorX;
   G4double HPWtest = 0;
   if (GetDefinition() == G4PionMinus::PionMinusDefinition()) aBeta = 1.;        
   if (GetDefinition() == G4PionPlus::PionPlusDefinition()) aBeta = 1.;     
   if (GetDefinition() == G4PionZero::PionZeroDefinition()) aBeta = 1.;       
   if (GetDefinition() == G4KaonPlus::KaonPlusDefinition()) aBeta = 0.;       
   if (GetDefinition() == G4KaonMinus::KaonMinusDefinition()) aBeta = 0.;       
   do
   {
     SumX = 0;
     nAttempt++;
     G4int    NumberOfUnsampledSeaQuarks = 2*nSeaPair;
     G4double beta1 = beta;
     if (abs(ColorEncoding) <= 1000 && abs(AntiColorEncoding) <= 1000) beta1 = 1.; //...  in a meson        
     ColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
     HPWtest = ColorX;
     while (ColorX < Xmin || ColorX > 1.|| 1. -  ColorX <= Xmin); 
     Color.last()->SetX(SumX = ColorX);// this is the valenz quark.
     for(G4int aPair = 0; aPair < nSeaPair; aPair++) 
     {
       NumberOfUnsampledSeaQuarks--;
       ColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
       Color.at(aPair)->SetX(ColorX);
       SumX += ColorX; 
       NumberOfUnsampledSeaQuarks--;
       AntiColorX = SampleX(Xmin, NumberOfUnsampledSeaQuarks, 2*nSeaPair, aBeta);
       AntiColor.at(aPair)->SetX(AntiColorX); // the 'sea' partons
       SumX += AntiColorX;
       if (1. - SumX <= Xmin)  break;
     }
   } 
   while (1. - SumX <= Xmin); 
   AntiColor.last()->SetX(1. - SumX); // the di-quark takes the rest, then go to momentum
   G4double lightCone = ((!Direction) ? Get4Momentum().minus() : Get4Momentum().plus());
   for(aSeaPair = 0; aSeaPair < nSeaPair+1; aSeaPair++) 
   {
     G4Parton* aParton = Color.at(aSeaPair);
     aParton->DefineMomentumInZ(lightCone, Direction);

     aParton = AntiColor.at(aSeaPair); 
     aParton->DefineMomentumInZ(lightCone, Direction);
   }  
//--DEBUG--   cout <<G4endl<<"XSAMPLE "<<HPWtest<<G4endl;
   return;
} 

void G4QGSMSplitableHadron::GetValenceQuarkFlavors(const G4ParticleDefinition * aPart, G4Parton *& Parton1, G4Parton *& Parton2)
{
   // Note! convention aEnd = q or qqbar and bEnd = qbar or qq.
  G4int aEnd;
  G4int bEnd;
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
}
 
    
G4ThreeVector G4QGSMSplitableHadron::GaussianPt(G4double widthSquare, G4double maxPtSquare)
{
  G4double R;
  while((R = -widthSquare*log(G4UniformRand())) > maxPtSquare);
  R = sqrt(R);
  G4double phi = twopi*G4UniformRand();
  return G4ThreeVector (R*cos(phi), R*sin(phi), 0.);    
}

G4Parton * G4QGSMSplitableHadron::
BuildSeaQuark(G4bool isAntiQuark, G4int aPDGCode, G4int nSeaPair)
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
  for(G4int ii=0; ii<100; ii++)
  {
    G4double y = pow(1./G4double(ii), alpha);
    y *= pow( pow(1-anXmin-totalSea*anXmin, alpha+1) - pow(anXmin, alpha+1), nSea);
    y *= pow(1-anXmin-totalSea*anXmin, aBeta+1) - pow(anXmin, aBeta+1);
    if(y>ymax) ymax = y;
  }
  G4double y;
  do
  {
    x1 = -1.;
    while(x1<anXmin||x1>=1-(totalSea+1)*anXmin) x1 = G4UniformRand();
    y = pow(x1, alpha);
    y *= pow( pow(1-x1-totalSea*anXmin, alpha+1) - pow(anXmin, alpha+1), nSea);
    y *= pow(1-x1-totalSea*anXmin, aBeta+1) - pow(anXmin, aBeta+1);  
    x2 = ymax*G4UniformRand();
  }
  while(x2>y);
  result = x1;
  return result;  
}
