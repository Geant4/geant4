// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first partial implementation, Maxim Komogorov, 9-Oct-1998
//
//      new OO design and implementation, largely replaces the 
//      initial prototype, H.P. Wellisch, Dec 1999.
// -----------------------------------------------------------------------------
#include "G4PartonStringAnnihilator.hh" 


G4PartonStringAnnihilator::G4PartonStringAnnihilator()
{
  widthOfPtSquare = 0.1*GeV*GeV; 
}


G4bool G4PartonStringAnnihilator::isMeson(G4int Encoding)
{
  return (abs(Encoding)%10000 < 1000) && abs(Encoding) > 100;
}

G4bool G4PartonStringAnnihilator::isBarion(G4int Encoding)
{
  return (abs(Encoding)%10000 > 1000) && abs(Encoding)%100 > 10;
}

G4bool G4PartonStringAnnihilator::isAntiParticle(G4int Encoding)
{
  return (Encoding < 0);
}


G4bool G4PartonStringAnnihilator::GetStringEnds(G4int Projectile, G4int Target, G4int* Left, G4int* Right)
{
   G4int LeftProjectile, RightProjectile;
   G4int LeftTarget,     RightTarget;
   G4bool result = FALSE;
   
   if (isBarion(Projectile) && isMeson(Target))
       G4SwapObj(&Projectile, &Target);
 
    // meson - meson annihilation
   if (isMeson(Projectile) && isMeson(Target))
   {
     theMesonSplitter.SplitMeson(Projectile, &LeftProjectile, &RightProjectile); 
     theMesonSplitter.SplitMeson(Target,     &LeftTarget,     &RightTarget); 
     if  (LeftProjectile == -RightTarget)
     {
       if (RightProjectile == -LeftTarget && G4UniformRand() < 0.5) 
       {
         *Left = LeftProjectile; 
         *Right = RightTarget;
       }
       else
       { 
         *Left  = LeftTarget;  
         *Right = RightProjectile;
       }
       result = TRUE;
     }        
     else if (-RightProjectile == LeftTarget)
     {
       *Left = LeftProjectile; 
       *Right = RightTarget;
       result = TRUE;
     }
   }
   
   // meson - barion(antibarion) annihilation
   else if (isMeson(Projectile) && isBarion(Target))
   {
     theMesonSplitter.SplitMeson (Projectile, &LeftProjectile, &RightProjectile); 
     G4int Diquark;
     if (isAntiParticle(Target))
     {
       if (theBaryonSplitter.FindDiquark(Target, LeftProjectile, &Diquark)) // assumes that LeftProjectile from mesonsplitter is quark.
       {
         *Left  = -Diquark;
         *Right = RightProjectile;
         result = TRUE;
       }            
     }
     else if (theBaryonSplitter.FindDiquark(Target, RightProjectile, &Diquark))
     {
       *Left  = LeftProjectile;
       *Right = Diquark;
       result = TRUE;
     }            
   }
         
   // barion-antibarion annihilation
   if (isAntiParticle(Projectile))
       G4SwapObj(&Projectile, &Target);
   if (!(isAntiParticle(Projectile) || !isAntiParticle(Target))) result = TRUE;
       
   if(!result) return FALSE;

   
   const G4SPBaryon & theProjectile = theBaryonSplitter.GetSPBaryon(Projectile);
   const G4SPBaryon & theTarget     = theBaryonSplitter.GetSPBaryon(Target);

   // Assumes that di-quarks always annihilate. Good only for nucleon
   // annihilation; not for strange baryons, etc.. (mass is higher)
   G4int activeDiquark=0;
   
   *Left = theProjectile.MatchDiQuarkAndGetQuark(theTarget, activeDiquark);
   
   *Right = theTarget.FindQuark(activeDiquark);
   
   return result;
}


//**************************************************************************************************************************

G4ExcitedString* G4PartonStringAnnihilator::GetString(G4KineticTrack& Target, G4KineticTrack& Projectile)
{
  G4LorentzVector Mom = Target.Get4Momentum() +  Projectile.Get4Momentum();
  G4ThreeVector   Pos = Target.GetPosition();
  
  // Convert the hadrons into one highmass colorneutral object with axis.

  G4int LeftEncoding;
  G4int RightEncoding;
  G4bool canAnnihilate = GetStringEnds(Projectile.GetDefinition()->GetPDGEncoding(), 
                                       Target.GetDefinition()->GetPDGEncoding(), 
                                      &LeftEncoding, &RightEncoding);
  if (!canAnnihilate) return 0;
  
  G4Parton* Left = new G4Parton(LeftEncoding);
  Left->SetPosition(Pos);

  G4Parton* Right = new G4Parton(RightEncoding);
  Right->SetPosition(Pos);

  // momenta of string ends
  G4double pt2 = Mom.perp2();
  G4double transverseMass2 = Mom.plus()*Mom.minus();
  G4double maxAvailMomentum2 = sqr(sqrt(transverseMass2) - sqrt(pt2));

  G4double R;
  while((R = -widthOfPtSquare*log(G4UniformRand())) > maxAvailMomentum2);
  R = sqrt(R);
  G4double phi = twopi*G4UniformRand();
  G4ThreeVector pt(R*cos(phi), R*sin(phi), 0.);    

  G4LorentzVector LeftMom(pt, 0.);
  G4LorentzVector RightMom;
  RightMom.setPx(Mom.px() - pt.x());
  RightMom.setPy(Mom.py() - pt.y());

  G4double Local1 = Mom.minus() + (RightMom.perp2() - LeftMom.perp2())/Mom.plus();
  G4double Local2 = sqrt(G4std::max(0., sqr(Local1) - 4.*RightMom.perp2()*Mom.minus()/Mom.plus()));
  G4double RightMinus   = 0.5*(Local1 + Local2);
  G4double LeftMinus = Mom.minus() - RightMinus;

  G4double LeftPlus  = LeftMom.perp2()/LeftMinus;
  G4double RightPlus = Mom.plus() - LeftPlus;
  LeftMom.setPz(0.5*(LeftPlus - LeftMinus));
  LeftMom.setE (0.5*(LeftPlus + LeftMinus));
  RightMom.setPz(0.5*(RightPlus - RightMinus));
  RightMom.setE (0.5*(RightPlus + RightMinus));

  Left->Set4Momentum(LeftMom);
  Right->Set4Momentum(RightMom);
  G4int aDirection=1;
  if(G4UniformRand()<0.5) aDirection = -1;
  return new G4ExcitedString(Left, Right, aDirection); 
}

//**************************************************************************************************************************







