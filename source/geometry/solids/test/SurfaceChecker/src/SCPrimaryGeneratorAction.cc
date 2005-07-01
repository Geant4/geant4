//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: SCPrimaryGeneratorAction.cc,v 1.3 2005-07-01 12:13:51 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Randomize.hh"

#include "SCPrimaryGeneratorAction.hh"
#include "SCDetectorConstruction.hh"

#include "SCSurfacePoint.hh" 
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SurfacePoint spoint ;

SCPrimaryGeneratorAction::SCPrimaryGeneratorAction(
                                               SCDetectorConstruction* myDC)
:myDetector(myDC)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
  
  particleGun->SetParticleDefinition(particle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCPrimaryGeneratorAction::~SCPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* aParticleDefinition 
    = particleTable->FindParticle(particleName="geantino");

  // create a new vertex in a random position 

  G4int i = anEvent->GetEventID() ;

  G4ThreeVector VertexPosition (GetRandomPosition());
  G4PrimaryVertex* aVertex = new G4PrimaryVertex( VertexPosition, 0);

  G4double u,v ;

  G4ThreeVector aSurfacePoint (GetSurfacePoint(u,v)) ;  // point A on surface
  G4ThreeVector aNormal(1,0,0) ;  // normal vector

  G4ThreeVector m = (aSurfacePoint-VertexPosition).unit() ;  // direction

  G4double distance = ( aSurfacePoint - VertexPosition ).mag() ;
  G4double theta = std::acos(m.unit()*aNormal) ;

  spoint.SetSurfacePoint(aSurfacePoint) ;
  
  G4PrimaryParticle* aPrimaryParticle =
    new G4PrimaryParticle(aParticleDefinition, m.x(), m.y(), m.z());
  aPrimaryParticle->SetMass (0.);

  aVertex->SetPrimary( aPrimaryParticle );

  G4cout  << "Event "     << i << G4endl 
	  << "Vertex "  << VertexPosition << " " << u << " " << v <<  G4endl 
	  << "Surface "  << aSurfacePoint << G4endl 
	  << "Distance " << distance << G4endl 
	  << "Momentum "  << m << G4endl 
          << "Angle " << theta 
	  << G4endl ;

  anEvent->AddPrimaryVertex( aVertex );


}


G4ThreeVector SCPrimaryGeneratorAction::GetRandomDirection() {

  G4ThreeVector retval;

  G4double CosTheta;
  G4double SinTheta;

  G4double Phi;
  G4double SinPhi;
  G4double CosPhi;

  G4double rand;

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}

G4ThreeVector SCPrimaryGeneratorAction::GetRandomPosition() 
{

  G4double a = 0.5*myDetector->GetWorldFullLength();

  G4double x = ( G4UniformRand()*2 - 1 )*a;
  G4double y = ( G4UniformRand()*2 - 1 )*a;
  G4double z = ( G4UniformRand()*2 - 1 )*a;

  G4ThreeVector retval (x, y, z);

  return retval;
}

G4ThreeVector SCPrimaryGeneratorAction::GetSurfacePoint(G4double &u, G4double &v)
{

  G4String val = myDetector->GetDetectorType() ;

  G4ThreeVector retval  ;

  if (val == "Sphere")
  {

  // get parameters

    G4double r1 = myDetector->GetTrackerR1() ;   // inner surface
    G4double r2 = myDetector->GetTrackerR2() ;   // outer surface
    G4double r = r1 ;    // we check the inner surface

    if ( r == 0 )       // but if no inner surface, then check outer instead.
      r = r2 ;  

    G4double phistart = myDetector->GetPhi() ;
    G4double phiseg   = myDetector->GetPhiSegment() ;

    G4double thetastart = myDetector->GetTheta() ;
    G4double thetaseg   = myDetector->GetThetaSegment() ;

    // Attention: theta goes from 0 to 180deg, but 
    // cos from 1 to -1. Thats why I inverse the intervall for 
    // the random number generation.
    G4double CosTheta1 = std::cos(thetastart+thetaseg) ;
    G4double CosTheta2   = std::cos(thetastart) ;
    G4double CosThetaSeg   = CosTheta2 - CosTheta1 ;  

    G4double CosTheta;
    G4double SinTheta;

    G4double Phi;
    G4double SinPhi;
    G4double CosPhi;

    G4double rand;

    rand = G4UniformRand();

    CosTheta = CosTheta1 + rand * CosThetaSeg;
    SinTheta = sqrt (1.-CosTheta*CosTheta);

    rand = G4UniformRand();
    Phi = phistart + rand * phiseg ;

    SinPhi = std::sin (Phi);
    CosPhi = std::cos (Phi);

    retval.setX(r*SinTheta*CosPhi);
    retval.setY(r*SinTheta*SinPhi);
    retval.setZ(r*CosTheta);

    u = Phi ;        // "rename" parameters
    v = CosTheta ;

#if 0 
    G4cout << "Sphere: " << G4endl <<
      " r1,r2 = " << r1 << ", " << r2 << G4endl <<
      " phistart, dphi = " << phistart  << ", " << phiseg << G4endl <<
      " theta, dtheta  = " << thetastart     << ", " << thetaseg << G4endl <<
      " CosTheta1, CosTheta2 = " << CosTheta1 << ", " << CosTheta2 << G4endl <<
      " Phi, CosTheta = " << Phi  << ", " << CosTheta << G4endl ;
#endif

  }
  else if (val == "Orb")
  {

    G4double r = myDetector->GetTrackerR() ;

    G4double SinTheta;
    G4double CosTheta;

    G4double rand;
    G4double Phi ;

    rand = G4UniformRand();

    CosTheta = 2.0*rand -1.0;
    SinTheta = sqrt (1.-CosTheta*CosTheta);

    rand = G4UniformRand();
    Phi = twopi*rand;

    retval.setX(r*SinTheta*std::cos(Phi));
    retval.setY(r*SinTheta*std::sin(Phi));
    retval.setZ(r*CosTheta);

    u = Phi ;     // rename parameteres
    v = CosTheta ;

  }
  else if (val == "Box") 
  {    

    u  = (G4UniformRand()*2 -1)*myDetector->GetTrackerpDy1() ;
    v   = (G4UniformRand()*2 -1)*myDetector->GetTrackerpDz() ;

    G4double fx  = myDetector->GetTrackerpDx1() ;

    retval.setX(fx) ;
    retval.setY(u) ;
    retval.setZ(v) ;

  }
  else if (val == "Cone")
  {        


    G4double r1 = myDetector->GetTrackerR1() ;
    G4double r2 = myDetector->GetTrackerR2() ;
    G4double pdz = myDetector->GetTrackerpDz() ;


    G4double Phi = twopi*G4UniformRand();
    G4double zpos = (2*G4UniformRand()-1)*pdz ;

    G4double r = 0.5*(r1+r2)-zpos*(r1-r2)/(2*pdz) ;

    retval.setX(r*std::cos(Phi));
    retval.setY(r*std::sin(Phi));
    retval.setZ(zpos);
    
    u = Phi ;
    v = zpos ;

  }
  else if (val == "manyCons")
  {        
      G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "Tube")
  {

    G4double r   = myDetector->GetTrackerR() ;
    G4double pdz = myDetector->GetTrackerpDz() ;

    G4double z = pdz*(2*G4UniformRand()-1);
    G4double phi = twopi*G4UniformRand();

    retval.setX(r*cos(phi)) ;
    retval.setY(r*sin(phi)) ;
    retval.setZ(z)  ;
    
    u = phi ;    // rename parameters
    v = z ;

  }
  else if (val == "Hype")
  {
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "Torus")
  {

    u = twopi*G4UniformRand();
    v = twopi*G4UniformRand();

    G4double c = myDetector->GetTrackerR() ;
    G4double a = myDetector->GetTrackerR2() ;   // take outer surface 

    retval.setX((c + a*std::cos(v))* std::cos(u));
    retval.setY((c + a*std::cos(v))* std::sin(u));
    retval.setZ( a * std::sin(v) );

  }
  else if (val == "Para")
  {
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "Trd")
  {
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Ub2") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Ib2") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Sb2") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Ib1") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Ub1") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "b1Sb1") 
  {         
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Not implemented!");
  }
  else if (val == "TwistedTubs")
  {
   
    // G4double r1 = myDetector->GetTrackerR1() ;  // not used 
    G4double r2 = myDetector->GetTrackerR2() ;
    G4double pdz = myDetector->GetTrackerpDz() ;
    G4double twistangle = myDetector->GetTwistAngle() ;
    G4double phisegment = myDetector->GetPhi() ;

    G4double zpos = (2*G4UniformRand()-1)*pdz ;
  // attention: twisted surface -> intervall of phi is z-dependent !

    // Corner x0 and direction d of boundary (from lower to uper endcap)
    // between hyperbolic and twisted surface of G4TwistedTubs
    // Corner0Min, Corner0Max are corners at lower endcap
    // Corner1Min, Corner1Max are corners at upper endcap

    G4ThreeVector Corner0Min( r2 * std::cos(-phisegment/2.-twistangle/2.) ,
			 r2 * std::sin(-phisegment/2.-twistangle/2.) , 
			 -pdz )  ; 
    G4ThreeVector Corner0Max( r2 * std::cos(phisegment/2.-twistangle/2.) , 
			 r2 * std::sin(phisegment/2.-twistangle/2.) ,
			 -pdz )  ;
    G4ThreeVector Corner1Min( r2 * std::cos(-phisegment/2.+twistangle/2.) ,
			 r2 * std::sin(-phisegment/2.+twistangle/2.) , 
			 pdz )  ;
    G4ThreeVector Corner1Max(r2 * std::cos(phisegment/2.+twistangle/2.) , 
			r2 * std::sin(phisegment/2.+twistangle/2.) ,
			 pdz )  ;


    G4ThreeVector dMin = (Corner1Min-Corner0Min).unit()  ;    // direction vectors
    G4ThreeVector dMax = (Corner1Max-Corner0Max).unit()  ;

    G4ThreeVector XMin =  ( zpos - Corner0Min.z() )/dMin.z() * dMin + Corner0Min ;  // coordinates of limits
    G4ThreeVector XMax =  ( zpos - Corner0Max.z() )/dMax.z() * dMax + Corner0Max ;  // at a given z position 

    G4double PhiMin  =  std::atan2( XMin.y(), XMin.x() ) ;  // limits in angle phi at given z position
    G4double PhiMax  =  std::atan2( XMax.y(), XMax.x() ) ;

    G4double Phi = G4UniformRand()*(PhiMax-PhiMin)+PhiMin ;  // generate a random angle between limits

    G4double tanAlpha = r2/pdz * std::sin(twistangle/2) ;
    G4double r0 = std::sqrt(r2*r2 - pdz*pdz*tanAlpha*tanAlpha) ;
    G4double r = std::sqrt( r0*r0 + zpos*zpos*tanAlpha*tanAlpha) ;

#if 0
    G4cout << "Corner0Min = " << Corner0Min    << G4endl ;
    G4cout << "z position = " << zpos          << G4endl ;
    G4cout << "XMin       = " << XMin          << G4endl ;
    G4cout << "XMax       = " << XMax          << G4endl ;
#endif

    retval.setX(r*std::cos(Phi));
    retval.setY(r*std::sin(Phi));
    retval.setZ(zpos);

    u = Phi ;   // rename
    v = zpos ;
    
  }

  else if (val == "TwistedBox")
  {

    G4double fDx1 = myDetector->GetTrackerpDx1() ;
    G4double fDx2 = fDx1 ;  // regular 
    G4double fDx3 = fDx1 ;
    G4double fDx4 = fDx1 ;

    G4double fDy1 = myDetector->GetTrackerpDy1() ;
    G4double fDy2 = fDy1 ;  // regular

    G4double fDz  = myDetector->GetTrackerpDz() ;
    G4double fPhiTwist = myDetector->GetTwistAngle() ;

    G4double fTAlph = 0.*deg ;

    G4double fDx4plus2  = fDx4 + fDx2 ;
    G4double fDx4minus2 = fDx4 - fDx2 ;
    G4double fDx3plus1  = fDx3 + fDx1 ; 
    G4double fDx3minus1 = fDx3 - fDx1 ;

    G4double  fDy2plus1  = fDy2 + fDy1 ;
    G4double  fDy2minus1 = fDy2 - fDy1 ;

    G4double fdeltaX = 0.  ;  // regular
    G4double fdeltaY = 0.  ;  // 


  // generate random point on surface
    G4double phi = ( G4UniformRand()*2 - 1 )*fPhiTwist/2;
    G4double Aphi = fDx4plus2 + fDx4minus2  * ( 2 * phi ) / fPhiTwist ;
    G4double Dphi = fDx3plus1 + fDx3minus1 * ( 2 * phi ) / fPhiTwist  ;
    G4double Bphi =  fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ;
    
    u = ( G4UniformRand()*2 - 1 ) * 0.5* Bphi  ;

  // calculate the cartesian position
    G4double w = Aphi/2. + (Dphi-Aphi)/4. - u*( ( Dphi-Aphi ) / ( 2 * Bphi ) + fTAlph )   ;
    G4double fx =  (w + fdeltaX*phi/fPhiTwist )*std::cos(phi) - (u+fdeltaY*phi/fPhiTwist)*std::sin(phi) ;
    G4double fy =  (w + fdeltaX*phi/fPhiTwist )*std::sin(phi) + (u+fdeltaY*phi/fPhiTwist)*std::cos(phi) ;
    G4double fz =  2*fDz*phi/fPhiTwist   ;

    retval.setX(fx) ;
    retval.setY(fy) ;
    retval.setZ(fz) ;

    v = phi ; //   "rename" parameter 


  }
  else if (val == "TwistedTrd")
  {

    G4double fDx1 = myDetector->GetTrackerpDx1() ;
    G4double fDx2 = fDx1 ;   // regular
    G4double fDx3 = myDetector->GetTrackerpDx2() ;  // attention !
    G4double fDx4 = fDx3 ;           


    G4double fDy1 = myDetector->GetTrackerpDy1() ;
    G4double fDy2 = myDetector->GetTrackerpDy2() ;

    G4double fDz  = myDetector->GetTrackerpDz() ;
  
    G4double fTAlph = 0 ;  // regular case
    G4double fPhi   = 0 ;
    G4double fTheta = 0 ;
 
    G4double fPhiTwist = myDetector->GetTwistAngle() ;

    G4double fDx4plus2  = fDx4 + fDx2 ;
    G4double fDx4minus2 = fDx4 - fDx2 ;
    G4double fDx3plus1  = fDx3 + fDx1 ; 
    G4double fDx3minus1 = fDx3 - fDx1 ;

    G4double  fDy2plus1  = fDy2 + fDy1 ;
    G4double  fDy2minus1 = fDy2 - fDy1 ;

    G4double fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;  // dx in surface equation
    G4double fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;  // dy in surface equation


  // generate random point on surface
    G4double phi = ( G4UniformRand()*2 - 1 )*fPhiTwist/2;
    G4double Aphi = fDx4plus2 + fDx4minus2  * ( 2 * phi ) / fPhiTwist ;
    G4double Dphi = fDx3plus1 + fDx3minus1 * ( 2 * phi ) / fPhiTwist  ;
    G4double Bphi =  fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ;
    
    u = ( G4UniformRand()*2 - 1 ) * 0.5* Bphi  ;

  // calculate the cartesian position
    G4double w = Aphi/2. + (Dphi-Aphi)/4. - u*( ( Dphi-Aphi ) / ( 2 * Bphi ) + fTAlph )   ;
    G4double fx =  (w + fdeltaX*phi/fPhiTwist )*std::cos(phi) - (u+fdeltaY*phi/fPhiTwist)*std::sin(phi) ;
    G4double fy =  (w + fdeltaX*phi/fPhiTwist )*std::sin(phi) + (u+fdeltaY*phi/fPhiTwist)*std::cos(phi) ;
    G4double fz =  2*fDz*phi/fPhiTwist   ;

    retval.setX(fx) ;
    retval.setY(fy) ;
    retval.setZ(fz) ;

    v = phi ; //   "rename" parameter 

  }
  else if (val == "TwistedTrap")
  {

    G4double fDx1 = myDetector->GetTrackerpDx1() ;
    G4double fDx2 = myDetector->GetTrackerpDx2() ;
    G4double fDx3 = fDx1 ;   // equal sized endcaps
    G4double fDx4 = fDx2 ;

    G4double fDy1 = myDetector->GetTrackerpDy1() ;
    G4double fDy2 = fDy1 ;   // equal sized endcaps

    G4double fDz  = myDetector->GetTrackerpDz() ;
  
    G4double fTAlph = 0 ;
    G4double fPhi   = 0 ;
    G4double fTheta = 0 ; // regular case
 
    G4double fPhiTwist = myDetector->GetTwistAngle() ;

    G4double fDx4plus2  = fDx4 + fDx2 ;
    G4double fDx4minus2 = fDx4 - fDx2 ;
    G4double fDx3plus1  = fDx3 + fDx1 ; 
    G4double fDx3minus1 = fDx3 - fDx1 ;

    G4double  fDy2plus1  = fDy2 + fDy1 ;
    G4double  fDy2minus1 = fDy2 - fDy1 ;

    G4double fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;  // dx in surface equation
    G4double fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;  // dy in surface equation


  // generate random point on surface
    G4double phi = ( G4UniformRand()*2 - 1 )*fPhiTwist/2;
    G4double Aphi = fDx4plus2 + fDx4minus2  * ( 2 * phi ) / fPhiTwist ;
    G4double Dphi = fDx3plus1 + fDx3minus1 * ( 2 * phi ) / fPhiTwist  ;
    G4double Bphi =  fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ;
    
    u = ( G4UniformRand()*2 - 1 ) * 0.5* Bphi  ;

  // calculate the cartesian position
    G4double w = Aphi/2. + (Dphi-Aphi)/4. - u*( ( Dphi-Aphi ) / ( 2 * Bphi ) + fTAlph )   ;
    G4double fx =  (w + fdeltaX*phi/fPhiTwist )*std::cos(phi) - (u+fdeltaY*phi/fPhiTwist)*std::sin(phi) ;
    G4double fy =  (w + fdeltaX*phi/fPhiTwist )*std::sin(phi) + (u+fdeltaY*phi/fPhiTwist)*std::cos(phi) ;
    G4double fz =  2*fDz*phi/fPhiTwist   ;

    retval.setX(fx) ;
    retval.setY(fy) ;
    retval.setZ(fz) ;

    v = phi ; //   "rename" parameter 



  }
  else if ( val == "TwistedTrap2") 
  {

    G4double fDx1 = myDetector->GetTrackerpDx1() ;
    G4double fDx2 = myDetector->GetTrackerpDx2() ;
    G4double fDx3 = myDetector->GetTrackerpDx3() ;
    G4double fDx4 = myDetector->GetTrackerpDx4() ;

    G4double fDy1 = myDetector->GetTrackerpDy1() ;
    G4double fDy2 = myDetector->GetTrackerpDy2() ;

    G4double fDz  = myDetector->GetTrackerpDz() ;
  
    G4double fAlph = -myDetector->GetAlpha() ;  // minus sign for surface equation
    G4double fTAlph = std::tan(fAlph) ;
    G4double fPhi   = myDetector->GetPhi() ;
    G4double fTheta = myDetector->GetTheta() ;
 
    G4double fPhiTwist = myDetector->GetTwistAngle() ;

    G4double fDx4plus2  = fDx4 + fDx2 ;
    G4double fDx4minus2 = fDx4 - fDx2 ;
    G4double fDx3plus1  = fDx3 + fDx1 ; 
    G4double fDx3minus1 = fDx3 - fDx1 ;

    G4double  fDy2plus1  = fDy2 + fDy1 ;
    G4double  fDy2minus1 = fDy2 - fDy1 ;

    G4double fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;  // dx in surface equation
    G4double fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;  // dy in surface equation


  // generate random point on surface
    G4double phi = ( G4UniformRand()*2 - 1 )*fPhiTwist/2;
    G4double Aphi = fDx4plus2 + fDx4minus2  * ( 2 * phi ) / fPhiTwist ;
    G4double Dphi = fDx3plus1 + fDx3minus1 * ( 2 * phi ) / fPhiTwist  ;
    G4double Bphi =  fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ;
    
    u = ( G4UniformRand()*2 - 1 ) * 0.5* Bphi  ;

  // calculate the cartesian position
    G4double w = Aphi/2. + (Dphi-Aphi)/4. - u*( ( Dphi-Aphi ) / ( 2 * Bphi ) + fTAlph )   ;
    G4double fx =  (w + fdeltaX*phi/fPhiTwist )*std::cos(phi) - (u+fdeltaY*phi/fPhiTwist)*std::sin(phi) ;
    G4double fy =  (w + fdeltaX*phi/fPhiTwist )*std::sin(phi) + (u+fdeltaY*phi/fPhiTwist)*std::cos(phi) ;
    G4double fz =  2*fDz*phi/fPhiTwist   ;

    retval.setX(fx) ;
    retval.setY(fy) ;
    retval.setZ(fz) ;

    v = phi ; //   "rename" parameter 

  }
  else if ( val == "TwistedTrap3") 
  {

    G4double fDx1 = myDetector->GetTrackerpDx1() ;
    G4double fDx2 = myDetector->GetTrackerpDx2() ;
    G4double fDx3 = myDetector->GetTrackerpDx3() ;
    G4double fDx4 = myDetector->GetTrackerpDx4() ;

    G4double fDy1 = myDetector->GetTrackerpDy1() ;
    G4double fDy2 = myDetector->GetTrackerpDy2() ;

    G4double fDz  = myDetector->GetTrackerpDz() ;
  
    G4double fAlph = -myDetector->GetAlpha() ;  // minus sign for surface equation
    G4double fTAlph = std::tan(fAlph) ;
    G4double fPhi   = myDetector->GetPhi() ;
    G4double fTheta = myDetector->GetTheta() ;
 
    G4double fPhiTwist = myDetector->GetTwistAngle() ;

    G4double fDx4plus2  = fDx4 + fDx2 ;
    G4double fDx4minus2 = fDx4 - fDx2 ;
    G4double fDx3plus1  = fDx3 + fDx1 ; 
    G4double fDx3minus1 = fDx3 - fDx1 ;

    G4double  fDy2plus1  = fDy2 + fDy1 ;
    G4double  fDy2minus1 = fDy2 - fDy1 ;

    G4double fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;  // dx in surface equation
    G4double fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;  // dy in surface equation


  // generate random point on surface
    G4double phi = ( G4UniformRand()*2 - 1 )*fPhiTwist/2;
    G4double Aphi = fDx4plus2 + fDx4minus2  * ( 2 * phi ) / fPhiTwist ;
    G4double Dphi = fDx3plus1 + fDx3minus1 * ( 2 * phi ) / fPhiTwist  ;
    G4double Bphi =  fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ;
    
    u = ( G4UniformRand()*2 - 1 ) * 0.5* Bphi  ;

  // calculate the cartesian position
    G4double w = Aphi/2. + (Dphi-Aphi)/4. - u*( ( Dphi-Aphi ) / ( 2 * Bphi ) + fTAlph )   ;
    G4double fx =  (w + fdeltaX*phi/fPhiTwist )*std::cos(phi) - (u+fdeltaY*phi/fPhiTwist)*std::sin(phi) ;
    G4double fy =  (w + fdeltaX*phi/fPhiTwist )*std::sin(phi) + (u+fdeltaY*phi/fPhiTwist)*std::cos(phi) ;
    G4double fz =  2*fDz*phi/fPhiTwist   ;

    retval.setX(fx) ;
    retval.setY(fy) ;
    retval.setZ(fz) ;

    v = phi ; //   "rename" parameter 

  }

  else
  {
     G4Exception("SCPrimaryGeneratorAction::GetSurfacePoint() - Invalid Shape!");
  }
 


  return retval ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

