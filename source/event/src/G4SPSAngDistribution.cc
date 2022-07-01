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
// G4SPSAngDistribution class implementation
//
// Author: Fan Lei, QinetiQ ltd. - 05/02/2004
// Customer: ESA/ESTEC
// Revisions: Andrea Dotti, SLAC
// --------------------------------------------------------------------

#include "G4SPSAngDistribution.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

G4SPSAngDistribution::G4SPSAngDistribution() 
{
  // Angular distribution Variables
  G4ThreeVector zero;
  particle_momentum_direction = G4ParticleMomentum(0,0,-1);

  AngDistType = "planar"; 
  AngRef1 = CLHEP::HepXHat;
  AngRef2 = CLHEP::HepYHat;
  AngRef3 = CLHEP::HepZHat;
  MinTheta = 0.;
  MaxTheta = pi;
  MinPhi = 0.;
  MaxPhi = twopi;
  DR = 0.;
  DX = 0.;
  DY = 0.;
  FocusPoint = G4ThreeVector(0., 0., 0.);
  UserDistType = "NULL";
  UserWRTSurface = true;
  UserAngRef = false;
  IPDFThetaExist = false;
  IPDFPhiExist = false;
  verbosityLevel = 0;

  G4MUTEXINIT(mutex);
}

G4SPSAngDistribution::~G4SPSAngDistribution()
{
    G4MUTEXDESTROY(mutex);
}

void G4SPSAngDistribution::SetAngDistType(const G4String& atype)
{
  G4AutoLock l(&mutex);
  if(atype != "iso" && atype != "cos" && atype != "user" && atype != "planar"
     && atype != "beam1d" && atype != "beam2d"  && atype != "focused")
  {
    G4cout << "Error, distribution must be iso, cos, planar, beam1d, beam2d, focused or user"
           << G4endl;
  }
  else
  {
    AngDistType = atype;
  }
  if (AngDistType == "cos")  { MaxTheta = pi/2.; }
  if (AngDistType == "user")
  {
    UDefThetaH = IPDFThetaH = ZeroPhysVector;
    IPDFThetaExist = false;
    UDefPhiH = IPDFPhiH = ZeroPhysVector;
    IPDFPhiExist = false;
  }
}

void G4SPSAngDistribution::DefineAngRefAxes(const G4String& refname,
                                            const G4ThreeVector& ref)
{
  G4AutoLock l(&mutex);
  if (refname == "angref1")
    AngRef1 = ref.unit(); // x'
  else if (refname == "angref2")
    AngRef2 = ref.unit(); // vector in x'y' plane

  // User defines x' (AngRef1) and a vector in the x'y'
  // plane (AngRef2). Then, AngRef1 x AngRef2 = AngRef3
  // the z' vector. Then, AngRef3 x AngRef1 = AngRef2
  // which will now be y'.

  AngRef3 = AngRef1.cross(AngRef2); // z'
  AngRef2 = AngRef3.cross(AngRef1); // y'
  UserAngRef = true ;
  if(verbosityLevel == 2)
  {
    G4cout << "Angular distribution rotation axes " << AngRef1
           << " " << AngRef2 << " " << AngRef3 << G4endl;
  }
}

void G4SPSAngDistribution::SetMinTheta(G4double mint)
{
  G4AutoLock l(&mutex);
  MinTheta = mint;
}

void G4SPSAngDistribution::SetMinPhi(G4double minp)
{
  G4AutoLock l(&mutex);
  MinPhi = minp;
}

void G4SPSAngDistribution::SetMaxTheta(G4double maxt)
{
  G4AutoLock l(&mutex);
  MaxTheta = maxt;
}

void G4SPSAngDistribution::SetMaxPhi(G4double maxp)
{
  G4AutoLock l(&mutex);
  MaxPhi = maxp;
}

void G4SPSAngDistribution::SetBeamSigmaInAngR(G4double r)
{
  G4AutoLock l(&mutex);
  DR = r;
}

void G4SPSAngDistribution::SetBeamSigmaInAngX(G4double r)
{
  G4AutoLock l(&mutex);
  DX = r;
}

void G4SPSAngDistribution::SetBeamSigmaInAngY(G4double r)
{
  G4AutoLock l(&mutex);
  DY = r;
}

void G4SPSAngDistribution::
SetParticleMomentumDirection(const G4ParticleMomentum& aMomentumDirection)
{
  G4AutoLock l(&mutex);
  particle_momentum_direction = aMomentumDirection.unit();
}

void G4SPSAngDistribution::SetPosDistribution(G4SPSPosDistribution* a)
{
  G4AutoLock l(&mutex);
  posDist = a;
}

void G4SPSAngDistribution::SetBiasRndm(G4SPSRandomGenerator* a)
{
  G4AutoLock l(&mutex);
  angRndm = a;
}

void G4SPSAngDistribution::SetVerbosity(G4int a)
{
  G4AutoLock l(&mutex);
  verbosityLevel = a;
}

void G4SPSAngDistribution::UserDefAngTheta(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  if(UserDistType == "NULL") UserDistType = "theta";
  if(UserDistType == "phi") UserDistType = "both";  
  G4double thi, val;
  thi = input.x();
  val = input.y();
  if(verbosityLevel >= 1) G4cout << "In UserDefAngTheta" << G4endl;
  UDefThetaH.InsertValues(thi, val);
}

G4String G4SPSAngDistribution::GetDistType()
{
  G4AutoLock l(&mutex);
  return AngDistType;
}

G4double G4SPSAngDistribution::GetMinTheta()
{
  G4AutoLock l(&mutex);
  return MinTheta;
}

G4double G4SPSAngDistribution::GetMaxTheta()
{
  G4AutoLock l(&mutex);
  return MaxTheta;
}

G4double G4SPSAngDistribution::GetMinPhi()
{
  G4AutoLock l(&mutex);
  return MinPhi;
}

G4double G4SPSAngDistribution::GetMaxPhi()
{
  G4AutoLock l(&mutex);
  return MaxPhi;
}

G4ThreeVector G4SPSAngDistribution::GetDirection()
{
  G4AutoLock l(&mutex);
  return particle_momentum_direction;
}

void G4SPSAngDistribution::UserDefAngPhi(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  if(UserDistType == "NULL") UserDistType = "phi";
  if(UserDistType == "theta") UserDistType = "both";  
  G4double phhi, val;
  phhi = input.x();
  val = input.y();
  if(verbosityLevel >= 1) G4cout << "In UserDefAngPhi" << G4endl;
  UDefPhiH.InsertValues(phhi, val); 
}

void G4SPSAngDistribution::SetFocusPoint(const G4ThreeVector& input)
{
  G4AutoLock l(&mutex);
  FocusPoint = input;
}

void G4SPSAngDistribution::SetUserWRTSurface(G4bool wrtSurf)
{
  G4AutoLock l(&mutex);

  // if UserWRTSurface = true then the user wants momenta with respect
  // to the surface normals.
  // When doing this theta has to be 0-90 only otherwise there will be
  // errors, which currently are flagged anywhere.
  //
  UserWRTSurface = wrtSurf;
}

void G4SPSAngDistribution::SetUseUserAngAxis(G4bool userang)
{
  G4AutoLock l(&mutex);

  // if UserAngRef = true  the angular distribution is defined wrt 
  // the user defined coordinates
  //
  UserAngRef = userang;
}

void G4SPSAngDistribution::GenerateBeamFlux(G4ParticleMomentum& mom)
{
  G4double theta, phi;
  G4double px, py, pz;
  if (AngDistType == "beam1d")
  { 
    theta = G4RandGauss::shoot(0.0,DR);
    phi = twopi * G4UniformRand();
  }
  else 
  {
    px = G4RandGauss::shoot(0.0,DX);
    py = G4RandGauss::shoot(0.0,DY);
    theta = std::sqrt (px*px + py*py);
    if (theta != 0.)
    { 
      phi = std::acos(px/theta);
      if ( py < 0.) phi = -phi;
    }
    else
    {
      phi = 0.0;
    }
  }
  px = -std::sin(theta) * std::cos(phi);
  py = -std::sin(theta) * std::sin(phi);
  pz = -std::cos(theta);
  G4double finx, finy,  finz;
  finx=px, finy=py, finz=pz;
  if (UserAngRef)
  {
    // Apply Angular Rotation Matrix
    // x * AngRef1, y * AngRef2 and z * AngRef3
    finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
    finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
    finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    G4double ResMag = std::sqrt((finx*finx) + (finy*finy) + (finz*finz));
    finx = finx/ResMag;
    finy = finy/ResMag;
    finz = finz/ResMag;
  }
  mom.setX(finx);
  mom.setY(finy);
  mom.setZ(finz);

  // particle_momentum_direction now holds unit momentum vector

  if(verbosityLevel >= 1)
  {
    G4cout << "Generating beam vector: " << mom << G4endl;
  }
}

void G4SPSAngDistribution::GenerateFocusedFlux(G4ParticleMomentum& mom)
{
  mom = (FocusPoint - posDist->GetParticlePos()).unit();

  // particle_momentum_direction now holds unit momentum vector.

  if(verbosityLevel >= 1)
  {
    G4cout << "Generating focused vector: " << mom << G4endl;
  }
}

void G4SPSAngDistribution::GenerateIsotropicFlux(G4ParticleMomentum& mom)
{
  // generates isotropic flux.
  // No vectors are needed.

  G4double rndm, rndm2;
  G4double px, py, pz;

  G4double sintheta, sinphi,costheta,cosphi;
  rndm = angRndm->GenRandTheta();
  costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta)
                                        - std::cos(MaxTheta));
  sintheta = std::sqrt(1. - costheta*costheta);
  
  rndm2 = angRndm->GenRandPhi();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = std::sin(Phi);
  cosphi = std::cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  // For volume and point source use mother or user defined coordinates
  // for plane and surface source user surface-normal or user-defined
  // coordinates
  //
  G4double finx, finy, finz;
  if (posDist->GetSourcePosType() == "Point"
   || posDist->GetSourcePosType() == "Volume")
  {
    if (UserAngRef)
    {
      // Apply Rotation Matrix
      // x * AngRef1, y * AngRef2 and z * AngRef3
      finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
      finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
      finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    }
    else
    {
      finx = px;
      finy = py;
      finz = pz;
    }
  }
  else
  {    // for plane and surface source   
    if (UserAngRef)
    {
      // Apply Rotation Matrix
      // x * AngRef1, y * AngRef2 and z * AngRef3
      finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
      finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
      finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    }
    else
    {
      finx = (px*posDist->GetSideRefVec1().x())
           + (py*posDist->GetSideRefVec2().x())
           + (pz*posDist->GetSideRefVec3().x());
      finy = (px*posDist->GetSideRefVec1().y())
           + (py*posDist->GetSideRefVec2().y())
           + (pz*posDist->GetSideRefVec3().y());
      finz = (px*posDist->GetSideRefVec1().z())
           + (py*posDist->GetSideRefVec2().z())
           + (pz*posDist->GetSideRefVec3().z());
    }
  }
  G4double ResMag = std::sqrt((finx*finx) + (finy*finy) + (finz*finz));
  finx = finx/ResMag;
  finy = finy/ResMag;
  finz = finz/ResMag;

  mom.setX(finx);
  mom.setY(finy);
  mom.setZ(finz);

  // particle_momentum_direction now holds unit momentum vector.

  if(verbosityLevel >= 1)
  {
    G4cout << "Generating isotropic vector: " << mom << G4endl;
  }
}

void G4SPSAngDistribution::GenerateCosineLawFlux(G4ParticleMomentum& mom)
{
  // Method to generate flux distributed with a cosine law

  G4double px, py, pz;
  G4double rndm, rndm2;
 
  G4double sintheta, sinphi,costheta,cosphi;
  rndm = angRndm->GenRandTheta();
  sintheta = std::sqrt( rndm * (std::sin(MaxTheta)*std::sin(MaxTheta)
                              - std::sin(MinTheta)*std::sin(MinTheta) ) 
                      + std::sin(MinTheta)*std::sin(MinTheta) );
  costheta = std::sqrt(1. -sintheta*sintheta);
  
  rndm2 = angRndm->GenRandPhi();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = std::sin(Phi);
  cosphi = std::cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  // for volume and point source use mother or user defined coordinates
  // for plane and surface source user surface-normal or userdefined
  // coordinates
  //
  G4double finx, finy, finz;
  if (posDist->GetSourcePosType() == "Point"
   || posDist->GetSourcePosType() == "Volume")
  {
    if (UserAngRef)
    {
      // Apply Rotation Matrix
      finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
      finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
      finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    }
    else
    {
      finx = px;
      finy = py;
      finz = pz;
    }
  }
  else
  {    // for plane and surface source   
    if (UserAngRef)
    {
      // Apply Rotation Matrix
      finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
      finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
      finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    }
    else
    {
      finx = (px*posDist->GetSideRefVec1().x())
           + (py*posDist->GetSideRefVec2().x())
           + (pz*posDist->GetSideRefVec3().x());
      finy = (px*posDist->GetSideRefVec1().y())
           + (py*posDist->GetSideRefVec2().y())
           + (pz*posDist->GetSideRefVec3().y());
      finz = (px*posDist->GetSideRefVec1().z())
           + (py*posDist->GetSideRefVec2().z())
           + (pz*posDist->GetSideRefVec3().z());
    }
  }
  G4double ResMag = std::sqrt((finx*finx) + (finy*finy) + (finz*finz));
  finx = finx/ResMag;
  finy = finy/ResMag;
  finz = finz/ResMag;

  mom.setX(finx);
  mom.setY(finy);
  mom.setZ(finz);

  // particle_momentum_direction now contains unit momentum vector.

  if(verbosityLevel >= 1)
  {
    G4cout << "Resultant cosine-law unit momentum vector " << mom << G4endl;
  }
}

void G4SPSAngDistribution::GeneratePlanarFlux(G4ParticleMomentum& mom)
{
  // particle_momentum_direction now contains unit momentum vector.
  // nothing need be done here as the m-directions have been set directly
  // under this option

  if(verbosityLevel >= 1)
  {
    G4cout << "Resultant Planar wave  momentum vector " << mom << G4endl;
  }
}

void G4SPSAngDistribution::GenerateUserDefFlux(G4ParticleMomentum& mom)
{
  G4double rndm, px, py, pz, pmag;

  if(UserDistType == "NULL")
  {
    G4cout << "Error: UserDistType undefined" << G4endl;
  }
  else if(UserDistType == "theta")
  {
    Theta = 10.;
    while(Theta > MaxTheta || Theta < MinTheta)
    {
      Theta = GenerateUserDefTheta();
    }
    Phi = 10.;
    while(Phi > MaxPhi || Phi < MinPhi)
    {
      rndm = angRndm->GenRandPhi();
      Phi = twopi * rndm;
    }
  }
  else if(UserDistType == "phi")
  {
    Theta = 10.;
    while(Theta > MaxTheta || Theta < MinTheta)
    {
      rndm = angRndm->GenRandTheta();
      Theta = std::acos(1. - (2. * rndm));
    }
    Phi = 10.;
    while(Phi > MaxPhi || Phi < MinPhi)
    {
      Phi = GenerateUserDefPhi();
    }
  }
  else if(UserDistType == "both")
  {
    Theta = 10.;
    while(Theta > MaxTheta || Theta < MinTheta)
    {
      Theta = GenerateUserDefTheta();
    }
    Phi = 10.;
    while(Phi > MaxPhi || Phi < MinPhi)
    {
      Phi = GenerateUserDefPhi();
    }
  }
  px = -std::sin(Theta) * std::cos(Phi);
  py = -std::sin(Theta) * std::sin(Phi);
  pz = -std::cos(Theta);

  pmag = std::sqrt((px*px) + (py*py) + (pz*pz));

  if(!UserWRTSurface)
  {
    G4double finx, finy, finz;      
    if (UserAngRef)
    {
      // Apply Rotation Matrix
      // x * AngRef1, y * AngRef2 and z * AngRef3
      finx = (px * AngRef1.x()) + (py * AngRef2.x()) + (pz * AngRef3.x());
      finy = (px * AngRef1.y()) + (py * AngRef2.y()) + (pz * AngRef3.y());
      finz = (px * AngRef1.z()) + (py * AngRef2.z()) + (pz * AngRef3.z());
    }
    else    // use mother coordinates
    {
      finx = px;
      finy = py;
      finz = pz;
    }
    G4double ResMag = std::sqrt((finx*finx) + (finy*finy) + (finz*finz));
    finx = finx/ResMag;
    finy = finy/ResMag;
    finz = finz/ResMag;
    
    mom.setX(finx);
    mom.setY(finy);
    mom.setZ(finz);
  } 
  else    // UserWRTSurface = true
  {
    G4double pxh = px/pmag;
    G4double pyh = py/pmag;
    G4double pzh = pz/pmag;
    if(verbosityLevel > 1)
    {
      G4cout << "SideRefVecs " << posDist->GetSideRefVec1()
             << posDist->GetSideRefVec2() << posDist->GetSideRefVec3()
             << G4endl;
      G4cout << "Raw Unit vector " << pxh
             << "," << pyh << "," << pzh << G4endl;
    }
    G4double resultx = (pxh*posDist->GetSideRefVec1().x())
                     + (pyh*posDist->GetSideRefVec2().x())
                     + (pzh*posDist->GetSideRefVec3().x());
    
    G4double resulty = (pxh*posDist->GetSideRefVec1().y())
                     + (pyh*posDist->GetSideRefVec2().y())
                     + (pzh*posDist->GetSideRefVec3().y());
    
    G4double resultz = (pxh*posDist->GetSideRefVec1().z())
                     + (pyh*posDist->GetSideRefVec2().z())
                     + (pzh*posDist->GetSideRefVec3().z());
    
    G4double ResMag = std::sqrt((resultx*resultx)
                              + (resulty*resulty)
                              + (resultz*resultz));
    resultx = resultx/ResMag;
    resulty = resulty/ResMag;
    resultz = resultz/ResMag;
    
    mom.setX(resultx);
    mom.setY(resulty);
    mom.setZ(resultz);
  }
  
  // particle_momentum_direction now contains unit momentum vector.

  if(verbosityLevel > 0 )
  {
    G4cout << "Final User Defined momentum vector "
           << particle_momentum_direction << G4endl;
  }
}

G4double G4SPSAngDistribution::GenerateUserDefTheta()
{
  // Create cumulative histogram if not already done so.
  // Then use RandFlat::shoot to generate the output Theta value.

  if(UserDistType == "NULL" || UserDistType == "phi")
  {
    // No user defined theta distribution
    G4cout << "Error ***********************" << G4endl;
    G4cout << "UserDistType = " << UserDistType << G4endl;
    return (0.);
  }
  
  // UserDistType = theta or both and so a theta distribution
  // is defined. This should be integrated if not already done.
  G4AutoLock l(&mutex);
  if(!IPDFThetaExist)
  {
    // IPDF has not been created, so create it
    //
    G4double bins[1024],vals[1024], sum;
    G4int ii;
    G4int maxbin = G4int(UDefThetaH.GetVectorLength());
    bins[0] = UDefThetaH.GetLowEdgeEnergy(std::size_t(0));
    vals[0] = UDefThetaH(std::size_t(0));
    sum = vals[0];
    for(ii=1; ii<maxbin; ++ii)
    {
      bins[ii] = UDefThetaH.GetLowEdgeEnergy(std::size_t(ii));
      vals[ii] = UDefThetaH(std::size_t(ii)) + vals[ii-1];
      sum = sum + UDefThetaH(std::size_t(ii));
    }
    for(ii=0; ii<maxbin; ++ii)
    {
      vals[ii] = vals[ii]/sum;
      IPDFThetaH.InsertValues(bins[ii], vals[ii]);
    }
      IPDFThetaExist = true;
  }
  l.unlock();

  // IPDF has been created so carry on
  //
  G4double rndm = G4UniformRand();
  return IPDFThetaH.GetEnergy(rndm);
}

G4double G4SPSAngDistribution::GenerateUserDefPhi()
{
  // Create cumulative histogram if not already done so.
  // Then use RandFlat::shoot to generate the output Theta value.

  if(UserDistType == "NULL" || UserDistType == "theta")
  {
    // No user defined phi distribution
    G4cout << "Error ***********************" << G4endl;
    G4cout << "UserDistType = " << UserDistType << G4endl;
    return(0.);
  }
  
  // UserDistType = phi or both and so a phi distribution
  // is defined. This should be integrated if not already done.
  G4AutoLock l(&mutex);
  if(!IPDFPhiExist)
  {
    // IPDF has not been created, so create it
    //
    G4double bins[1024],vals[1024], sum;
    G4int ii;
    G4int maxbin = G4int(UDefPhiH.GetVectorLength());
    bins[0] = UDefPhiH.GetLowEdgeEnergy(std::size_t(0));
    vals[0] = UDefPhiH(std::size_t(0));
    sum = vals[0];
    for(ii=1; ii<maxbin; ++ii)
    {
      bins[ii] = UDefPhiH.GetLowEdgeEnergy(std::size_t(ii));
      vals[ii] = UDefPhiH(std::size_t(ii)) + vals[ii-1];
      sum = sum + UDefPhiH(std::size_t(ii));
    }
    for(ii=0; ii<maxbin; ++ii)
    {
      vals[ii] = vals[ii]/sum;
      IPDFPhiH.InsertValues(bins[ii], vals[ii]);
    }
    IPDFPhiExist = true;
  }
  l.unlock();

  // IPDF has been create so carry on
  //
  G4double rndm = G4UniformRand();
  return IPDFPhiH.GetEnergy(rndm); 
}

void G4SPSAngDistribution::ReSetHist(const G4String& atype)
{
  G4AutoLock l(&mutex);
  if (atype == "theta")
  {
    UDefThetaH = IPDFThetaH = ZeroPhysVector ;
    IPDFThetaExist = false ;
  }
  else if (atype == "phi")
  {    
    UDefPhiH = IPDFPhiH = ZeroPhysVector ;
    IPDFPhiExist = false ;
  } 
  else
  {
    G4cout << "Error, histtype not accepted " << G4endl;
  }
}

G4ParticleMomentum G4SPSAngDistribution::GenerateOne()
{
  // Local copy for thread safety
  //
  G4ParticleMomentum localM = particle_momentum_direction;

  // Angular stuff
  //
  if(AngDistType == "iso")
    GenerateIsotropicFlux(localM);
  else if(AngDistType == "cos")
    GenerateCosineLawFlux(localM);
  else if(AngDistType == "planar")
    GeneratePlanarFlux(localM);
  else if(AngDistType == "beam1d" || AngDistType == "beam2d" )
    GenerateBeamFlux(localM);
  else if(AngDistType == "user")
    GenerateUserDefFlux(localM);
  else if(AngDistType == "focused")
    GenerateFocusedFlux(localM);
  else
    G4cout << "Error: AngDistType has unusual value" << G4endl;
  return localM;
}
