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
// G4SPSPosDistribution class implementation
//
// Author: Fan Lei, QinetiQ ltd. - 05/02/2004
// Customer: ESA/ESTEC
// Revisions: Andrea Dotti, John Allison, Makoto Asai, Maxime Chauvin
// --------------------------------------------------------------------

#include "G4SPSPosDistribution.hh"

#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4AutoLock.hh"
#include "G4AutoDelete.hh"

G4SPSPosDistribution::thread_data_t::thread_data_t()
{
  CSideRefVec1 = G4ThreeVector(CLHEP::HepXHat);
  CSideRefVec2 = G4ThreeVector(CLHEP::HepYHat);
  CSideRefVec3 = G4ThreeVector(CLHEP::HepZHat);
  CParticlePos = G4ThreeVector(0,0,0);
}

G4SPSPosDistribution::G4SPSPosDistribution()
{
  SourcePosType = "Point";
  Shape = "NULL";
  CentreCoords = G4ThreeVector(0,0,0);
  Rotx = CLHEP::HepXHat;
  Roty = CLHEP::HepYHat;
  Rotz = CLHEP::HepZHat;
  halfx = 0.;
  halfy = 0.;
  halfz = 0.;
  Radius = 0.;
  Radius0 = 0.;
  SR = 0.;
  SX = 0.;
  SY = 0.;
  ParAlpha = 0.;
  ParTheta = 0.;
  ParPhi = 0.;
  VolName = "NULL";
  verbosityLevel = 0 ;
  G4MUTEXINIT(a_mutex);
}

G4SPSPosDistribution::~G4SPSPosDistribution()
{
  G4MUTEXDESTROY(a_mutex);
}

void G4SPSPosDistribution::SetPosDisType(const G4String& PosType)
{
  SourcePosType = PosType;
}

void G4SPSPosDistribution::SetPosDisShape(const G4String& shapeType)
{
  Shape = shapeType;
}

void G4SPSPosDistribution::SetCentreCoords(const G4ThreeVector& coordsOfCentre)
{
  CentreCoords = coordsOfCentre;
}

void G4SPSPosDistribution::SetPosRot1(const G4ThreeVector& posrot1)
{
  // This should be x'

  Rotx = posrot1;
  if(verbosityLevel == 2)
  {
    G4cout << "Vector x' " << Rotx << G4endl;
  }
  GenerateRotationMatrices();
}

void G4SPSPosDistribution::SetPosRot2(const G4ThreeVector& posrot2)
{
  // This is a vector in the plane x'y' but need not be y'

  Roty = posrot2;
  if(verbosityLevel == 2)
  {
    G4cout << "The vector in the x'-y' plane " << Roty << G4endl;
  }
  GenerateRotationMatrices();
}

void G4SPSPosDistribution::SetHalfX(G4double xhalf)
{
  halfx = xhalf;
}

void G4SPSPosDistribution::SetHalfY(G4double yhalf)
{
  halfy = yhalf;
}

void G4SPSPosDistribution::SetHalfZ(G4double zhalf)
{
  halfz = zhalf;
}

void G4SPSPosDistribution::SetRadius(G4double rds)
{
  Radius = rds;
}

void G4SPSPosDistribution::SetRadius0(G4double rds)
{
  Radius0 = rds;
}

void G4SPSPosDistribution::SetBeamSigmaInR(G4double r)
{
  SX = SY = r;
  SR = r;
}

void G4SPSPosDistribution::SetBeamSigmaInX(G4double r)
{
  SX = r;
}

void G4SPSPosDistribution::SetBeamSigmaInY(G4double r)
{
  SY = r;
}

void G4SPSPosDistribution::SetParAlpha(G4double paralp)
{
  ParAlpha = paralp;
}

void G4SPSPosDistribution::SetParTheta(G4double parthe)
{
  ParTheta = parthe;
}

void G4SPSPosDistribution::SetParPhi(G4double parphi)
{
  ParPhi = parphi;
}

const G4String& G4SPSPosDistribution::GetPosDisType() const
{
  return SourcePosType;
}

const G4String& G4SPSPosDistribution::GetPosDisShape() const
{
  return Shape;
}

const G4ThreeVector& G4SPSPosDistribution::GetCentreCoords() const
{
  return CentreCoords;
}

G4double G4SPSPosDistribution::GetHalfX() const
{
  return halfx;
}

G4double G4SPSPosDistribution::GetHalfY() const
{
  return halfy;
}

G4double G4SPSPosDistribution::GetHalfZ() const
{
  return halfz;
}

G4double G4SPSPosDistribution::GetRadius() const
{
  return Radius;
}

void G4SPSPosDistribution::SetBiasRndm (G4SPSRandomGenerator* a)
{
  G4AutoLock l(&a_mutex);
  PosRndm = a;
}

void G4SPSPosDistribution::SetVerbosity(G4int a)
{
  verbosityLevel = a;
}

const G4String& G4SPSPosDistribution::GetSourcePosType() const
{
  return SourcePosType;
}

const G4ThreeVector& G4SPSPosDistribution::GetParticlePos() const
{
  return ThreadData.Get().CParticlePos;
}

const G4ThreeVector& G4SPSPosDistribution::GetSideRefVec1() const
{
  return ThreadData.Get().CSideRefVec1;
}

const G4ThreeVector& G4SPSPosDistribution::GetSideRefVec2() const
{
  return ThreadData.Get().CSideRefVec2;
}

const G4ThreeVector& G4SPSPosDistribution::GetSideRefVec3() const
{
  return ThreadData.Get().CSideRefVec3;
}

void G4SPSPosDistribution::GenerateRotationMatrices()
{
  // This takes in 2 vectors, x' and one in the plane x'-y',
  // and from these takes a cross product to calculate z'.
  // Then a cross product is taken between x' and z' to give y'

  Rotx = Rotx.unit(); // x'
  Roty = Roty.unit(); // vector in x'y' plane
  Rotz = Rotx.cross(Roty); // z'
  Rotz = Rotz.unit();
  Roty =Rotz.cross(Rotx); // y'
  Roty =Roty.unit();
  if(verbosityLevel == 2)
  {
    G4cout << "The new axes, x', y', z' "
           << Rotx << " " << Roty << " " << Rotz << G4endl;
  }
}

void G4SPSPosDistribution::ConfineSourceToVolume(const G4String& Vname)
{
  VolName = Vname;
  if(verbosityLevel == 2) { G4cout << VolName << G4endl; }

  if(VolName=="NULL")
  {
    if(verbosityLevel >= 1)
    { G4cout << "Volume confinement is set off." << G4endl; }
    Confine = false;
    return;
  }

  G4VPhysicalVolume* tempPV = nullptr;
  G4PhysicalVolumeStore* PVStore = G4PhysicalVolumeStore::GetInstance();
  if(verbosityLevel == 2) { G4cout << PVStore->size() << G4endl; }

  tempPV = PVStore->GetVolume(VolName);

  // the volume exists else it doesn't
  //
  if (tempPV != nullptr)
  {
    if(verbosityLevel >= 1)
    {
      G4cout << "Volume " << VolName << " exists" << G4endl;
    }
    Confine = true;
  }
  else
  {
    G4cout << " **** Error: Volume <" << VolName
           << "> does not exist **** " << G4endl;
    G4cout << " Ignoring confine condition" << G4endl;
    Confine = false;
    VolName = "NULL";
  }
}

void G4SPSPosDistribution::GeneratePointSource(G4ThreeVector& pos)
{
  // Generates Points given the point source

  if(SourcePosType == "Point")
  {
    pos = CentreCoords;
  }
  else
  {
    if(verbosityLevel >= 1)
    {
      G4cerr << "Error SourcePosType is not set to Point" << G4endl;
    }
  }
}

void G4SPSPosDistribution::GeneratePointsInBeam(G4ThreeVector& pos)
{
  G4double x, y, z;

  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  z = 0.;
  
  // Private Method to create points in a plane
  //
  if(Shape == "Circle")
  {
    x = Radius + 100.;
    y = Radius + 100.;
    while(std::sqrt((x*x) + (y*y)) > Radius)
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();

      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
    }
    x += G4RandGauss::shoot(0.0,SX) ;
    y += G4RandGauss::shoot(0.0,SY) ;
  }  
  else
  {
    // All other cases default to Rectangle case
    //
    x = PosRndm->GenRandX();
    y = PosRndm->GenRandY();
    x = (x*2.*halfx) - halfx;
    y = (y*2.*halfy) - halfy;
    x += G4RandGauss::shoot(0.0,SX);
    y += G4RandGauss::shoot(0.0,SY);
  }

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  //
  if(verbosityLevel >= 2)
  {
    G4cout << "Raw position " << x << "," << y << "," << z << G4endl;
  }
  tempx = (x * Rotx.x()) + (y * Roty.x()) + (z * Rotz.x());
  tempy = (x * Rotx.y()) + (y * Roty.y()) + (z * Rotz.y());
  tempz = (x * Rotx.z()) + (y * Roty.z()) + (z * Rotz.z());
  
  RandPos.setX(tempx);
  RandPos.setY(tempy);
  RandPos.setZ(tempz);
  
  // Translate
  //
  pos = CentreCoords + RandPos;
  if(verbosityLevel >= 1)
  {
    if(verbosityLevel >= 2)
    {
      G4cout << "Rotated Position " << RandPos << G4endl;
    }
    G4cout << "Rotated and Translated position " << pos << G4endl;
  }
}

void G4SPSPosDistribution::GeneratePointsInPlane(G4ThreeVector& pos)
{
  G4double x, y, z;
  G4double expression;
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  x = y = z = 0.;
  thread_data_t& td = ThreadData.Get();

  if(SourcePosType != "Plane" && verbosityLevel >= 1)
  {
    G4cerr << "Error: SourcePosType is not Plane" << G4endl;
  }

  // Private Method to create points in a plane
  //
  if(Shape == "Circle")
  {
    x = Radius + 100.;
    y = Radius + 100.;
    while(std::sqrt((x*x) + (y*y)) > Radius)
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();

      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
    }
  }
  else if(Shape == "Annulus")
  {
    x = Radius + 100.;
    y = Radius + 100.;
    while(std::sqrt((x*x) + (y*y)) > Radius
       || std::sqrt((x*x) + (y*y)) < Radius0 )
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();

      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
    }
  }
  else if(Shape == "Ellipse")
  {
    expression = 20.;
    while(expression > 1.)
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();

      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;

      expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
    }
  }
  else if(Shape == "Square")
  {
    x = PosRndm->GenRandX();
    y = PosRndm->GenRandY();
    x = (x*2.*halfx) - halfx;
    y = (y*2.*halfy) - halfy;
  }
  else if(Shape == "Rectangle")
  {
    x = PosRndm->GenRandX();
    y = PosRndm->GenRandY();
    x = (x*2.*halfx) - halfx;
    y = (y*2.*halfy) - halfy;
  }
  else
  {
    G4cout << "Shape not one of the plane types" << G4endl;
  }

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  //
  if(verbosityLevel == 2)
  {
    G4cout << "Raw position " << x << "," << y << "," << z << G4endl;
  }
  tempx = (x * Rotx.x()) + (y * Roty.x()) + (z * Rotz.x());
  tempy = (x * Rotx.y()) + (y * Roty.y()) + (z * Rotz.y());
  tempz = (x * Rotx.z()) + (y * Roty.z()) + (z * Rotz.z());

  RandPos.setX(tempx);
  RandPos.setY(tempy);
  RandPos.setZ(tempz);

  // Translate
  //
  pos = CentreCoords + RandPos;
  if(verbosityLevel >= 1)
  {
    if(verbosityLevel == 2)
    {
      G4cout << "Rotated Position " << RandPos << G4endl;
    }
    G4cout << "Rotated and Translated position " << pos << G4endl;
  }

  // For Cosine-Law make SideRefVecs = to Rotation matrix vectors
  // Note: these need to be thread-local, use G4Cache
  //
  td.CSideRefVec1 = Rotx;
  td.CSideRefVec2 = Roty;
  td.CSideRefVec3 = Rotz;

  // If rotation matrix z' point to origin then invert the matrix
  // So that SideRefVecs point away
  //
  if((CentreCoords.x() > 0. && Rotz.x() < 0.)
     || (CentreCoords.x() < 0. && Rotz.x() > 0.)
     || (CentreCoords.y() > 0. && Rotz.y() < 0.)
     || (CentreCoords.y() < 0. && Rotz.y() > 0.)
     || (CentreCoords.z() > 0. && Rotz.z() < 0.)
     || (CentreCoords.z() < 0. && Rotz.z() > 0.))
  {
    // Invert y and z
    //
    td.CSideRefVec2 = - td.CSideRefVec2;
    td.CSideRefVec3 = - td.CSideRefVec3;
  }
  if(verbosityLevel == 2)
  {
    G4cout << "Reference vectors for cosine-law "
           << td.CSideRefVec1<< " " << td.CSideRefVec2
           << " " << td.CSideRefVec3 << G4endl;
  }
}

void G4SPSPosDistribution::GeneratePointsOnSurface(G4ThreeVector& pos)
{
  // Private method to create points on a surface
  //
  G4double theta, phi;
  G4double x, y, z;
  x = y = z = 0.;
  G4double expression;
  G4ThreeVector RandPos;

  thread_data_t& td = ThreadData.Get();
  if(SourcePosType != "Surface" && verbosityLevel >= 1)
  {
    G4cout << "Error SourcePosType not Surface" << G4endl;
  }

  if(Shape == "Sphere")
  {
    theta = PosRndm->GenRandPosTheta();
    phi = PosRndm->GenRandPosPhi();
    theta = std::acos(1. - 2.*theta); // theta isotropic
    phi = phi * 2. * pi;
      
    x = Radius * std::sin(theta) * std::cos(phi);
    y = Radius * std::sin(theta) * std::sin(phi);
    z = Radius * std::cos(theta);
      
    RandPos.setX(x);
    RandPos.setY(y);
    RandPos.setZ(z);

    // Cosine-law (not a good idea to use this here)
    //
    G4ThreeVector zdash(x,y,z);
    zdash = zdash.unit();
    G4ThreeVector xdash = Rotz.cross(zdash);
    G4ThreeVector ydash = xdash.cross(zdash);
    td.CSideRefVec1 = xdash.unit();
    td.CSideRefVec2 = ydash.unit();
    td.CSideRefVec3 = zdash.unit();
  }
  else if(Shape == "Ellipsoid")
  {
    G4double minphi, maxphi, middlephi;
    G4double answer, constant;

    constant = pi/(halfx*halfx) + pi/(halfy*halfy) + twopi/(halfz*halfz);
      
    // Simplified approach
    //
    theta = PosRndm->GenRandPosTheta();
    phi = PosRndm->GenRandPosPhi();
      
    theta = std::acos(1. - 2.*theta);
    minphi = 0.;
    maxphi = twopi;
    while(maxphi-minphi > 0.)
    {
      middlephi = (maxphi+minphi)/2.;
      answer = (1./(halfx*halfx))*(middlephi/2. + std::sin(2*middlephi)/4.)
             + (1./(halfy*halfy))*(middlephi/2. - std::sin(2*middlephi)/4.)
             + middlephi/(halfz*halfz);
      answer = answer/constant;
      if(answer > phi) maxphi = middlephi;
      if(answer < phi) minphi = middlephi;
      if(std::fabs(answer-phi) <= 0.00001)
      {
        minphi = maxphi +1;
        phi = middlephi;
      }
    }

    // x,y and z form a unit vector. Put this onto the ellipse.
    //
    x = std::sin(theta)*std::cos(phi);
    y = std::sin(theta)*std::sin(phi);
    z = std::cos(theta);

    G4double lhs;

    // Solve for x
    //
    G4double numYinX = y/x;
    G4double numZinX = z/x;
    G4double tempxvar;          
    tempxvar = 1./(halfx*halfx)+(numYinX*numYinX)/(halfy*halfy)
             + (numZinX*numZinX)/(halfz*halfz);
    tempxvar = 1./tempxvar;
    G4double coordx = std::sqrt(tempxvar);
  
    // Solve for y
    //
    G4double numXinY = x/y;
    G4double numZinY = z/y;
    G4double tempyvar;
    tempyvar = (numXinY*numXinY)/(halfx*halfx)+1./(halfy*halfy)
             + (numZinY*numZinY)/(halfz*halfz);
    tempyvar = 1./tempyvar;
    G4double coordy = std::sqrt(tempyvar);
   
    // Solve for z
    //
    G4double numXinZ = x/z;
    G4double numYinZ = y/z;
    G4double tempzvar;
    tempzvar = (numXinZ*numXinZ)/(halfx*halfx)
             + (numYinZ*numYinZ)/(halfy*halfy)+1./(halfz*halfz);
    tempzvar = 1./tempzvar;
    G4double coordz = std::sqrt(tempzvar);

    lhs = std::sqrt((coordx*coordx)/(halfx*halfx) +
                    (coordy*coordy)/(halfy*halfy) +
                    (coordz*coordz)/(halfz*halfz));
      
    if(std::fabs(lhs-1.) > 0.001 && verbosityLevel >= 1)
    {
      G4cout << "Error: theta, phi not really on ellipsoid" << G4endl;
    }

    // coordx, coordy and coordz are all positive
    //
    G4double TestRandVar = G4UniformRand();
    if(TestRandVar > 0.5)
    {
      coordx = -coordx;
    }
    TestRandVar = G4UniformRand();
    if(TestRandVar > 0.5)
    {
      coordy = -coordy;
    }
    TestRandVar = G4UniformRand();
    if(TestRandVar > 0.5)
    {
      coordz = -coordz;
    }

    RandPos.setX(coordx);
    RandPos.setY(coordy);
    RandPos.setZ(coordz);

    // Cosine-law (not a good idea to use this here)
    //
    G4ThreeVector zdash(coordx,coordy,coordz);
    zdash = zdash.unit();
    G4ThreeVector xdash = Rotz.cross(zdash);
    G4ThreeVector ydash = xdash.cross(zdash);
    td.CSideRefVec1 = xdash.unit();
    td.CSideRefVec2 = ydash.unit();
    td.CSideRefVec3 = zdash.unit();
  }
  else if(Shape == "Cylinder")
  {
    G4double AreaTop, AreaBot, AreaLat;
    G4double AreaTotal, prob1, prob2, prob3;
    G4double testrand;

    // User giver Radius and z-half length
    // Calculate surface areas, maybe move this to 
    // a different method

    AreaTop = pi * Radius * Radius;
    AreaBot = AreaTop;
    AreaLat = 2. * pi * Radius * 2. * halfz;
    AreaTotal = AreaTop + AreaBot + AreaLat;
      
    prob1 = AreaTop / AreaTotal;
    prob2 = AreaBot / AreaTotal;
    prob3 = 1.00 - prob1 - prob2;
    if(std::fabs(prob3 - (AreaLat/AreaTotal)) >= 0.001)
    {
      if(verbosityLevel >= 1)
      {
        G4cout << AreaLat/AreaTotal << " " << prob3 << G4endl;
      }
      G4cout << "Error in prob3" << G4endl;
    }

    // Decide surface to calculate point on.

    testrand = G4UniformRand();
    if(testrand <= prob1)  // Point on Top surface
    {
      z = halfz;
      x = Radius + 100.;
      y = Radius + 100.;
      while(((x*x)+(y*y)) > (Radius*Radius))
      {
         x = PosRndm->GenRandX();
         y = PosRndm->GenRandY();

         x = x * 2. * Radius;
         y = y * 2. * Radius;
         x = x - Radius;
         y = y - Radius;
      }
      // Cosine law
      //
      td.CSideRefVec1 = Rotx;
      td.CSideRefVec2 = Roty;
      td.CSideRefVec3 = Rotz;
    }
    else if((testrand > prob1) && (testrand <= (prob1 + prob2)))
    {                          // Point on Bottom surface
      z = -halfz;
      x = Radius + 100.;
      y = Radius + 100.;
      while(((x*x)+(y*y)) > (Radius*Radius))
      {
        x = PosRndm->GenRandX();
        y = PosRndm->GenRandY();

        x = x * 2. * Radius;
        y = y * 2. * Radius;
        x = x - Radius;
        y = y - Radius;
      }
      // Cosine law
      //
      td.CSideRefVec1 = Rotx;
      td.CSideRefVec2 = -Roty;
      td.CSideRefVec3 = -Rotz;
    }
    else if(testrand > (prob1+prob2))  // Point on Lateral Surface
    {
      G4double rand;

      rand = PosRndm->GenRandPosPhi();
      rand = rand * 2. * pi;

      x = Radius * std::cos(rand);
      y = Radius * std::sin(rand);

      z = PosRndm->GenRandZ();

      z = z * 2. *halfz;
      z = z - halfz;
          
      // Cosine law
      //
      G4ThreeVector zdash(x,y,0.);
      zdash = zdash.unit();
      G4ThreeVector xdash = Rotz.cross(zdash);
      G4ThreeVector ydash = xdash.cross(zdash);
      td.CSideRefVec1 = xdash.unit();
      td.CSideRefVec2 = ydash.unit();
      td.CSideRefVec3 = zdash.unit();
    }
    else
    {
      G4cout << "Error: testrand " << testrand << G4endl;
    }

    RandPos.setX(x);
    RandPos.setY(y);
    RandPos.setZ(z);
  }
  else if(Shape == "EllipticCylinder")
  {
    G4double AreaTop, AreaBot, AreaLat, AreaTotal;
    G4double h;
    G4double prob1, prob2, prob3;
    G4double testrand;

    // User giver x, y and z-half length
    // Calculate surface areas, maybe move this to
    // a different method

    AreaTop = pi * halfx * halfy;
    AreaBot = AreaTop;

    // Using circumference approximation by Ramanujan (order h^3)
    // AreaLat = 2*halfz * pi*( 3*(halfx + halfy)
    //         - std::sqrt((3*halfx + halfy) * (halfx + 3*halfy)) );
    // Using circumference approximation by Ramanujan (order h^5)
    //
    h = ((halfx - halfy)*(halfx - halfy)) / ((halfx + halfy)*(halfx + halfy));
    AreaLat = 2*halfz * pi*(halfx + halfy)
            * (1 + (3*h)/(10 + std::sqrt(4 - 3*h)));
    AreaTotal = AreaTop + AreaBot + AreaLat;

    prob1 = AreaTop / AreaTotal;
    prob2 = AreaBot / AreaTotal;
    prob3 = 1.00 - prob1 - prob2;
    if(std::fabs(prob3 - (AreaLat/AreaTotal)) >= 0.001)
    {
      if(verbosityLevel >= 1)
      {
        G4cout << AreaLat/AreaTotal << " " << prob3 << G4endl;
      }
      G4cout << "Error in prob3" << G4endl;
    }

    // Decide surface to calculate point on

    testrand = G4UniformRand();
    if(testrand <= prob1)  // Point on Top surface
    {
      z = halfz;
      expression = 20.;
      while(expression > 1.)
      {
        x = PosRndm->GenRandX();
        y = PosRndm->GenRandY();

        x = (x * 2. * halfx) - halfx;
        y = (y * 2. * halfy) - halfy;

        expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
      }
    }
    else if((testrand > prob1) && (testrand <= (prob1 + prob2)))
    {                          // Point on Bottom surface
      z = -halfz;
      expression = 20.;
      while(expression > 1.)
      {
        x = PosRndm->GenRandX();
        y = PosRndm->GenRandY();

        x = (x * 2. * halfx) - halfx;
        y = (y * 2. * halfy) - halfy;

        expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
      }
    }
    else if(testrand > (prob1+prob2))  // Point on Lateral Surface
    {
      G4double rand;

      rand = PosRndm->GenRandPosPhi();
      rand = rand * 2. * pi;

      x = halfx * std::cos(rand);
      y = halfy * std::sin(rand);

      z = PosRndm->GenRandZ();

      z = (z * 2. * halfz) - halfz;
    }
    else
    {
      G4cout << "Error: testrand " << testrand << G4endl;
    }

    RandPos.setX(x);
    RandPos.setY(y);
    RandPos.setZ(z);

    // Cosine law
    //
    G4ThreeVector zdash(x,y,z);
    zdash = zdash.unit();
    G4ThreeVector xdash = Rotz.cross(zdash);
    G4ThreeVector ydash = xdash.cross(zdash);
    td.CSideRefVec1 = xdash.unit();
    td.CSideRefVec2 = ydash.unit();
    td.CSideRefVec3 = zdash.unit();
  }
  else if(Shape == "Para")
  {
    // Right Parallelepiped.
    // User gives x,y,z half lengths and ParAlpha
    // ParTheta and ParPhi
    // +x = <1, -x >1 & <2, +y >2 & <3, -y >3 &<4
    // +z >4 & < 5, -z >5 &<6

    G4double testrand = G4UniformRand();
    G4double AreaX = halfy * halfz * 4.;
    G4double AreaY = halfx * halfz * 4.;
    G4double AreaZ = halfx * halfy * 4.;
    G4double AreaTotal = 2*(AreaX + AreaY + AreaZ);
    G4double Probs[6];
    Probs[0] = AreaX/AreaTotal;
    Probs[1] = Probs[0] + AreaX/AreaTotal;
    Probs[2] = Probs[1] + AreaY/AreaTotal;
    Probs[3] = Probs[2] + AreaY/AreaTotal;
    Probs[4] = Probs[3] + AreaZ/AreaTotal;
    Probs[5] = Probs[4] + AreaZ/AreaTotal;
      
    x = PosRndm->GenRandX();
    y = PosRndm->GenRandY();
    z = PosRndm->GenRandZ();
      
    x = x * halfx * 2.;
    x = x - halfx;
    y = y * halfy * 2.;
    y = y - halfy;
    z = z * halfz * 2.;
    z = z - halfz;

    // Pick a side first
    //
    if(testrand < Probs[0])
    {
      // side is +x

      x = halfx + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
      y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
      // z = z;

      // Cosine-law
      //
      G4ThreeVector xdash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
                          halfz*std::tan(ParTheta)*std::sin(ParPhi),
                          halfz/std::cos(ParPhi));
      G4ThreeVector ydash(halfy*std::tan(ParAlpha), -halfy, 0.0);
      xdash = xdash.unit();
      ydash = ydash.unit();
      G4ThreeVector zdash = xdash.cross(ydash);
      td.CSideRefVec1 = xdash.unit();
      td.CSideRefVec2 = ydash.unit();
      td.CSideRefVec3 = zdash.unit();
    }
    else if(testrand >= Probs[0] && testrand < Probs[1])
    {
      // side is -x

      x = -halfx + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
      y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
      // z = z;

      // Cosine-law
      //
      G4ThreeVector xdash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
                          halfz*std::tan(ParTheta)*std::sin(ParPhi),
                          halfz/std::cos(ParPhi));
      G4ThreeVector ydash(halfy*std::tan(ParAlpha), halfy, 0.0);
      xdash = xdash.unit();
      ydash = ydash.unit();
      G4ThreeVector zdash = xdash.cross(ydash);
      td.CSideRefVec1 = xdash.unit();
      td.CSideRefVec2 = ydash.unit();
      td.CSideRefVec3 = zdash.unit();
    }
    else if(testrand >= Probs[1] && testrand < Probs[2])
    {
      // side is +y

      x = x + z*std::tan(ParTheta)*std::cos(ParPhi) + halfy*std::tan(ParAlpha);
      y = halfy + z*std::tan(ParTheta)*std::sin(ParPhi);
      // z = z;

      // Cosine-law
      //
      G4ThreeVector ydash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
                          halfz*std::tan(ParTheta)*std::sin(ParPhi),
                          halfz/std::cos(ParPhi));
      ydash = ydash.unit();
      G4ThreeVector xdash = Roty.cross(ydash);
      G4ThreeVector zdash = xdash.cross(ydash);
      td.CSideRefVec1 = xdash.unit();
      td.CSideRefVec2 = -ydash.unit();
      td.CSideRefVec3 = -zdash.unit();
    }
    else if(testrand >= Probs[2] && testrand < Probs[3])
    {
      // side is -y

      x = x + z*std::tan(ParTheta)*std::cos(ParPhi) - halfy*std::tan(ParAlpha);
      y = -halfy + z*std::tan(ParTheta)*std::sin(ParPhi);
      // z = z;

      // Cosine-law
      //
      G4ThreeVector ydash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
                          halfz*std::tan(ParTheta)*std::sin(ParPhi),
                          halfz/std::cos(ParPhi));
      ydash = ydash.unit();
      G4ThreeVector xdash = Roty.cross(ydash);
      G4ThreeVector zdash = xdash.cross(ydash);
      td.CSideRefVec1 = xdash.unit();
      td.CSideRefVec2 = ydash.unit();
      td.CSideRefVec3 = zdash.unit();
    }
    else if(testrand >= Probs[3] && testrand < Probs[4])
    {
      // side is +z

      z = halfz;
      y = y + halfz*std::sin(ParPhi)*std::tan(ParTheta);
      x = x + halfz*std::cos(ParPhi)*std::tan(ParTheta) + y*std::tan(ParAlpha);

      // Cosine-law
      //
      td.CSideRefVec1 = Rotx;
      td.CSideRefVec2 = Roty;
      td.CSideRefVec3 = Rotz;
    }
    else if(testrand >= Probs[4] && testrand < Probs[5])
    {
      // side is -z

      z = -halfz;
      y = y - halfz*std::sin(ParPhi)*std::tan(ParTheta);
      x = x - halfz*std::cos(ParPhi)*std::tan(ParTheta) + y*std::tan(ParAlpha);

      // Cosine-law
      //
      td.CSideRefVec1 = Rotx;
      td.CSideRefVec2 = -Roty;
      td.CSideRefVec3 = -Rotz;
    }
    else
    {
      G4cout << "Error: testrand out of range" << G4endl;
      if(verbosityLevel >= 1)
      {
        G4cout << "testrand=" << testrand << " Probs[5]=" << Probs[5] <<G4endl;
      }
    }

    RandPos.setX(x);
    RandPos.setY(y);
    RandPos.setZ(z);
  }

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  //
  if(verbosityLevel == 2)
  {
    G4cout << "Raw position " << RandPos << G4endl;
  }
  x=(RandPos.x()*Rotx.x())+(RandPos.y()*Roty.x())+(RandPos.z()*Rotz.x());
  y=(RandPos.x()*Rotx.y())+(RandPos.y()*Roty.y())+(RandPos.z()*Rotz.y());
  z=(RandPos.x()*Rotx.z())+(RandPos.y()*Roty.z())+(RandPos.z()*Rotz.z());
  
  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);

  // Translate
  //
  pos = CentreCoords + RandPos;

  if(verbosityLevel >= 1)
  {
    if(verbosityLevel == 2)
    {
      G4cout << "Rotated position " << RandPos << G4endl;
    }
    G4cout << "Rotated and translated position " << pos << G4endl;
  }
  if(verbosityLevel == 2)
  {
    G4cout << "Reference vectors for cosine-law " << td.CSideRefVec1
           << " " << td.CSideRefVec2 << " " << td.CSideRefVec3 << G4endl;
  }
}

void G4SPSPosDistribution::GeneratePointsInVolume(G4ThreeVector& pos)
{
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  G4double x, y, z;
  G4double expression;
  x = y = z = 0.;

  if(SourcePosType != "Volume" && verbosityLevel >= 1)
  {
    G4cout << "Error SourcePosType not Volume" << G4endl;
  }

  // Private method to create points in a volume
  //
  if(Shape == "Sphere")
  {
    x = Radius*2.;
    y = Radius*2.;
    z = Radius*2.;
    while(((x*x)+(y*y)+(z*z)) > (Radius*Radius))
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();
      z = PosRndm->GenRandZ();

      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*Radius) - Radius;
    }
  }
  else if(Shape == "Ellipsoid")
  {
    G4double temp;
    temp = 100.;
    while(temp > 1.)
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();
      z = PosRndm->GenRandZ();

      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
      z = (z*2.*halfz) - halfz;
          
      temp = ((x*x)/(halfx*halfx))+((y*y)/(halfy*halfy))+((z*z)/(halfz*halfz));
    }
  }
  else if(Shape == "Cylinder")
  {
    x = Radius*2.;
    y = Radius*2.;
    while(((x*x)+(y*y)) > (Radius*Radius))
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();
      z = PosRndm->GenRandZ();

      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*halfz) - halfz;
    }
  }
  else if(Shape == "EllipticCylinder")
  {
    expression = 20.;
    while(expression > 1.)
    {
      x = PosRndm->GenRandX();
      y = PosRndm->GenRandY();
      z = PosRndm->GenRandZ();

      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
      z = (z*2.*halfz) - halfz;

      expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
    }
  }
  else if(Shape == "Para")
  {
    x = PosRndm->GenRandX();
    y = PosRndm->GenRandY();
    z = PosRndm->GenRandZ();
    x = (x*2.*halfx) - halfx;
    y = (y*2.*halfy) - halfy;
    z = (z*2.*halfz) - halfz;
    x = x + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
    y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
    // z = z;
  }
  else
  {
    G4cout << "Error: Volume Shape does not exist" << G4endl;
  }

  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  //
  tempx = (x * Rotx.x()) + (y * Roty.x()) + (z * Rotz.x());
  tempy = (x * Rotx.y()) + (y * Roty.y()) + (z * Rotz.y());
  tempz = (x * Rotx.z()) + (y * Roty.z()) + (z * Rotz.z());

  RandPos.setX(tempx);
  RandPos.setY(tempy);
  RandPos.setZ(tempz);

  // Translate
  //
  pos = CentreCoords + RandPos;

  if(verbosityLevel == 2)
  {
    G4cout << "Raw position " << x << "," << y << "," << z << G4endl;
    G4cout << "Rotated position " << RandPos << G4endl;
  }
  if(verbosityLevel >= 1)
  {
    G4cout << "Rotated and translated position " << pos << G4endl;
  }

  // Cosine-law (not a good idea to use this here)
  //
  G4ThreeVector zdash(tempx,tempy,tempz);
  zdash = zdash.unit();
  G4ThreeVector xdash = Rotz.cross(zdash);
  G4ThreeVector ydash = xdash.cross(zdash);
  thread_data_t& td = ThreadData.Get();
  td.CSideRefVec1 = xdash.unit();
  td.CSideRefVec2 = ydash.unit();
  td.CSideRefVec3 = zdash.unit();
  if(verbosityLevel == 2)
  {
    G4cout << "Reference vectors for cosine-law " << td.CSideRefVec1
           << " " << td.CSideRefVec2 << " " << td.CSideRefVec3 << G4endl;
  } 
}

G4bool G4SPSPosDistribution::IsSourceConfined(G4ThreeVector& pos)
{
  // Method to check point is within the volume specified

  if(!Confine)
  {
    G4cout << "Error: Confine is false" << G4endl;
  }
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector* ptr = &null;

  // Check position is within VolName, if so true, else false
  //
  G4VPhysicalVolume* theVolume;
  G4Navigator* gNavigator = G4TransportationManager::GetTransportationManager()
                          ->GetNavigatorForTracking();
  theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
  if(theVolume == nullptr) { return false; }
  G4String theVolName = theVolume->GetName();

  if(theVolName == VolName)
  {
    if(verbosityLevel >= 1)
    {
      G4cout << "Particle is in volume " << VolName << G4endl;
    }
    return true;
  }
  
  return false;
 
}

G4ThreeVector G4SPSPosDistribution::GenerateOne()
{
  G4ThreeVector localP;
  G4bool srcconf = false;
  G4int LoopCount = 0;
  while(!srcconf)
  {
    if(SourcePosType == "Point")
      GeneratePointSource(localP);
    else if(SourcePosType == "Beam")
      GeneratePointsInBeam(localP);
    else if(SourcePosType == "Plane")
      GeneratePointsInPlane(localP);
    else if(SourcePosType == "Surface")
      GeneratePointsOnSurface(localP);
    else if(SourcePosType == "Volume")
      GeneratePointsInVolume(localP);
    else
    {
      G4ExceptionDescription msg;
      msg << "Error: SourcePosType undefined\n";
      msg << "Generating point source\n";
      G4Exception("G4SPSPosDistribution::GenerateOne()",
                  "G4GPS001", JustWarning, msg);
      GeneratePointSource(localP);
    }
    if(Confine)
    {
      srcconf = IsSourceConfined(localP);
      // if source in confined srcconf = true terminating the loop
      // if source isnt confined srcconf = false and loop continues
    }
    else if(!Confine)
    {
      srcconf = true; // terminate loop
    }
    ++LoopCount;
    if(LoopCount == 100000)
    {
      G4ExceptionDescription msg;
      msg << "LoopCount = 100000\n";
      msg << "Either the source distribution >> confinement\n";
      msg << "or any confining volume may not overlap with\n";
      msg << "the source distribution or any confining volumes\n";
      msg << "may not exist\n"<< G4endl;
      msg << "If you have set confine then this will be ignored\n";
      msg << "for this event.\n" << G4endl;
      G4Exception("G4SPSPosDistribution::GenerateOne()",
                  "G4GPS001", JustWarning, msg);
      srcconf = true; // Avoids an infinite loop
    }
  }
  ThreadData.Get().CParticlePos = localP;
  return localP;
}

