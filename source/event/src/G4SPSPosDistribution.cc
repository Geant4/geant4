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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4SPSPosDistribution.cc
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//    Based on the G4GeneralParticleSource class in Geant4 v6.0
//
///////////////////////////////////////////////////////////////////////////////
//
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4SPSPosDistribution.hh"

G4SPSPosDistribution::G4SPSPosDistribution()
{

  // Initialise all variables
  // Position distribution Variables

  SourcePosType = "Point";
  Shape = "NULL";
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
  CentreCoords = G4ThreeVector(0., 0., 0.);
  Rotx = CLHEP::HepXHat;
  Roty = CLHEP::HepYHat;
  Rotz = CLHEP::HepZHat;
  Confine = false; //If true confines source distribution to VolName
  VolName = "NULL";
  SideRefVec1 = CLHEP::HepXHat; // x-axis
  SideRefVec2 = CLHEP::HepYHat; // y-axis
  SideRefVec3 = CLHEP::HepZHat; // z-axis
  verbosityLevel = 0 ;
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}

G4SPSPosDistribution::~G4SPSPosDistribution()
{
}

void G4SPSPosDistribution::SetPosDisType(G4String PosType)
{
  SourcePosType = PosType;
}

void G4SPSPosDistribution::SetPosDisShape(G4String shapeType)
{
  Shape = shapeType;
}

void G4SPSPosDistribution::SetCentreCoords(G4ThreeVector coordsOfCentre)
{
  CentreCoords = coordsOfCentre;
}

void G4SPSPosDistribution::SetPosRot1(G4ThreeVector posrot1)
{
  // This should be x'
  Rotx = posrot1;
  if(verbosityLevel == 2)
    {
      G4cout << "Vector x' " << Rotx << G4endl;
    }
  GenerateRotationMatrices();
}

void G4SPSPosDistribution::SetPosRot2(G4ThreeVector posrot2)
{
  // This is a vector in the plane x'y' but need not
  // be y'
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

void G4SPSPosDistribution::SetRadius(G4double rad)
{
  Radius = rad;
}

void G4SPSPosDistribution::SetRadius0(G4double rad)
{
  Radius0 = rad;
}

void G4SPSPosDistribution::SetBeamSigmaInR(G4double r)
{
  SR = r;
  SX = SY = r/std::sqrt(2.);
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

void G4SPSPosDistribution::GenerateRotationMatrices()
{
  // This takes in 2 vectors, x' and one in the plane x'-y',
  // and from these takes a cross product to calculate z'.
  // Then a cross product is taken between x' and z' to give
  // y'.
  Rotx = Rotx.unit(); // x'
  Roty = Roty.unit(); // vector in x'y' plane
  Rotz = Rotx.cross(Roty); // z'
  Roty = Rotz.cross(Rotx); // y'
  if(verbosityLevel == 2)
    {
      G4cout << "The new axes, x', y', z' " << Rotx << " " << Roty << " " << Rotz << G4endl;
    }
}

void G4SPSPosDistribution::ConfineSourceToVolume(G4String Vname)
{
  VolName = Vname;
  if(verbosityLevel == 2)
    G4cout << VolName << G4endl;
  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = 0;
  G4String theRequiredVolumeName = VolName;
  PVStore      = G4PhysicalVolumeStore::GetInstance();
  G4int      i = 0;
  G4bool found = false;
  if(verbosityLevel == 2)
    G4cout << PVStore->size() << G4endl;
  while (!found && i<G4int(PVStore->size())) {
    tempPV = (*PVStore)[i];
    found  = tempPV->GetName() == theRequiredVolumeName;
    if(verbosityLevel == 2)
      G4cout << i << " " << " " << tempPV->GetName() << " " << theRequiredVolumeName << " " << found << G4endl;
    if (!found)
      {i++;}
  }
  // found = true then the volume exists else it doesnt.
  if(found == true)
    {
      if(verbosityLevel >= 1)
	G4cout << "Volume " << VolName << " exists" << G4endl;
      Confine = true;
    }
  else
    {
      G4cout << " **** Error: Volume does not exist **** " << G4endl;
      G4cout << " Ignoring confine condition" << G4endl;
      Confine = false;
      VolName = "NULL";
    }

}

void G4SPSPosDistribution::GeneratePointSource()
{
  // Generates Points given the point source.
  if(SourcePosType == "Point")
    particle_position = CentreCoords;
  else
    if(verbosityLevel >= 1)
      G4cout << "Error SourcePosType is not set to Point" << G4endl;
}

void G4SPSPosDistribution::GeneratePointsInBeam()
{
  G4double x, y, z;

  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  z = 0.;
  
  // Private Method to create points in a plane
  if(Shape == "Circle")
    {
      x = Radius + 100.;
      y = Radius + 100.;
      while(std::sqrt((x*x) + (y*y)) > Radius)
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	}
      x += G4RandGauss::shoot(0.0,SX) ;
      y += G4RandGauss::shoot(0.0,SY) ;
    }  
  else
    {
      // all other cases default to Rectangle case
      x = posRndm->GenRandX();
      y = posRndm->GenRandY();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
      x += G4RandGauss::shoot(0.0,SX);
      y += G4RandGauss::shoot(0.0,SY);
    } 
  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
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
  particle_position = CentreCoords + RandPos;
  if(verbosityLevel >= 1)
    {
      if(verbosityLevel >= 2)
	{
	  G4cout << "Rotated Position " << RandPos << G4endl;
	}
      G4cout << "Rotated and Translated position " << particle_position << G4endl;
    }
}

void G4SPSPosDistribution::GeneratePointsInPlane()
{
  G4double x, y, z;
  G4double expression;
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  x = y = z = 0.;

  if(SourcePosType != "Plane" && verbosityLevel >= 1)
    G4cout << "Error: SourcePosType is not Plane" << G4endl;

  // Private Method to create points in a plane
  if(Shape == "Circle")
    {
      x = Radius + 100.;
      y = Radius + 100.;
      while(std::sqrt((x*x) + (y*y)) > Radius)
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	}
    }
  else if(Shape == "Annulus")
    {
      x = Radius + 100.;
      y = Radius + 100.;
      while(std::sqrt((x*x) + (y*y)) > Radius || std::sqrt((x*x) + (y*y)) < Radius0 )
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	}
    }
  else if(Shape == "Ellipse")
    {
      expression = 20.;
      while(expression > 1.)
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();

	  x = (x*2.*halfx) - halfx;
	  y = (y*2.*halfy) - halfy;

	  expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
	}
    }
  else if(Shape == "Square")
    {
      x = posRndm->GenRandX();
      y = posRndm->GenRandY();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
    }
  else if(Shape == "Rectangle")
    {
      x = posRndm->GenRandX();
      y = posRndm->GenRandY();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
    }
  else
    G4cout << "Shape not one of the plane types" << G4endl;

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
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
  particle_position = CentreCoords + RandPos;
  if(verbosityLevel >= 1)
    {
      if(verbosityLevel == 2)
	{
	  G4cout << "Rotated Position " << RandPos << G4endl;
	}
      G4cout << "Rotated and Translated position " << particle_position << G4endl;
    }

  // For Cosine-Law make SideRefVecs = to Rotation matrix vectors
  SideRefVec1 = Rotx;
  SideRefVec2 = Roty;
  SideRefVec3 = Rotz;
  // If rotation matrix z' point to origin then invert the matrix
  // So that SideRefVecs point away.
  if((CentreCoords.x() > 0. && Rotz.x() < 0.)
     || (CentreCoords.x() < 0. && Rotz.x() > 0.)
     || (CentreCoords.y() > 0. && Rotz.y() < 0.)
     || (CentreCoords.y() < 0. && Rotz.y() > 0.)
     || (CentreCoords.z() > 0. && Rotz.z() < 0.)
     || (CentreCoords.z() < 0. && Rotz.z() > 0.))
    {
      // Invert y and z.
      SideRefVec2 = -SideRefVec2;
      SideRefVec3 = -SideRefVec3;
    }
  if(verbosityLevel == 2)
    {
      G4cout << "Reference vectors for cosine-law " << SideRefVec1 << " " << SideRefVec2 << " " << SideRefVec3 << G4endl;
    }
}

void G4SPSPosDistribution::GeneratePointsOnSurface()
{
  //Private method to create points on a surface
  G4double theta, phi;
  G4double x, y, z;
  x = y = z = 0.;
  G4ThreeVector RandPos;
  //  G4double tempx, tempy, tempz;

  if(SourcePosType != "Surface" && verbosityLevel >= 1)
    G4cout << "Error SourcePosType not Surface" << G4endl;

  if(Shape == "Sphere")
    {
      G4double tantheta;
      theta = posRndm->GenRandPosTheta();
      phi = posRndm->GenRandPosPhi();
      theta = std::acos(1. - 2.*theta); // theta isotropic
      phi = phi * 2. * pi;
      tantheta = std::tan(theta);
      
      x = Radius * std::sin(theta) * std::cos(phi);
      y = Radius * std::sin(theta) * std::sin(phi);
      z = Radius * std::cos(theta);
      
      RandPos.setX(x);
      RandPos.setY(y);
      RandPos.setZ(z);

      // Cosine-law (not a good idea to use this here)
      G4ThreeVector zdash(x,y,z);
      zdash = zdash.unit();
      G4ThreeVector xdash = Rotz.cross(zdash);
      G4ThreeVector ydash = xdash.cross(zdash);
      SideRefVec1 = xdash.unit();
      SideRefVec2 = ydash.unit();
      SideRefVec3 = zdash.unit();
    }
  else if(Shape == "Ellipsoid")
    {
      G4double theta, phi, minphi, maxphi, middlephi;
      G4double answer, constant;

      constant = pi/(halfx*halfx) + pi/(halfy*halfy) + 
	twopi/(halfz*halfz);
      
      // simplified approach
      theta = posRndm->GenRandPosTheta();
      phi = posRndm->GenRandPosPhi();
      
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

      x = std::sin(theta)*std::cos(phi);
      y = std::sin(theta)*std::sin(phi);
      z = std::cos(theta);
      // x,y and z form a unit vector. Put this onto the ellipse.
      G4double lhs;
      // solve for x
      G4double numYinX = y/x;
      G4double numZinX = z/x;
      G4double tempxvar;	  
      tempxvar= 1./(halfx*halfx)+(numYinX*numYinX)/(halfy*halfy)
	+ (numZinX*numZinX)/(halfz*halfz);

      tempxvar = 1./tempxvar;
      G4double coordx = std::sqrt(tempxvar);
  
      //solve for y
      G4double numXinY = x/y;
      G4double numZinY = z/y;
      G4double tempyvar;
      tempyvar=(numXinY*numXinY)/(halfx*halfx)+1./(halfy*halfy)
	+(numZinY*numZinY)/(halfz*halfz);
      tempyvar = 1./tempyvar;
      G4double coordy = std::sqrt(tempyvar);
      
      //solve for z
      G4double numXinZ = x/z;
      G4double numYinZ = y/z;
      G4double tempzvar;
      tempzvar=(numXinZ*numXinZ)/(halfx*halfx)
	+(numYinZ*numYinZ)/(halfy*halfy)+1./(halfz*halfz);
      tempzvar = 1./tempzvar;
      G4double coordz = std::sqrt(tempzvar);

      lhs = std::sqrt((coordx*coordx)/(halfx*halfx) + 
		 (coordy*coordy)/(halfy*halfy) + 
		 (coordz*coordz)/(halfz*halfz));
      
      if(std::fabs(lhs-1.) > 0.001 && verbosityLevel >= 1)
	G4cout << "Error: theta, phi not really on ellipsoid" << G4endl;

      // coordx, coordy and coordz are all positive
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
      G4ThreeVector zdash(coordx,coordy,coordz);
      zdash = zdash.unit();
      G4ThreeVector xdash = Rotz.cross(zdash);
      G4ThreeVector ydash = xdash.cross(zdash);
      SideRefVec1 = xdash.unit();
      SideRefVec2 = ydash.unit();
      SideRefVec3 = zdash.unit();
    }
  else if(Shape == "Cylinder")
    {
      G4double AreaTop, AreaBot, AreaLat;
      G4double AreaTotal, prob1, prob2, prob3;
      G4double testrand;

      // User giver Radius and z-half length
      // Calculate surface areas, maybe move this to 
      // a different routine.

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
	    G4cout << AreaLat/AreaTotal << " " << prob3<<G4endl;
	  G4cout << "Error in prob3" << G4endl;
	}

      // Decide surface to calculate point on.

      testrand = G4UniformRand();
      if(testrand <= prob1)
	{
	  //Point on Top surface
	  z = halfz;
	  x = Radius + 100.;
	  y = Radius + 100.;
	  while(((x*x)+(y*y)) > (Radius*Radius))
	    {
	      x = posRndm->GenRandX();
	      y = posRndm->GenRandY();

	      x = x * 2. * Radius;
	      y = y * 2. * Radius;
	      x = x - Radius;
	      y = y - Radius;
	    }
	  // Cosine law
	  SideRefVec1 = Rotx;
	  SideRefVec2 = Roty;
	  SideRefVec3 = Rotz;
	}
      else if((testrand > prob1) && (testrand <= (prob1 + prob2)))
	{
	  //Point on Bottom surface
	  z = -halfz;
	  x = Radius + 100.;
	  y = Radius + 100.;
	  while(((x*x)+(y*y)) > (Radius*Radius))
	    {
	      x = posRndm->GenRandX();
	      y = posRndm->GenRandY();

	      x = x * 2. * Radius;
	      y = y * 2. * Radius;
	      x = x - Radius;
	      y = y - Radius;
	    }
	  // Cosine law
	  SideRefVec1 = Rotx;
	  SideRefVec2 = -Roty;
	  SideRefVec3 = -Rotz;
	}
      else if(testrand > (prob1+prob2))
	{
	  G4double rand;
	  //Point on Lateral Surface

	  rand = posRndm->GenRandPosPhi();
	  rand = rand * 2. * pi;

	  x = Radius * std::cos(rand);
	  y = Radius * std::sin(rand);

	  z = posRndm->GenRandZ();

	  z = z * 2. * halfz;
	  z = z - halfz;
	  
	  // Cosine law
	  G4ThreeVector zdash(x,y,0.);
	  zdash = zdash.unit();
	  G4ThreeVector xdash = Rotz.cross(zdash);
	  G4ThreeVector ydash = xdash.cross(zdash);
	  SideRefVec1 = xdash.unit();
	  SideRefVec2 = ydash.unit();
	  SideRefVec3 = zdash.unit();
	}
      else
	G4cout << "Error: testrand " << testrand << G4endl;

      RandPos.setX(x);
      RandPos.setY(y);
      RandPos.setZ(z);

    }
  else if(Shape == "Para")
    {
      G4double testrand;
      //Right Parallelepiped.
      // User gives x,y,z half lengths and ParAlpha
      // ParTheta and ParPhi
      // +x = <1, -x >1 & <2, +y >2 & <3, -y >3 &<4
      // +z >4 & < 5, -z >5 &<6.
      testrand = G4UniformRand();
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
      
      x = posRndm->GenRandX();
      y = posRndm->GenRandY();
      z = posRndm->GenRandZ();
      
      x = x * halfx * 2.;
      x = x - halfx;
      y = y * halfy * 2.;
      y = y - halfy;
      z = z * halfz * 2.;
      z = z - halfz;
      // Pick a side first
      if(testrand < Probs[0])
	{
	  // side is +x
	  x = halfx + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
	  y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector xdash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
			      halfz*std::tan(ParTheta)*std::sin(ParPhi), 
			      halfz/std::cos(ParPhi));
	  G4ThreeVector ydash(halfy*std::tan(ParAlpha), -halfy, 0.0);
	  xdash = xdash.unit();
	  ydash = ydash.unit();
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash.unit();
	  SideRefVec2 = ydash.unit();
	  SideRefVec3 = zdash.unit();
	}
      else if(testrand >= Probs[0] && testrand < Probs[1])
	{
	  // side is -x
	  x = -halfx + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
	  y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector xdash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
			      halfz*std::tan(ParTheta)*std::sin(ParPhi), 
			      halfz/std::cos(ParPhi));
	  G4ThreeVector ydash(halfy*std::tan(ParAlpha), halfy, 0.0);
	  xdash = xdash.unit();
	  ydash = ydash.unit();
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash.unit();
	  SideRefVec2 = ydash.unit();
	  SideRefVec3 = zdash.unit();
	}
      else if(testrand >= Probs[1] && testrand < Probs[2])
	{
	  // side is +y
	  x = x + z*std::tan(ParTheta)*std::cos(ParPhi) + halfy*std::tan(ParAlpha);
	  y = halfy + z*std::tan(ParTheta)*std::sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector ydash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
			      halfz*std::tan(ParTheta)*std::sin(ParPhi), 
			      halfz/std::cos(ParPhi));
	  ydash = ydash.unit();
	  G4ThreeVector xdash = Roty.cross(ydash);
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash.unit();
	  SideRefVec2 = -ydash.unit();
	  SideRefVec3 = -zdash.unit();
	}
      else if(testrand >= Probs[2] && testrand < Probs[3])
	{
	  // side is -y
	  x = x + z*std::tan(ParTheta)*std::cos(ParPhi) - halfy*std::tan(ParAlpha);
	  y = -halfy + z*std::tan(ParTheta)*std::sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector ydash(halfz*std::tan(ParTheta)*std::cos(ParPhi),
			      halfz*std::tan(ParTheta)*std::sin(ParPhi), 
			      halfz/std::cos(ParPhi));
	  ydash = ydash.unit();
	  G4ThreeVector xdash = Roty.cross(ydash);
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash.unit();
	  SideRefVec2 = ydash.unit();
	  SideRefVec3 = zdash.unit();
	}
      else if(testrand >= Probs[3] && testrand < Probs[4])
	{
	  // side is +z
	  z = halfz;
	  y = y + halfz*std::sin(ParPhi)*std::tan(ParTheta);
	  x = x + halfz*std::cos(ParPhi)*std::tan(ParTheta) + y*std::tan(ParAlpha);
	  // Cosine-law
	  SideRefVec1 = Rotx;
	  SideRefVec2 = Roty;
	  SideRefVec3 = Rotz;
	}
      else if(testrand >= Probs[4] && testrand < Probs[5])
	{
	  // side is -z
	  z = -halfz;
	  y = y - halfz*std::sin(ParPhi)*std::tan(ParTheta);
	  x = x - halfz*std::cos(ParPhi)*std::tan(ParTheta) + y*std::tan(ParAlpha);
	  // Cosine-law
	  SideRefVec1 = Rotx;
	  SideRefVec2 = -Roty;
	  SideRefVec3 = -Rotz;
	}
      else
	{
	  G4cout << "Error: testrand out of range" << G4endl;
	  if(verbosityLevel >= 1)
	    G4cout << "testrand=" << testrand << " Probs[5]=" << Probs[5] <<G4endl;
	}

      RandPos.setX(x);
      RandPos.setY(y);
      RandPos.setZ(z);
    }

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  if(verbosityLevel == 2)
    G4cout << "Raw position " << RandPos << G4endl;

  x=(RandPos.x()*Rotx.x())+(RandPos.y()*Roty.x())+(RandPos.z()*Rotz.x());
  y=(RandPos.x()*Rotx.y())+(RandPos.y()*Roty.y())+(RandPos.z()*Rotz.y());
  z=(RandPos.x()*Rotx.z())+(RandPos.y()*Roty.z())+(RandPos.z()*Rotz.z());
  
  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);

  // Translate
  particle_position = CentreCoords + RandPos;

  if(verbosityLevel >= 1)
    {
      if(verbosityLevel == 2)
	G4cout << "Rotated position " << RandPos << G4endl;
      G4cout << "Rotated and translated position " << particle_position << G4endl;
    }
  if(verbosityLevel == 2)
    {
      G4cout << "Reference vectors for cosine-law " << SideRefVec1 << " " << SideRefVec2 << " " << SideRefVec3 << G4endl;
    }
}

void G4SPSPosDistribution::GeneratePointsInVolume()
{
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  G4double x, y, z;
  x = y = z = 0.;
  if(SourcePosType != "Volume" && verbosityLevel >= 1)
    G4cout << "Error SourcePosType not Volume" << G4endl;
  //Private method to create points in a volume
  if(Shape == "Sphere")
    {
      x = Radius*2.;
      y = Radius*2.;
      z = Radius*2.;
      while(((x*x)+(y*y)+(z*z)) > (Radius*Radius))
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();
	  z = posRndm->GenRandZ();

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
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();
	  z = posRndm->GenRandZ();

	  x = (x*2.*halfx) - halfx;
	  y = (y*2.*halfy) - halfy;
	  z = (z*2.*halfz) - halfz;
	  
	  temp = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy))
	    + ((z*z)/(halfz*halfz));
	}
    }
  else if(Shape == "Cylinder")
    {
      x = Radius*2.;
      y = Radius*2.;
      while(((x*x)+(y*y)) > (Radius*Radius))
	{
	  x = posRndm->GenRandX();
	  y = posRndm->GenRandY();
	  z = posRndm->GenRandZ();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	  z = (z*2.*halfz) - halfz;
	}
    }
  else if(Shape == "Para")
    {
      x = posRndm->GenRandX();
      y = posRndm->GenRandY();
      z = posRndm->GenRandZ();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
      z = (z*2.*halfz) - halfz;
      x = x + z*std::tan(ParTheta)*std::cos(ParPhi) + y*std::tan(ParAlpha);
      y = y + z*std::tan(ParTheta)*std::sin(ParPhi);
      z = z;
    }
  else
    G4cout << "Error: Volume Shape Doesnt Exist" << G4endl;

  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);

  // Apply Rotation Matrix
  // x * Rotx, y * Roty and z * Rotz
  tempx = (x * Rotx.x()) + (y * Roty.x()) + (z * Rotz.x());
  tempy = (x * Rotx.y()) + (y * Roty.y()) + (z * Rotz.y());
  tempz = (x * Rotx.z()) + (y * Roty.z()) + (z * Rotz.z());

  RandPos.setX(tempx);
  RandPos.setY(tempy);
  RandPos.setZ(tempz);

  // Translate
  particle_position = CentreCoords + RandPos;

  if(verbosityLevel == 2)
    {
      G4cout << "Raw position " << x << "," << y << "," << z << G4endl;
      G4cout << "Rotated position " << RandPos << G4endl;
    }
  if(verbosityLevel >= 1)
    G4cout << "Rotated and translated position " << particle_position << G4endl;

  // Cosine-law (not a good idea to use this here)
  G4ThreeVector zdash(tempx,tempy,tempz);
  zdash = zdash.unit();
  G4ThreeVector xdash = Rotz.cross(zdash);
  G4ThreeVector ydash = xdash.cross(zdash);
  SideRefVec1 = xdash.unit();
  SideRefVec2 = ydash.unit();
  SideRefVec3 = zdash.unit();

  if(verbosityLevel == 2)
    {
      G4cout << "Reference vectors for cosine-law " << SideRefVec1 << " " << SideRefVec2 << " " << SideRefVec3 << G4endl;
    } 
}

G4bool G4SPSPosDistribution::IsSourceConfined()
{
  // Method to check point is within the volume specified
  if(Confine == false)
    G4cout << "Error: Confine is false" << G4endl;
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr;
  ptr = &null;

  // Check particle_position is within VolName, if so true, 
  // else false
  G4VPhysicalVolume *theVolume;
  theVolume=gNavigator->LocateGlobalPointAndSetup(particle_position,ptr,true);
  G4String theVolName = theVolume->GetName();
  if(theVolName == VolName)
    {
      if(verbosityLevel >= 1)
	G4cout << "Particle is in volume " << VolName << G4endl;
      return(true);
    }
  else
    return(false);
}

G4ThreeVector G4SPSPosDistribution::GenerateOne()
{
  //
  G4bool srcconf = false;
  G4int LoopCount = 0;
  while(srcconf == false)
    {
      if(SourcePosType == "Point")
	GeneratePointSource();
      else if(SourcePosType == "Beam")
	GeneratePointsInBeam();
      else if(SourcePosType == "Plane")
	GeneratePointsInPlane();
      else if(SourcePosType == "Surface")
	GeneratePointsOnSurface();
      else if(SourcePosType == "Volume")
	GeneratePointsInVolume();
      else
	{
	  G4cout << "Error: SourcePosType undefined" << G4endl;
	  G4cout << "Generating point source" << G4endl;
	  GeneratePointSource();
	}
      if(Confine == true)
	{
	  srcconf = IsSourceConfined();
	  // if source in confined srcconf = true terminating the loop
	  // if source isnt confined srcconf = false and loop continues
	}
      else if(Confine == false)
	srcconf = true; // terminate loop
      LoopCount++;
      if(LoopCount == 100000)
	{
	  G4cout << "*************************************" << G4endl;
	  G4cout << "LoopCount = 100000" << G4endl;
	  G4cout << "Either the source distribution >> confinement" << G4endl;
	  G4cout << "or any confining volume may not overlap with" << G4endl;
	  G4cout << "the source distribution or any confining volumes" << G4endl;
	  G4cout << "may not exist"<< G4endl;
	  G4cout << "If you have set confine then this will be ignored" <<G4endl;
	  G4cout << "for this event." << G4endl;
	  G4cout << "*************************************" << G4endl;
	  srcconf = true; //Avoids an infinite loop
	}
    }
  return particle_position;
}




