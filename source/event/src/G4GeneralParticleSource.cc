///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSource.cc
//
// Version:      1.1
// Date:         18/10/00
// Author:       C Ferguson, F Lei, P Truscott
// Organisation: University of Southampton / DERA
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 28 February 2000, C Ferguson, Created.
//
// Version 1.1, 18 October 2000, Modified to inherit from G4VPrimaryGenerator.
// New name at the request of M. Asai.
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <math.h>
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

#include "G4GeneralParticleSource.hh"

G4GeneralParticleSource::G4GeneralParticleSource()
{
  // Initialise all variables
  // Position distribution Variables

  NumberOfParticlesToBeGenerated = 1;
  particle_definition = NULL;
  G4ThreeVector zero;
  particle_momentum_direction = (G4ParticleMomentum)zero;
  particle_energy = 0.0;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = 0.0;

  SourcePosType = "Point";
  Shape = "NULL";
  halfx = 0.;
  halfy = 0.;
  halfz = 0.;
  Radius = 0.;
  Radius0 = 0.;
  ParAlpha = 0.;
  ParTheta = 0.;
  ParPhi = 0.;
  CentreCoords = (0.,0.,0.);
  Rotx = HepXHat;
  Roty = HepYHat;
  Rotz = HepZHat;
  Confine = false; //If true confines source distribution to VolName
  VolName = "NULL";
  SideRefVec1 = HepXHat; // x-axis
  SideRefVec2 = HepYHat; // y-axis
  SideRefVec3 = HepZHat; // z-axis
  UserWRTSurface = false; // Any user-defined distribution is wrt co-ordinate
  // axes.

  // Angular distribution variables.
  AngDistType = "iso"; 
  AngRef1 = (1.,0.,0.);
  AngRef2 = (0.,1.,0.);
  AngRef3 = (0.,0.,1.);
  MinTheta = 0.;
  MaxTheta = pi;
  MinPhi = 0.;
  MaxPhi = twopi;
  UserDistType = "NULL";
  IPDFThetaExist = false;
  IPDFPhiExist = false;

  // Energy Distribution variables
  EnergyDisType = "Mono";
  MonoEnergy = 1*MeV;
  Emin = 0.;
  Emax = 0.;
  alpha = 0.;
  Ezero = 0.;
  Temp = 0.;
  grad = 0.;
  cept = 0.;
  EnergySpec = true; // true - energy spectra, false - momentum spectra
  DiffSpec = true;  // true - differential spec, false integral spec
  IntType = "NULL"; // Interpolation type
  IPDFEnergyExist = false;
  IPDFArbExist = false;

  // Bias variables
  XBias = false;
  IPDFXBias = false;
  YBias = false;
  IPDFYBias = false;
  ZBias = false;
  IPDFZBias = false;
  ThetaBias = false;
  IPDFThetaBias = false;
  PhiBias = false;
  IPDFPhiBias = false;
  EnergyBias = false;
  IPDFEnergyBias = false;

  ArbEmin = 0.;
  ArbEmax = 0.;

  // verbosity
  verbosityLevel = 0;

  theMessenger = new G4GeneralParticleSourceMessenger(this);
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}

G4GeneralParticleSource::~G4GeneralParticleSource()
{
  delete theMessenger;
}

void G4GeneralParticleSource::SetPosDisType(G4String PosType)
{
  SourcePosType = PosType;
}

void G4GeneralParticleSource::SetPosDisShape(G4String shapeType)
{
  Shape = shapeType;
}

void G4GeneralParticleSource::SetCentreCoords(G4ThreeVector coordsOfCentre)
{
  CentreCoords = coordsOfCentre;
}

void G4GeneralParticleSource::SetPosRot1(G4ThreeVector posrot1)
{
  // This should be x'
  Rotx = posrot1;
  if(verbosityLevel == 2)
    {
      G4cout << "Vector x' " << Rotx << G4endl;
    }
  GenerateRotationMatrices();
}

void G4GeneralParticleSource::SetPosRot2(G4ThreeVector posrot2)
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

void G4GeneralParticleSource::SetHalfX(G4double xhalf)
{
  halfx = xhalf;
}

void G4GeneralParticleSource::SetHalfY(G4double yhalf)
{
  halfy = yhalf;
}

void G4GeneralParticleSource::SetHalfZ(G4double zhalf)
{
  halfz = zhalf;
}

void G4GeneralParticleSource::SetRadius(G4double rad)
{
  Radius = rad;
}

void G4GeneralParticleSource::SetRadius0(G4double rad)
{
  Radius0 = rad;
}

void G4GeneralParticleSource::SetParAlpha(G4double paralp)
{
  ParAlpha = paralp;
}

void G4GeneralParticleSource::SetParTheta(G4double parthe)
{
  ParTheta = parthe;
}

void G4GeneralParticleSource::SetParPhi(G4double parphi)
{
  ParPhi = parphi;
}

void G4GeneralParticleSource::GenerateRotationMatrices()
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

void G4GeneralParticleSource::ConfineSourceToVolume(G4String Vname)
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
    G4cout << PVStore->length() << G4endl;
  while (!found && i<PVStore->length()) {
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

void G4GeneralParticleSource::GeneratePointSource()
{
  // Generates Points given the point source.
  if(SourcePosType == "Point")
    particle_position = CentreCoords;
  else
    if(verbosityLevel >= 1)
      G4cout << "Error SourcePosType is not set to Point" << G4endl;
}

void G4GeneralParticleSource::GeneratePointsInPlane()
{
  G4double x, y, z;
  G4double expression;
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;

  z = 0.;

  if(SourcePosType != "Plane" && verbosityLevel >= 1)
    G4cout << "Error: SourcePosType is not Plane" << G4endl;

  // Private Method to create points in a plane
  if(Shape == "Circle")
    {
      x = Radius + 100.;
      y = Radius + 100.;
      while(sqrt((x*x) + (y*y)) > Radius)
	{
	  x = GenRandX();
	  y = GenRandY();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	}
    }
  else if(Shape == "Annulus")
    {
      x = Radius + 100.;
      y = Radius + 100.;
      while(sqrt((x*x) + (y*y)) > Radius || sqrt((x*x) + (y*y)) < Radius0 )
	{
	  x = GenRandX();
	  y = GenRandY();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	}
    }
  else if(Shape == "Ellipse")
    {
      expression = 20.;
      while(expression > 1.)
	{
	  x = GenRandX();
	  y = GenRandY();

	  x = (x*2.*halfx) - halfx;
	  y = (y*2.*halfy) - halfy;

	  expression = ((x*x)/(halfx*halfx)) + ((y*y)/(halfy*halfy));
	}
    }
  else if(Shape == "Square")
    {
      x = GenRandX();
      y = GenRandY();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
    }
  else if(Shape == "Rectangle")
    {
      x = GenRandX();
      y = GenRandY();
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

void G4GeneralParticleSource::GeneratePointsOnSurface()
{
  //Private method to create points on a surface
  G4double theta, phi;
  G4double x, y, z;
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;

  if(SourcePosType != "Surface" && verbosityLevel >= 1)
    G4cout << "Error SourcePosType not Surface" << G4endl;

  if(Shape == "Sphere")
    {
      G4double tantheta, tanphi;
      theta = GenRandTheta();
      phi = GenRandPhi();

      theta = acos(1. - 2.*theta); // theta isotropic
      phi = phi * 2. * pi;
      tantheta = tan(theta);
      
      x = Radius * sin(theta) * cos(phi);
      y = Radius * sin(theta) * sin(phi);
      z = Radius * cos(theta);
      
      RandPos.setX(x);
      RandPos.setY(y);
      RandPos.setZ(z);

      // Cosine-law (not a good idea to use this here)
      G4ThreeVector zdash(x,y,z);
      zdash = zdash.unit();
      G4ThreeVector xdash = Rotz.cross(zdash);
      G4ThreeVector ydash = xdash.cross(zdash);
      SideRefVec1 = xdash;
      SideRefVec2 = ydash;
      SideRefVec3 = zdash;
    }
  else if(Shape == "Ellipsoid")
    {
      G4double testrand, theta, phi, minphi, maxphi, middlephi;
      G4double answer, constant;

      constant = pi/(halfx*halfx) + pi/(halfy*halfy) + 
	twopi/(halfz*halfz);
      
      // simplified approach
      theta = GenRandTheta();
      phi = GenRandPhi();

      theta = acos(1. - 2.*theta);
      minphi = 0.;
      maxphi = twopi;
      while(maxphi-minphi > 0.)
	{
	  middlephi = (maxphi+minphi)/2.;
	  answer = (1./(halfx*halfx))*(middlephi/2. + sin(2*middlephi)/4.)
	    + (1./(halfy*halfy))*(middlephi/2. - sin(2*middlephi)/4.)
	       + middlephi/(halfz*halfz);
	  answer = answer/constant;
	  if(answer > phi) maxphi = middlephi;
	  if(answer < phi) minphi = middlephi;
	  if(fabs(answer-phi) <= 0.00001)
	    {
	      minphi = maxphi +1;
	      phi = middlephi;
	    }
	}

      x = sin(theta)*cos(phi);
      y = sin(theta)*sin(phi);
      z = cos(theta);
      // x,y and z form a unit vector. Put this onto the ellipse.
      G4double lhs;
      // solve for x
      G4double numYinX = y/x;
      G4double numZinX = z/x;
      G4double tempxvar;	  
      tempxvar= 1./(halfx*halfx)+(numYinX*numYinX)/(halfy*halfy)
	+ (numZinX*numZinX)/(halfz*halfz);

      tempxvar = 1./tempxvar;
      G4double coordx = sqrt(tempxvar);
  
      //solve for y
      G4double numXinY = x/y;
      G4double numZinY = z/y;
      G4double tempyvar;
      tempyvar=(numXinY*numXinY)/(halfx*halfx)+1./(halfy*halfy)
	+(numZinY*numZinY)/(halfz*halfz);
      tempyvar = 1./tempyvar;
      G4double coordy = sqrt(tempyvar);
      
      //solve for z
      G4double numXinZ = x/z;
      G4double numYinZ = y/z;
      G4double tempzvar;
      tempzvar=(numXinZ*numXinZ)/(halfx*halfx)
	+(numYinZ*numYinZ)/(halfy*halfy)+1./(halfz*halfz);
      tempzvar = 1./tempzvar;
      G4double coordz = sqrt(tempzvar);

      lhs = sqrt((coordx*coordx)/(halfx*halfx) + 
		 (coordy*coordy)/(halfy*halfy) + 
		 (coordz*coordz)/(halfz*halfz));
      
      if(fabs(lhs-1.) > 0.001 && verbosityLevel >= 1)
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
      SideRefVec1 = xdash;
      SideRefVec2 = ydash;
      SideRefVec3 = zdash;
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
      if(fabs(prob3 - (AreaLat/AreaTotal)) >= 0.001)
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
	      x = GenRandX();
	      y = GenRandY();

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
	      x = GenRandX();
	      y = GenRandY();

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

	  rand = G4UniformRand();
	  rand = rand * 2. * pi;

	  x = Radius * cos(rand);
	  y = Radius * sin(rand);

	  z = GenRandZ();

	  z = z * 2. * halfz;
	  z = z - halfz;
	  
	  // Cosine law
	  G4ThreeVector zdash(x,y,0.);
	  zdash = zdash.unit();
	  G4ThreeVector xdash = Rotz.cross(zdash);
	  G4ThreeVector ydash = xdash.cross(zdash);
	  SideRefVec1 = xdash;
	  SideRefVec2 = ydash;
	  SideRefVec3 = zdash;
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
      
      x = GenRandX();
      y = GenRandY();
      z = GenRandZ();
      
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
	  x = halfx + z*tan(ParTheta)*cos(ParPhi) + y*tan(ParAlpha);
	  y = y + z*tan(ParTheta)*sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector xdash(halfz*tan(ParTheta)*cos(ParPhi),
			      halfz*tan(ParTheta)*sin(ParPhi), 
			      halfz/cos(ParPhi));
	  G4ThreeVector ydash(halfy*tan(ParAlpha), halfy, 0.0);
	  xdash = xdash.unit();
	  ydash = ydash.unit();
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash;
	  SideRefVec2 = ydash;
	  SideRefVec3 = zdash;
	}
      else if(testrand >= Probs[0] && testrand < Probs[1])
	{
	  // side is -x
	  x = -halfx + z*tan(ParTheta)*cos(ParPhi) + y*tan(ParAlpha);
	  y = y + z*tan(ParTheta)*sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector xdash(halfz*tan(ParTheta)*cos(ParPhi),
			      halfz*tan(ParTheta)*sin(ParPhi), 
			      halfz/cos(ParPhi));
	  G4ThreeVector ydash(halfy*tan(ParAlpha), halfy, 0.0);
	  xdash = xdash.unit();
	  ydash = ydash.unit();
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash;
	  SideRefVec2 = -ydash;
	  SideRefVec3 = -zdash;
	}
      else if(testrand >= Probs[1] && testrand < Probs[2])
	{
	  // side is +y
	  x = x + z*tan(ParTheta)*cos(ParPhi) + halfy*tan(ParAlpha);
	  y = halfy + z*tan(ParTheta)*sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector ydash(halfz*tan(ParTheta)*cos(ParPhi),
			      halfz*tan(ParTheta)*sin(ParPhi), 
			      halfz/cos(ParPhi));
	  ydash = ydash.unit();
	  G4ThreeVector xdash = Roty.cross(ydash);
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash;
	  SideRefVec2 = ydash;
	  SideRefVec3 = zdash;
	}
      else if(testrand >= Probs[2] && testrand < Probs[3])
	{
	  // side is -y
	  x = x + z*tan(ParTheta)*cos(ParPhi) - halfy*tan(ParAlpha);
	  y = -halfy + z*tan(ParTheta)*sin(ParPhi);
	  z = z;
	  // Cosine-law
	  G4ThreeVector ydash(halfz*tan(ParTheta)*cos(ParPhi),
			      halfz*tan(ParTheta)*sin(ParPhi), 
			      halfz/cos(ParPhi));
	  ydash = ydash.unit();
	  G4ThreeVector xdash = Roty.cross(ydash);
	  G4ThreeVector zdash = xdash.cross(ydash);
	  SideRefVec1 = xdash;
	  SideRefVec2 = -ydash;
	  SideRefVec3 = -zdash;
	}
      else if(testrand >= Probs[3] && testrand < Probs[4])
	{
	  // side is +z
	  z = halfz;
	  y = y + halfz*sin(ParPhi)*tan(ParTheta);
	  x = x + halfz*cos(ParPhi)*tan(ParTheta) + y*tan(ParAlpha);
	  // Cosine-law
	  SideRefVec1 = Rotx;
	  SideRefVec2 = Roty;
	  SideRefVec3 = Rotz;
	}
      else if(testrand >= Probs[4] && testrand < Probs[5])
	{
	  // side is -z
	  z = -halfz;
	  y = y - halfz*sin(ParPhi)*tan(ParTheta);
	  x = x - halfz*cos(ParPhi)*tan(ParTheta) + y*tan(ParAlpha);
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

void G4GeneralParticleSource::GeneratePointsInVolume()
{
  G4ThreeVector RandPos;
  G4double tempx, tempy, tempz;
  G4double x, y, z;

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
	  x = GenRandX();
	  y = GenRandY();
	  z = GenRandZ();

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
	  x = GenRandX();
	  y = GenRandY();
	  z = GenRandZ();

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
	  x = GenRandX();
	  y = GenRandY();
	  z = GenRandZ();

	  x = (x*2.*Radius) - Radius;
	  y = (y*2.*Radius) - Radius;
	  z = (z*2.*halfz) - halfz;
	}
    }
  else if(Shape == "Para")
    {
      x = GenRandX();
      y = GenRandY();
      z = GenRandZ();
      x = (x*2.*halfx) - halfx;
      y = (y*2.*halfy) - halfy;
      z = (z*2.*halfz) - halfz;
      x = x + z*tan(ParTheta)*cos(ParPhi) + y*tan(ParAlpha);
      y = y + z*tan(ParTheta)*sin(ParPhi);
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
  SideRefVec1 = xdash;
  SideRefVec2 = ydash;
  SideRefVec3 = zdash;

  if(verbosityLevel == 2)
    {
      G4cout << "Reference vectors for cosine-law " << SideRefVec1 << " " << SideRefVec2 << " " << SideRefVec3 << G4endl;
    } 
}

G4bool G4GeneralParticleSource::IsSourceConfined()
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

// Now follows the angular distribution methods.

void G4GeneralParticleSource::SetAngDistType(G4String atype)
{
  if(atype != "iso" && atype != "cos" && atype != "user")
    G4cout << "Error, distribution must be iso, cos or user" << G4endl;
  else
    AngDistType = atype;
}

void G4GeneralParticleSource::DefineAngRefAxes(G4String refname, G4ThreeVector ref)
{

  if(refname == "angref1")
    AngRef1 = ref; // x'
  else if(refname == "angref2")
    AngRef2 = ref; // vector in x'y' plane

  // User defines x' (AngRef1) and a vector in the x'y'
  // plane (AngRef2). Then, AngRef1 x AngRef2 = AngRef3
  // the z' vector. Then, AngRef3 x AngRef1 = AngRef2
  // which will now be y'.

  AngRef3 = AngRef1.cross(AngRef2); // z'
  AngRef2 = AngRef3.cross(AngRef1); // y'
  if(verbosityLevel == 2)
    {
      G4cout << "Angular distribution rotation axes " << AngRef1 << " " << AngRef2 << " " << AngRef3 << G4endl;
    }
}

void G4GeneralParticleSource::SetMinTheta(G4double mint)
{
  MinTheta = mint;
}

void G4GeneralParticleSource::SetMinPhi(G4double minp)
{
  MinPhi = minp;
}

void G4GeneralParticleSource::SetMaxTheta(G4double maxt)
{
  MaxTheta = maxt;
}

void G4GeneralParticleSource::SetMaxPhi(G4double maxp)
{
  MaxPhi = maxp;
}

void G4GeneralParticleSource::UserDefAngTheta(G4ThreeVector input)
{
  if(UserDistType == "NULL") UserDistType = "theta";
  if(UserDistType == "phi") UserDistType = "both";  
  G4double thi, val;
  thi = input.x();
  val = input.y();
  if(verbosityLevel >= 1)
    G4cout << "In UserDefAngTheta" << G4endl;
  UDefThetaH.InsertValues(thi, val);
}

void G4GeneralParticleSource::UserDefAngPhi(G4ThreeVector input)
{
  if(UserDistType == "NULL") UserDistType = "phi";
  if(UserDistType == "theta") UserDistType = "both";  
  G4double phhi, val;
  phhi = input.x();
  val = input.y();
  if(verbosityLevel >= 1)
    G4cout << "In UserDefAngPhi" << G4endl;
  UDefPhiH.InsertValues(phhi, val); 
}

void G4GeneralParticleSource::SetUserWRTSurface(G4bool wrtSurf)
{
  // if UserWRTSurface = true then the user wants momenta with respect
  // to the surface normals.
  // When doing this theta has to be 0-90 only otherwise there will be
  // errors, which currently are flagged anywhere.
  UserWRTSurface = wrtSurf;
}

void G4GeneralParticleSource::GenerateIsotropicFlux()
{
  // generates isotropic flux.
  // No vectors are needed.
  G4double rndm, rndm2;
  G4double px, py, pz, pmag;

  // generate rand nos. but make sure in theta/phi limits
  Theta = 10.;  // This is well above the pi upper limit.
  Phi = 10.;    // Still well above 2pi upper limit.

  while(Theta > MaxTheta || Theta < MinTheta)
    {
      rndm = GenRandTheta();
      Theta = acos(1. - 2.*rndm);
    }

  while(Phi > MaxPhi || Phi < MinPhi)
    {
      rndm2 = GenRandPhi();
      Phi = twopi * rndm2;
    }

  px = sin(Theta) * cos(Phi);
  py = sin(Theta) * sin(Phi);
  pz = cos(Theta);

  px = -px;
  py = -py;
  pz = -pz;

  pmag = sqrt((px*px) + (py*py) + (pz*pz));

  px = px/pmag;
  py = py/pmag;
  pz = pz/pmag;

  //Need to rotate the particle_momentum_direction round such that the 
  //start position is at the north pole, so that all particles go into
  // the centre of the geomtries.
  // particle_position stores the particle's position.
  // use particle position to calculate theta and phi - posthe and posphi

  G4double posthe, posphi;
  //  G4cout << "particle_position " << particle_position << G4endl;
  G4double tx, ty, tz, tt;
  tx = particle_position.x();
  ty = particle_position.y();
  tz = particle_position.z();
  tt = sqrt(tx*tx + ty*ty + tz*tz);
  tx = tx/tt;
  ty = ty/tt;
  tz = tz/tt;
  //  G4cout << "unit position " << tx << " " << ty << " " << tz << G4endl;
  posthe = acos(tz);
  posphi = acos(tx/sin(posthe));
  //  G4cout << "Posthe and posphi " << posthe << " " << posphi << G4endl;
  G4double finx, finy, finz;
  finx = (px*cos(posthe)*cos(posphi)) - (py*sin(posphi)) + (pz*sin(posthe)*cos(posphi));
  finy = (px*cos(posthe)*sin(posphi)) + (py*cos(posphi)) + (pz*sin(posthe)*sin(posphi));
  finz = (-px*sin(posthe)) + (pz*cos(posthe));
  //  G4cout << "finx... " << finx << " " << finy << " " << finz << G4endl;

  particle_momentum_direction.setX(finx);
  particle_momentum_direction.setY(finy);
  particle_momentum_direction.setZ(finz);

  // particle_momentum_direction now holds unit momentum vector.
  if(verbosityLevel >= 1)
    G4cout << "Generating isotropic vector: " << particle_momentum_direction << G4endl;
}

void G4GeneralParticleSource::GenerateCosineLawFlux()
{
  // Method to generate flux distributed with a cosine law
  // such as is in TIMM.
  G4double px, py, pz, pmag;
  G4double rndm, rndm2;
  G4double resultx, resulty, resultz;

  Theta = 10.; // Well above MaxTheta
  Phi = 10.;   // Well above MaxPhi

  while(Theta > MaxTheta || Theta < MinTheta)
    {
      rndm = GenRandTheta();
      Theta = asin(sqrt(rndm));
    }

  while(Phi > MaxPhi || Phi < MinPhi)
    {
      rndm2 = GenRandPhi();
      Phi = twopi * rndm2;
    }

  px = sin(Theta) * cos(Phi);
  py = sin(Theta) * sin(Phi);
  pz = cos(Theta);

  px = -px;
  py = -py;
  pz = -pz;
  pmag = sqrt((px*px) + (py*py) + (pz*pz));
  G4double pxh = px/pmag;
  G4double pyh = py/pmag;
  G4double pzh = pz/pmag;

  if(verbosityLevel == 2)
    {
      G4cout <<"SideRefVecs " <<SideRefVec1<<SideRefVec2<<SideRefVec3<<G4endl;
      G4cout <<"Raw Unit vector "<<pxh<<","<<pyh<<","<<pzh<<G4endl;
    }
  resultx = (pxh*SideRefVec1.x()) + (pyh*SideRefVec2.x()) + 
    (pzh*SideRefVec3.x());

  resulty = (pxh*SideRefVec1.y()) + (pyh*SideRefVec2.y()) + 
    (pzh*SideRefVec3.y());

  resultz = (pxh*SideRefVec1.z()) + (pyh*SideRefVec2.z()) + 
    (pzh*SideRefVec3.z());

  G4double ResMag = sqrt((resultx*resultx) + (resulty*resulty) + (resultz*resultz));
  resultx = resultx/ResMag;
  resulty = resulty/ResMag;
  resultz = resultz/ResMag;
  
  particle_momentum_direction.setX(resultx);
  particle_momentum_direction.setY(resulty);
  particle_momentum_direction.setZ(resultz);

  // particle_momentum_direction now contains unit momentum vector.
  if(verbosityLevel >= 1)
    {
      G4cout << "Resultant cosine-law unit momentum vector " << particle_momentum_direction << G4endl;
    }
}

void G4GeneralParticleSource::GenerateUserDefFlux()
{
  G4double rndm, px, py, pz, pmag;

  if(UserDistType == "NULL")
    G4cout << "Error: UserDistType undefined" << G4endl;
  else if(UserDistType == "theta")
    {
      Theta = GenerateUserDefTheta();
      Phi = 10.;
      while(Phi > MaxPhi || Phi < MinPhi)
	{
	  rndm = GenRandPhi();
	  Phi = twopi * rndm;
	}
    }
  else if(UserDistType == "phi")
    {
      Theta = 10.;
      while(Theta > MaxTheta || Theta < MinTheta)
	{
	  rndm = GenRandTheta();
	  Theta = acos(1. - (2. * rndm));
	}
      Phi = GenerateUserDefPhi();
    }
  else if(UserDistType == "both")
    {
      Theta = GenerateUserDefTheta();
      Phi = GenerateUserDefPhi();
    }

  px = -sin(Theta) * cos(Phi);
  py = -sin(Theta) * sin(Phi);
  pz = -cos(Theta);

  pmag = sqrt((px*px) + (py*py) + (pz*pz));

  if(UserWRTSurface == false)
    {
      particle_momentum_direction.setX(px/pmag);
      particle_momentum_direction.setY(py/pmag);
      particle_momentum_direction.setZ(pz/pmag);
    }
  else if(UserWRTSurface == true)
    {
      G4double pxh = px/pmag;
      G4double pyh = py/pmag;
      G4double pzh = pz/pmag;

      if(verbosityLevel == 2)
	{
	  G4cout <<"SideRefVecs " <<SideRefVec1<<SideRefVec2<<SideRefVec3<<G4endl;
	  G4cout <<"Raw Unit vector "<<pxh<<","<<pyh<<","<<pzh<<G4endl;
	}
      G4double resultx = (pxh*SideRefVec1.x()) + (pyh*SideRefVec2.x()) + 
	(pzh*SideRefVec3.x());

      G4double resulty = (pxh*SideRefVec1.y()) + (pyh*SideRefVec2.y()) + 
	(pzh*SideRefVec3.y());

      G4double resultz = (pxh*SideRefVec1.z()) + (pyh*SideRefVec2.z()) + 
	(pzh*SideRefVec3.z());

      G4double ResMag = sqrt((resultx*resultx) + (resulty*resulty) + (resultz*resultz));
      resultx = resultx/ResMag;
      resulty = resulty/ResMag;
      resultz = resultz/ResMag;
      
      particle_momentum_direction.setX(resultx);
      particle_momentum_direction.setY(resulty);
      particle_momentum_direction.setZ(resultz);
    }

  // particle_momentum_direction now contains unit momentum vector.
  if(verbosityLevel >= 1)
    {
      G4cout << "Final User Defined momentum vector " << particle_momentum_direction << G4endl;
    }
}

G4double G4GeneralParticleSource::GenerateUserDefTheta()
{
  // Create cumulative histogram if not already done so. Then use RandFlat
  //::shoot to generate the output Theta value.
  if(UserDistType == "NULL" || UserDistType == "phi")
    {
      // No user defined theta distribution
      G4cout << "Error ***********************" << G4endl;
      G4cout << "UserDistType = " << UserDistType << G4endl;
    }
  else
    {
      // UserDistType = theta or both and so a theta distribution
      // is defined. This should be integrated if not already done.
      if(IPDFThetaExist == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(UDefThetaH.GetVectorLength());
	  bins[0] = UDefThetaH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = UDefThetaH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = UDefThetaH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = UDefThetaH(size_t(ii)) + vals[ii-1];
	      sum = sum + UDefThetaH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFThetaH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFThetaExist = true
	  IPDFThetaExist = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      return(IPDFThetaH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenerateUserDefPhi()
{
  // Create cumulative histogram if not already done so. Then use RandFlat
  //::shoot to generate the output Theta value.

  if(UserDistType == "NULL" || UserDistType == "theta")
    {
      // No user defined phi distribution
      G4cout << "Error ***********************" << G4endl;
      G4cout << "UserDistType = " << UserDistType << G4endl;
    }
  else
    {
      // UserDistType = phi or both and so a phi distribution
      // is defined. This should be integrated if not already done.
      if(IPDFPhiExist == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(UDefPhiH.GetVectorLength());
	  bins[0] = UDefPhiH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = UDefPhiH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = UDefPhiH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = UDefPhiH(size_t(ii)) + vals[ii-1];
	      sum = sum + UDefPhiH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFPhiH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFPhiExist = true
	  IPDFPhiExist = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      return(IPDFPhiH.GetEnergy(rndm));
    }
}

// Now follows Energy distribution Methods

void G4GeneralParticleSource::SetEnergyDisType(G4String DisType)
{
  EnergyDisType = DisType;
}

void G4GeneralParticleSource::SetEmin(G4double emi)
{
  Emin = emi;
}

void G4GeneralParticleSource::SetEmax(G4double ema)
{
  Emax = ema;
}

void G4GeneralParticleSource::SetMonoEnergy(G4double menergy)
{
  MonoEnergy = menergy;
  Emin = menergy;
  Emax = menergy;
}

void G4GeneralParticleSource::SetAlpha(G4double alp)
{
  alpha = alp;
}

void G4GeneralParticleSource::SetTemp(G4double tem)
{
  Temp = tem;
}

void G4GeneralParticleSource::SetEzero(G4double eze)
{
  Ezero = eze;
}

void G4GeneralParticleSource::SetGradient(G4double gr)
{
  grad = gr;
}

void G4GeneralParticleSource::SetInterCept(G4double c)
{
  cept = c;
}

void G4GeneralParticleSource::UserEnergyHisto(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  if(verbosityLevel == 2)
    G4cout << " " << ehi << " " << val << G4endl;
  if(verbosityLevel == 2)
    G4cout << "In UserEnergyHisto" << G4endl;
  UDefEnergyH.InsertValues(ehi, val);
  Emax = ehi;
}

void G4GeneralParticleSource::ArbEnergyHisto(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  if(verbosityLevel == 2)
    G4cout << "In ArbEnergyHisto" << G4endl;
  ArbEnergyH.InsertValues(ehi, val);
}

void G4GeneralParticleSource::EpnEnergyHisto(G4ThreeVector input)
{
  G4double elo, ehi, val;
  ehi = input.x();
  val = input.y();
  if(verbosityLevel == 2)
    G4cout << "In EpnEnergyHisto" << G4endl;
  EpnEnergyH.InsertValues(ehi, val);
  Emax = ehi;
  Epnflag = true;
}

void G4GeneralParticleSource::Calculate()
{
  if(EnergyDisType == "Cdg")
    CalculateCdgSpectrum();
  else if(EnergyDisType == "Bbody")
    CalculateBbodySpectrum();
}

void G4GeneralParticleSource::CalculateCdgSpectrum()
{
  // This uses the spectrum from The INTEGRAL Mass Model (TIMM)
  // to generate a Cosmic Diffuse X/gamma ray spectrum.
  G4double pfact[2] = {8.5, 112};
  G4double spind[2] = {1.4, 2.3};
  G4double ene_line[3] = {1.*keV, 18.*keV, 1E6*keV};
  G4int n_par;

  ene_line[0] = Emin;
  if(Emin < 18*keV)
    {
      n_par = 2;
      ene_line[2] = Emax;
      if(Emax < 18*keV)
	{
	  n_par = 1;
	  ene_line[1] = Emax;
	}
    }
  else
    {
      n_par = 1;
      pfact[0] = 112.;
      spind[0] = 2.3;
      ene_line[1] = Emax;
    }
  
  // Create a cumulative histogram.
  CDGhist[0] = 0.;
  G4double omalpha;
  G4int i = 0;

  while(i < n_par)
    {
      omalpha = 1. - spind[i];
      CDGhist[i+1] = CDGhist[i] + (pfact[i]/omalpha)*
	(pow(ene_line[i+1],omalpha)-pow(ene_line[i],omalpha));
      i++;
    }
  
  // Normalise histo and divide by 1000 to make MeV.
  i = 0;
  while(i < n_par)
    {
      CDGhist[i+1] = CDGhist[i+1]/CDGhist[n_par];
      //      G4cout << CDGhist[i] << CDGhist[n_par] << G4endl;
      i++;
    }
}

void G4GeneralParticleSource::CalculateBbodySpectrum()
{
  // create bbody spectrum
  // Proved very hard to integrate indefinitely, so different
  // method. User inputs emin, emax and T. These are used to
  // create a 10,000 bin histogram.
  // Use photon density spectrum = 2 nu**2/c**2 * (exp(h nu/kT)-1)
  // = 2 E**2/h**2c**2 times the exponential
  G4double erange = Emax - Emin;
  G4double steps = erange/10000.;
  G4double Bbody_y[10000];
  G4double k = 8.6181e-11; //Boltzmann const in MeV/K
  G4double h = 4.1362e-21; // Plancks const in MeV s
  G4double c = 3e8; // Speed of light
  G4double h2 = h*h;
  G4double c2 = c*c;
  G4int count = 0;
  G4double sum = 0.;
  BBHist[0] = 0.;
  while(count < 10000)
    {
      Bbody_x[count] = Emin + G4double(count*steps);
      Bbody_y[count] = (2.*pow(Bbody_x[count],2.))/
      	(h2*c2*(exp(Bbody_x[count]/(k*Temp)) - 1.));
      sum = sum + Bbody_y[count];
      BBHist[count+1] = BBHist[count] + Bbody_y[count];
      count++;
    }

  Bbody_x[10000] = Emax;
  // Normalise cumulative histo.
  count = 0;
  while(count<10001)
    {
      BBHist[count] = BBHist[count]/sum;
      count++;
    }
}

void G4GeneralParticleSource::InputEnergySpectra(G4bool value)
{
  // Allows user to specifiy spectrum is momentum
  EnergySpec = value; // false if momentum
  if(verbosityLevel == 2)
    G4cout << "EnergySpec has value " << EnergySpec << G4endl;
}

void G4GeneralParticleSource::InputDifferentialSpectra(G4bool value)
{
  // Allows user to specify integral or differential spectra
  DiffSpec = value; // true = differential, false = integral
  if(verbosityLevel == 2)
    G4cout << "Diffspec has value " << DiffSpec << G4endl;
}

void G4GeneralParticleSource::ArbInterpolate(G4String IType)
{
  if(EnergyDisType != "Arb")
    G4cout << "Error: this is for arbitrary distributions" << G4endl;
  IntType = IType;
  G4int i=0;

  // Calcuate Emin and Emax, mainly for use in debugging
  G4int len = G4int(ArbEnergyH.GetVectorLength());
  ArbEmax = ArbEnergyH.GetLowEdgeEnergy(size_t(len-1));
  ArbEmin = ArbEnergyH.GetLowEdgeEnergy(size_t(0));
  G4cout << "ArbEmin, ArbEmax " << ArbEmin << " " << ArbEmax << G4endl;
  
  // Now interpolate points
  if(IntType == "Lin")
    LinearInterpolation();
  if(IntType == "Log")
    LogInterpolation();
  if(IntType == "Exp")
    ExpInterpolation();
  if(IntType == "Spline")
    SplineInterpolation();  
}

void G4GeneralParticleSource::LinearInterpolation()
{
  // Method to do linear interpolation on the Arb points
  // Calculate equation of each line segment, max 256.
  // Calculate Area under each segment
  // Create a cumulative array which is then normalised Arb_Cum_Area
  G4double Area_seg[256]; // Stores area under each segment
  G4double sum = 0., Arb_x[256], Arb_y[256], Arb_Cum_Area[256];
  G4int i;
  G4int maxi = ArbEnergyH.GetVectorLength();
  for(i=0;i<maxi;i++)
    {
      Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
      Arb_y[i] = ArbEnergyH(size_t(i));
    }
  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy.
  if(EnergySpec == false)
    {
      // change currently stored values (emin etc) which are actually momenta
      // to energies.
      if(particle_definition == NULL)
	G4cout << "Error: particle not defined" << G4endl;
      else
	{
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent.
	  G4double mass = particle_definition->GetPDGMass();

	  // multiply the function (Arb_y) up by the bin width
	  // to make the function counts/s (i.e. get rid of momentum
	  // dependence).
	  for(int count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count] * (Arb_x[count+1] - Arb_x[count]);
	    }
	  // Change Arb_x to energy, plus divide by energy bin width
	  // to make Arb_y counts/s/energy
	  for(count=0;count<maxi;count++)
	    {
	      Arb_x[count] = sqrt((Arb_x[count]*Arb_x[count]) 
				  + (mass*mass));
	    }
	  for(count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count]/(Arb_x[count+1] - Arb_x[count]);
	    }
	}  
    }
  if(DiffSpec == false)
    {
      // Converts integral point-wise spectra to Differential
      for(int count=1;count<=maxi;count++)
	{
	  Arb_y[count] = Arb_y[count] - Arb_y[count-1];
	}
    }

  i=1;
  Arb_grad[0] = 0.;
  Arb_cept[0] = 0.;
  Area_seg[0] = 0.;
  Arb_Cum_Area[0] = 0.;
  while(i < maxi)
    {
      // calc gradient and intercept for each segment
      Arb_grad[i] = (Arb_y[i] - Arb_y[i-1]) / (Arb_x[i] - Arb_x[i-1]);
      if(verbosityLevel == 2)
	G4cout << Arb_grad[i] << G4endl;
      if(Arb_grad[i] > 0.)
	{
	  if(verbosityLevel == 2)
	    G4cout << "Arb_grad is positive" << G4endl;
	  Arb_cept[i] = Arb_y[i] - (Arb_grad[i] * Arb_x[i]);
	}
      else if(Arb_grad[i] < 0.)
	{
	  if(verbosityLevel == 2)
	    G4cout << "Arb_grad is negative" << G4endl;
	  Arb_cept[i] = Arb_y[i] + (-Arb_grad[i] * Arb_x[i]);
	}
      else
	{
	  if(verbosityLevel == 2)
	    G4cout << "Arb_grad is 0." << G4endl;
	  Arb_cept[i] = Arb_y[i];
	}

      Area_seg[i] = ((Arb_grad[i]/2)*(Arb_x[i]*Arb_x[i] - Arb_x[i-1]*Arb_x[i-1]) + Arb_cept[i]*(Arb_x[i] - Arb_x[i-1]));
      Arb_Cum_Area[i] = Arb_Cum_Area[i-1] + Area_seg[i];
      sum = sum + Area_seg[i];
      if(verbosityLevel == 2)
	G4cout << Arb_x[i] << Arb_y[i] << Area_seg[i] << sum << Arb_grad[i] << G4endl;
      i++;
    }

  i=0;
  while(i < maxi)
    {
      Arb_Cum_Area[i] = Arb_Cum_Area[i]/sum; // normalisation
      IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
      i++;
    }

  if(verbosityLevel >= 1)
    {
      G4cout << "Leaving LinearInterpolation" << G4endl;
      ArbEnergyH.DumpValues();
      IPDFArbEnergyH.DumpValues();
    }
}

void G4GeneralParticleSource::LogInterpolation()
{
  // Interpolation based on Logarithmic equations
  // Generate equations of line segments
  // y = Ax**alpha => log y = alpha*logx + logA
  // Find area under line segments
  // create normalised, cumulative array Arb_Cum_Area
  G4double Area_seg[256], Arb_x[256], Arb_y[256], Arb_Cum_Area[256];
  G4double alp, sum=0.;
  Arb_Cum_Area[0] = 0.;
  if(verbosityLevel == 2)
    ArbEnergyH.DumpValues();

  G4int i;
  G4int maxi = ArbEnergyH.GetVectorLength();
  for(i=0;i<maxi;i++)
    {
      Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
      Arb_y[i] = ArbEnergyH(size_t(i));
    }
  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy.
  if(EnergySpec == false)
    {
      // change currently stored values (emin etc) which are actually momenta
      // to energies.
      if(particle_definition == NULL)
	G4cout << "Error: particle not defined" << G4endl;
      else
	{
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent.
	  G4double mass = particle_definition->GetPDGMass();

	  // multiply the function (Arb_y) up by the bin width
	  // to make the function counts/s (i.e. get rid of momentum
	  // dependence).
	  for(int count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count] * (Arb_x[count+1] - Arb_x[count]);
	    }
	  // Change Arb_x to energy, plus divide by energy bin width
	  // to make Arb_y counts/s/energy
	  for(count=0;count<maxi;count++)
	    {
	      Arb_x[count] = sqrt((Arb_x[count]*Arb_x[count]) 
				  + (mass*mass));
	    }
	  for(count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count]/(Arb_x[count+1] - Arb_x[count]);
	    }
	}  
    }
  if(DiffSpec == false)
    {
      // Converts integral point-wise spectra to Differential
      for(int count=1;count<=maxi;count++)
	{
	  Arb_y[count] = Arb_y[count] - Arb_y[count-1];
	}
    }

  i=1;
  Arb_alpha[0] = 0.;
  Arb_Const[0] = 0.;
  Area_seg[0] = 0.;
  if(Arb_x[0] <= 0. || Arb_y[0] <= 0.)
    {
      G4cout << "You should not use log interpolation with points <= 0." << G4endl;
      G4cout << "These will be changed to 1e-20, which may cause problems" << G4endl;
      if(Arb_x[0] <= 0.)
	Arb_x[0] = 1e-20;
      if(Arb_y[0] <= 0.)
	Arb_y[0] = 1e-20;
    }
  
  while(i <maxi)
    {
      // Incase points are negative or zero
      if(Arb_x[i] <= 0. || Arb_y[i] <= 0.)
	{
	  G4cout << "You should not use log interpolation with points <= 0." << G4endl;
	  G4cout << "These will be changed to 1e-20, which may cause problems" << G4endl;
	  if(Arb_x[i] <= 0.)
	    Arb_x[i] = 1e-20;
	  if(Arb_y[i] <= 0.)
	    Arb_y[i] = 1e-20;
	}

      Arb_alpha[i] = (log10(Arb_y[i])-log10(Arb_y[i-1]))/(log10(Arb_x[i])-log10(Arb_x[i-1]));
      Arb_Const[i] = Arb_y[i]/(pow(Arb_x[i],Arb_alpha[i]));
      alp = Arb_alpha[i] + 1;
      Area_seg[i] = (Arb_Const[i]/alp) * (pow(Arb_x[i],alp) - pow(Arb_x[i-1],alp));
      sum = sum + Area_seg[i];
      Arb_Cum_Area[i] = Arb_Cum_Area[i-1] + Area_seg[i];
      if(verbosityLevel == 2)
	G4cout << Arb_alpha[i] << Arb_Const[i] << Area_seg[i] << G4endl;
      i++;
    }
  
  i=0;
  while(i<maxi)
    {
      Arb_Cum_Area[i] = Arb_Cum_Area[i]/sum;
      IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
      i++;
    }
  if(verbosityLevel >= 1)
    G4cout << "Leaving LogInterpolation " << G4endl;
}

void G4GeneralParticleSource::ExpInterpolation()
{
  // Interpolation based on Exponential equations
  // Generate equations of line segments
  // y = Ae**-(x/e0) => ln y = -x/e0 + lnA
  // Find area under line segments
  // create normalised, cumulative array Arb_Cum_Area
  G4double Area_seg[256], Arb_x[256], Arb_y[256], Arb_Cum_Area[256];
  G4double sum=0.;
  Arb_Cum_Area[0] = 0.;
  if(verbosityLevel == 2)
    ArbEnergyH.DumpValues();

  G4int i;
  G4int maxi = ArbEnergyH.GetVectorLength();
  for(i=0;i<maxi;i++)
    {
      Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
      Arb_y[i] = ArbEnergyH(size_t(i));
    }
  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy.
  if(EnergySpec == false)
    {
      // change currently stored values (emin etc) which are actually momenta
      // to energies.
      if(particle_definition == NULL)
	G4cout << "Error: particle not defined" << G4endl;
      else
	{
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent.
	  G4double mass = particle_definition->GetPDGMass();

	  // multiply the function (Arb_y) up by the bin width
	  // to make the function counts/s (i.e. get rid of momentum
	  // dependence).
	  for(int count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count] * (Arb_x[count+1] - Arb_x[count]);
	    }
	  // Change Arb_x to energy, plus divide by energy bin width
	  // to make Arb_y counts/s/energy
	  for(count=0;count<maxi;count++)
	    {
	      Arb_x[count] = sqrt((Arb_x[count]*Arb_x[count]) 
				  + (mass*mass));
	    }
	  for(count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count]/(Arb_x[count+1] - Arb_x[count]);
	    }
	}  
    }
  if(DiffSpec == false)
    {
      // Converts integral point-wise spectra to Differential
      for(int count=1;count<=maxi;count++)
	{
	  Arb_y[count] = Arb_y[count] - Arb_y[count-1];
	}
    }

  i=1;
  Arb_ezero[0] = 0.;
  Arb_Const[0] = 0.;
  Area_seg[0] = 0.;
  while(i < maxi)
    {
      G4double test = log(Arb_y[i]) - log(Arb_y[i-1]);
      if(test > 0. || test < 0.)
	{
	  Arb_ezero[i] = -(Arb_x[i] - Arb_x[i-1])/(log(Arb_y[i]) - log(Arb_y[i-1]));
	  Arb_Const[i] = Arb_y[i]/(exp(-Arb_x[i]/Arb_ezero[i]));
	  Area_seg[i]=-(Arb_Const[i]*Arb_ezero[i])*(exp(-Arb_x[i]/Arb_ezero[i])
						    -exp(-Arb_x[i-1]/Arb_ezero[i]));
	}
      else 
	{
	  G4cout << "Flat line segment: problem" << G4endl;
	  Arb_ezero[i] = 0.;
	  Arb_Const[i] = 0.;
	  Area_seg[i] = 0.;
	}
      sum = sum + Area_seg[i];
      Arb_Cum_Area[i] = Arb_Cum_Area[i-1] + Area_seg[i];
      if(verbosityLevel == 2)
	G4cout << Arb_ezero[i] << Arb_Const[i] << Area_seg[i] << G4endl;
      i++;
    }
  
  i=0;
  while(i<maxi)
    {
      Arb_Cum_Area[i] = Arb_Cum_Area[i]/sum;
      IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
      i++;
    }
  if(verbosityLevel >= 1)
    G4cout << "Leaving ExpInterpolation " << G4endl;
}

void G4GeneralParticleSource::SplineInterpolation()
{
  // Interpolation using Splines.
  // Create Normalised arrays, make x 0->1 and y hold
  // the function (Energy)
  G4double Arb_x[256], Arb_y[256];
  G4double sum = 0.;
  G4int i;
  if(verbosityLevel == 2)
    ArbEnergyH.DumpValues();

  G4int maxi = ArbEnergyH.GetVectorLength();
  for(i=0;i<maxi;i++)
    {
      Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
      //      if(i == 0)
      Arb_y[i] = ArbEnergyH(size_t(i));
	//      else
	//	Arb_y[i] = Arb_y[i-1] + ArbEnergyH(size_t(i));
	//      sum = sum + ArbEnergyH(size_t(i));
    }
  // Points are now in x,y arrays. If the spectrum is integral it has to be
  // made differential and if momentum it has to be made energy.
  if(EnergySpec == false)
    {
      // change currently stored values (emin etc) which are actually momenta
      // to energies.
      if(particle_definition == NULL)
	G4cout << "Error: particle not defined" << G4endl;
      else
	{
      // Apply Energy**2 = p**2c**2 + m0**2c**4
      // p should be entered as E/c i.e. without the division by c
      // being done - energy equivalent.
	  G4double mass = particle_definition->GetPDGMass();

	  // multiply the function (Arb_y) up by the bin width
	  // to make the function counts/s (i.e. get rid of momentum
	  // dependence).
	  for(int count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count] * (Arb_x[count+1] - Arb_x[count]);
	    }
	  // Change Arb_x to energy, plus divide by energy bin width
	  // to make Arb_y counts/s/energy
	  for(count=0;count<maxi;count++)
	    {
	      Arb_x[count] = sqrt((Arb_x[count]*Arb_x[count]) 
				  + (mass*mass));
	    }
	  for(count=0;count<maxi;count++)
	    {
	      Arb_y[count] = Arb_y[count]/(Arb_x[count+1] - Arb_x[count]);
	    }
	}  
    }
  if(DiffSpec == false)
    {
      // Converts integral point-wise spectra to Differential
      for(int count=1;count<=maxi;count++)
	{
	  Arb_y[count] = Arb_y[count] - Arb_y[count-1];
	}
    }

  for(i=1;i<maxi;i++)
    {
      //      Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
      //      if(i == 0)
      //	Arb_y[i] = ArbEnergyH(size_t(i));
      //      else
      Arb_y[i] = Arb_y[i-1] + ArbEnergyH(size_t(i));
      sum = sum + ArbEnergyH(size_t(i));
    }

  for(i=0;i<maxi;i++)
    {
      Arb_y[i] = Arb_y[i]/sum;
    }

  if((Arb_x[0] != 0.) || (Arb_y[0] != 0.))
    {
      for(i=maxi;i>=0;i--)
	{
	  if(i == 0)
	    {
	      Arb_x[0] = 0.;
	      Arb_y[0] = 0.;
	    }
	  else
	    {
	      Arb_x[i] = Arb_x[i-1];
	      Arb_y[i] = Arb_y[i-1];
	    }
	}
    }

  for(i=0;i<=maxi;i++)
    {
      if(verbosityLevel == 2)
	G4cout << i <<" "<< Arb_x[i] << " " << Arb_y[i] << G4endl;
      IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_y[i]);
    }

  // Should now have normalised cumulative probabilities in Arb_y
  // and energy values in Arb_x.
  maxi = maxi + 1;
  // Put y into x and x into y. The spline interpolation will then
  // go through x-axis to find where to interpolate (cum probability)
  // then generate a y (which will now be energy).
  SplineInt = new G4DataInterpolation(Arb_y,Arb_x,maxi,1e30,1e30);
  if(verbosityLevel == 2)
    {
      G4cout << SplineInt << G4endl;
      G4cout << SplineInt->LocateArgument(1.0) << G4endl;
    }
  if(verbosityLevel >= 1)
    G4cout << "Leaving SplineInterpolation " << G4endl;
}

void G4GeneralParticleSource::GenerateMonoEnergetic()
{
  // Method to generate MonoEnergetic particles.
  particle_energy = MonoEnergy;
}

void G4GeneralParticleSource::GenerateLinearEnergies()
{
  G4double energy, rndm;
  G4double emaxsq = pow(Emax,2.); //Emax squared
  G4double eminsq = pow(Emin,2.); //Emin squared
  G4double intersq = pow(cept,2.); //cept squared

  rndm = GenRandEnergy();
  
  G4double bracket = ((grad/2.)*(emaxsq - eminsq) + cept*(Emax-Emin));
  bracket = bracket * rndm;
  bracket = bracket + (grad/2.)*eminsq + cept*Emin;
  // Now have a quad of form m/2 E**2 + cE - bracket = 0
  bracket = -bracket;
  //  G4cout << "BRACKET" << bracket << G4endl;
  if(grad != 0.)
    {
      G4double sqbrack = (intersq - 4*(grad/2.)*(bracket));
      //      G4cout << "SQBRACK" << sqbrack << G4endl;
      sqbrack = sqrt(sqbrack);
      G4double root1 = -cept + sqbrack; 
      root1 = root1/(2.*(grad/2.));
      
      G4double root2 = -cept - sqbrack;
      root2 = root2/(2.*(grad/2.));

      //      G4cout << root1 << " roots " << root2 << G4endl;

      if(root1 > Emin && root1 < Emax)
	particle_energy = root1;
      if(root2 > Emin && root2 < Emax)
	particle_energy = root2;
    }
  else if(grad == 0.)
    // have equation of form cE - bracket =0
    particle_energy = bracket/cept;

  if(particle_energy < 0.)
    particle_energy = -particle_energy;
  
  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GeneratePowEnergies()
{
  // Method to generate particle energies distributed as
  // a powerlaw

  G4double rndm;
  G4double emina, emaxa;

  emina = pow(Emin,alpha+1);
  emaxa = pow(Emax,alpha+1);

  rndm = GenRandEnergy();

  if(alpha != -1.)
    {
      particle_energy = ((rndm*(emaxa - emina)) + emina);
      particle_energy = pow(particle_energy,(1./(alpha+1.)));
    }
  else if(alpha == -1.)
    {
      particle_energy = (log(Emin) + rndm*(log(Emax) - log(Emin)));
      particle_energy = exp(particle_energy);
    }
  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GenerateExpEnergies()
{
  // Method to generate particle energies distributed according
  // to an exponential curve.
  G4double rndm;
  rndm = GenRandEnergy();

  particle_energy = -Ezero*(log(rndm*(exp(-Emax/Ezero) - exp(-Emin/Ezero)) + 
				exp(-Emin/Ezero)));
  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GenerateBremEnergies()
{
  // Method to generate particle energies distributed according 
  // to a Bremstrahlung equation of 
  // form I = const*((kT)**1/2)*E*(e**(-E/kT))
  
  G4double rndm;
  rndm = GenRandEnergy();
  G4double expmax, expmin, k;

  k = 8.6181e-11; // Boltzmann's const in MeV/K
  G4double ksq = pow(k,2.); // k squared
  G4double Tsq = pow(Temp,2.); // Temp squared

  expmax = exp(-Emax/(k*Temp));
  expmin = exp(-Emin/(k*Temp));

  // If either expmax or expmin are zero then this will cause problems
  // Most probably this will be because T is too low or E is too high

  if(expmax == 0.)
    G4cout << "*****EXPMAX=0. Choose different E's or Temp" << G4endl;
  if(expmin == 0.)
    G4cout << "*****EXPMIN=0. Choose different E's or Temp" << G4endl;

  G4double tempvar = rndm *((-k)*Temp*(Emax*expmax - Emin*expmin) -
    (ksq*Tsq*(expmax-expmin)));

  G4double bigc = (tempvar - k*Temp*Emin*expmin - ksq*Tsq*expmin)/(-k*Temp);

  // This gives an equation of form: Ee(-E/kT) + kTe(-E/kT) - C =0
  // Solve this iteratively, step from Emin to Emax in 1000 steps
  // and take the best solution.

  G4double erange = Emax - Emin;
  G4double steps = erange/1000.;
  G4int i;
  G4double etest, diff, err;
  
  err = 100000.;

  for(i=1; i<1000; i++)
    {
      etest = Emin + (i-1)*steps;
      
      diff = etest*(exp(-etest/(k*Temp))) + k*Temp*(exp(-etest/(k*Temp))) - bigc;

      if(diff < 0.)
	diff = -diff;

      if(diff < err)
	{
	  err = diff;
	  particle_energy = etest;
	}
    }  
  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GenerateBbodyEnergies()
{
  // BBody_x holds Energies, and BBHist holds the cumulative histo.
  // binary search to find correct bin then lin interpolation.
  // Use the earlier defined histogram + RandGeneral method to generate
  // random numbers following the histos distribution.
  G4double rndm;
  G4int nabove, nbelow = 0, middle;
  nabove = 10001;
  rndm = GenRandEnergy();

  // Binary search to find bin that rndm is in
  while(nabove-nbelow > 1)
    {
      middle = (nabove + nbelow)/2;
      if(rndm == BBHist[middle]) break;
      if(rndm < BBHist[middle]) nabove = middle;
      else nbelow = middle;
    }
  
  // Now interpolate in that bin to find the correct output value.
  G4double x1, x2, y1, y2, m, q;
  x1 = Bbody_x[nbelow];
  x2 = Bbody_x[nbelow+1];
  y1 = BBHist[nbelow];
  y2 = BBHist[nbelow+1];
  m = (y2-y1)/(x2-x1);
  q = y1 - m*x1;
  
  particle_energy = (rndm - q)/m;

  if(verbosityLevel >= 1)
    {
      G4cout << "Energy is " << particle_energy << G4endl;
    }
}

void G4GeneralParticleSource::GenerateCdgEnergies()
{
  // Gen random numbers, compare with values in cumhist
  // to find appropriate part of spectrum and then 
  // generate energy in the usual inversion way.
  //  G4double pfact[2] = {8.5, 112};
  // G4double spind[2] = {1.4, 2.3};
  // G4double ene_line[3] = {1., 18., 1E6};
  G4double rndm, rndm2;
  G4double ene_line[3];
  G4double omalpha[2];
  if(Emin < 18*keV && Emax < 18*keV)
    {
      omalpha[0] = 1. - 1.4;
      ene_line[0] = Emin;
      ene_line[1] = Emax;
    }
  if(Emin < 18*keV && Emax > 18*keV)
    {
      omalpha[0] = 1. - 1.4;
      omalpha[1] = 1. - 2.3;
      ene_line[0] = Emin;
      ene_line[1] = 18.;
      ene_line[2] = Emax;
    }
  if(Emin > 18*keV)
    {
      omalpha[0] = 1. - 2.3;
      ene_line[0] = Emin;
      ene_line[1] = Emax;
    }
  rndm = GenRandEnergy();
  rndm2 = GenRandEnergy();

  G4int i = 0;
  while( rndm >= CDGhist[i])
    {
      i++;
    }
  // Generate final energy.
  particle_energy = (pow(ene_line[i-1],omalpha[i-1]) + (pow(ene_line[i],omalpha[i-1])
					- pow(ene_line[i-1],omalpha[i-1]))*rndm2);
  particle_energy = pow(particle_energy,(1./omalpha[i-1]));

  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GenUserHistEnergies()
{
  // Histograms are DIFFERENTIAL.
  //  G4cout << "In GenUserHistEnergies " << G4endl;
  if(IPDFEnergyExist == false)
    {
      G4int ii;
      G4int maxbin = G4int(UDefEnergyH.GetVectorLength());
      G4double bins[256], vals[256], sum;
      //      UDefEnergyH.DumpValues();
      G4double mass = particle_definition->GetPDGMass();
      //      G4cout << mass << G4endl;
      //      G4cout << EnergySpec << " " << DiffSpec << G4endl;
      if((EnergySpec == false) && (particle_definition == NULL))
	G4cout << "Error: particle definition is NULL" << G4endl;
      
      if(maxbin > 256)
	{
	  G4cout << "Maxbin > 256" << G4endl;
	  G4cout << "Setting maxbin to 256, other bins are lost" << G4endl;
	}

      if(DiffSpec == false)
	G4cout << "Histograms are Differential!!! " << G4endl;
      else
	{
	  //	  G4cout << "Here 2" << G4endl;
	  bins[0] = UDefEnergyH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = UDefEnergyH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = UDefEnergyH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = UDefEnergyH(size_t(ii)) + vals[ii-1];
	      sum = sum + UDefEnergyH(size_t(ii));
	    }
	}

      if(EnergySpec == false)
	{
	  // multiply the function (vals) up by the bin width
	  // to make the function counts/s (i.e. get rid of momentum
	  // dependence).
	  for(ii=1;ii<maxbin;ii++)
	    {
	      //	      G4cout << vals[ii] << " " << bins[ii] << " " << bins[ii-1] << G4endl;
	      vals[ii] = vals[ii] * (bins[ii] - bins[ii-1]);
	    }
	  // Put energy bins into new histo, plus divide by energy bin width
	  // to make evals counts/s/energy
	  for(ii=0;ii<maxbin;ii++)
	    {
	      bins[ii] = sqrt((bins[ii]*bins[ii]) + (mass*mass));
	      //	      G4cout << bins[ii] << " " << mass << G4endl;
	    }
	  for(ii=1;ii<maxbin;ii++)
	    {
	      //	      G4cout << vals[ii] << " " << bins[ii] << " " << bins[ii-1] << G4endl;
	      vals[ii] = vals[ii]/(bins[ii] - bins[ii-1]);
	    }
	  sum = vals[maxbin-1];
	}

      for(ii=0;ii<maxbin;ii++)
	{
	  vals[ii] = vals[ii]/sum;
	  IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
	}

      // Make IPDFEnergyExist = true
      IPDFEnergyExist = true;
      if(verbosityLevel == 2)
	IPDFEnergyH.DumpValues();
      
    }

  // IPDF has been create so carry on
  G4double rndm = GenRandEnergy();
  particle_energy = IPDFEnergyH.GetEnergy(rndm);
  
  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::GenArbPointEnergies()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenArbPointEnergies" << G4endl;
  G4double rndm;
  rndm = GenRandEnergy();

  if(IntType != "Spline")
    {
      //      IPDFArbEnergyH.DumpValues();
      // Find the Bin
      // have x, y, no of points, and cumulative area distribution
      G4int nabove, nbelow = 0, middle;
      nabove = IPDFArbEnergyH.GetVectorLength();
      //      G4cout << nabove << G4endl;
      // Binary search to find bin that rndm is in
      while(nabove-nbelow > 1)
	{
	  middle = (nabove + nbelow)/2;
	  if(rndm == IPDFArbEnergyH(size_t(middle))) break;
	  if(rndm < IPDFArbEnergyH(size_t(middle))) nabove = middle;
	  else nbelow = middle;
	}
      if(IntType == "Lin")
	{
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow+1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  grad = Arb_grad[nbelow+1];
	  cept = Arb_cept[nbelow+1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << grad << " " << cept << G4endl;
	  GenerateLinearEnergies();
	}
      else if(IntType == "Log")
	{
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow+1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  alpha = Arb_alpha[nbelow+1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << alpha << G4endl;
	  GeneratePowEnergies();
	}
      else if(IntType == "Exp")
	{
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow+1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  Ezero = Arb_ezero[nbelow+1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << Ezero << G4endl;
	  GenerateExpEnergies();
	}
    }
  else if(IntType == "Spline")
    {
      if(verbosityLevel == 2)
	G4cout << "IntType = Spline " << rndm << G4endl;
      // in SplineInterpolation created SplineInt
      // Now generate a random number put it into CubicSplineInterpolation
      // and you should get out an energy!?!
      particle_energy = SplineInt->CubicSplineInterpolation(rndm);
      if(verbosityLevel >= 1)
	G4cout << "Energy is " << particle_energy << G4endl;
    }
  else
    G4cout << "Error: IntType unknown type" << G4endl;
}

void G4GeneralParticleSource::GenEpnHistEnergies()
{
  //  G4cout << "In GenEpnHistEnergies " << Epnflag << G4endl;
  
  // Firstly convert to energy if not already done.
  if(Epnflag == true)
    // epnflag = true means spectrum is epn, false means e.
    {
      // convert to energy by multiplying by A number
      ConvertEPNToEnergy();
      // EpnEnergyH will be replace by UDefEnergyH.
      //      UDefEnergyH.DumpValues();
    }

  //  G4cout << "Creating IPDFEnergy if not already done so" << G4endl;
  if(IPDFEnergyExist == false)
    {
      // IPDF has not been created, so create it
      G4double bins[256],vals[256], sum;
      G4int ii;
      G4int maxbin = G4int(UDefEnergyH.GetVectorLength());
      bins[0] = UDefEnergyH.GetLowEdgeEnergy(size_t(0));
      vals[0] = UDefEnergyH(size_t(0));
      sum = vals[0];
      for(ii=1;ii<maxbin;ii++)
	{
	  bins[ii] = UDefEnergyH.GetLowEdgeEnergy(size_t(ii));
	  vals[ii] = UDefEnergyH(size_t(ii)) + vals[ii-1];
	  sum = sum + UDefEnergyH(size_t(ii));
	}
      
      for(ii=0;ii<maxbin;ii++)
	{
	  vals[ii] = vals[ii]/sum;
	  IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
	}
      // Make IPDFEpnExist = true
      IPDFEnergyExist = true;
    }
  //  IPDFEnergyH.DumpValues();
  // IPDF has been create so carry on
  G4double rndm = GenRandEnergy();
  particle_energy = IPDFEnergyH.GetEnergy(rndm);

  if(verbosityLevel >= 1)
    G4cout << "Energy is " << particle_energy << G4endl;
}

void G4GeneralParticleSource::ConvertEPNToEnergy()
{
  // Use this before particle generation to convert  the
  // currently stored histogram from energy/nucleon 
  // to energy.
  //  G4cout << "In ConvertEpntoEnergy " << G4endl;
  if(particle_definition==NULL)
    G4cout << "Error: particle not defined" << G4endl;
  else
    {
      // Need to multiply histogram by the number of nucleons.
      // Baryon Number looks to hold the no. of nucleons.
      G4int Bary = particle_definition->GetBaryonNumber();
      //      G4cout << "Baryon No. " << Bary << G4endl;
      // Change values in histogram, Read it out, delete it, re-create it
      G4int count, maxcount;
      maxcount = G4int(EpnEnergyH.GetVectorLength());
      //      G4cout << maxcount << G4endl;
      G4double ebins[256],evals[256];
      if(maxcount > 256)
	{
	  G4cout << "Histogram contains more than 256 bins!" << G4endl;
	  G4cout << "Those above 256 will be ignored" << G4endl;
	  maxcount = 256;
	}
      for(count=0;count<maxcount;count++)
	{
	  // Read out
	  ebins[count] = EpnEnergyH.GetLowEdgeEnergy(size_t(count));
	  evals[count] = EpnEnergyH(size_t(count));
	}

      // Multiply the channels by the nucleon number to give energies
      for(count=0;count<maxcount;count++)
	{
	  ebins[count] = ebins[count] * Bary;
	}

      // Set Emin and Emax
      Emin = ebins[0];
      Emax = ebins[maxcount-1];

      // Put energy bins into new histogram - UDefEnergyH.
      for(count=0;count<maxcount;count++)
	{
	  UDefEnergyH.InsertValues(ebins[count], evals[count]);
	}
      Epnflag = false; //so that you dont repeat this method.
    }  
}

// Biasing methods

void G4GeneralParticleSource::SetXBias(G4ThreeVector input)
{  
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  XBiasH.InsertValues(ehi, val);
  XBias = true;
}

void G4GeneralParticleSource::SetYBias(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  YBiasH.InsertValues(ehi, val);
  YBias = true;
}

void G4GeneralParticleSource::SetZBias(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  ZBiasH.InsertValues(ehi, val);
  ZBias = true;
}

void G4GeneralParticleSource::SetThetaBias(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  ThetaBiasH.InsertValues(ehi, val);
  ThetaBias = true;
}

void G4GeneralParticleSource::SetPhiBias(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  PhiBiasH.InsertValues(ehi, val);
  PhiBias = true;
}

void G4GeneralParticleSource::SetEnergyBias(G4ThreeVector input)
{
  G4double ehi, val;
  ehi = input.x();
  val = input.y();
  EnergyBiasH.InsertValues(ehi, val);
  EnergyBias = true;
}

G4double G4GeneralParticleSource::GenRandX()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenRandX" << G4endl;
  if(XBias == false)
    {
      // X is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // X is biased
      if(IPDFXBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(XBiasH.GetVectorLength());
	  bins[0] = XBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = XBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = XBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = XBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + XBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFXBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFXBias = true
	  IPDFXBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();

      // Calculate the weighting: Find the bin that the determined
      // rndm is in and the weigthing will be the difference in the
      // natural probability (from the x-axis) divided by the 
      // difference in the biased probability (the area).
      size_t numberOfBin = IPDFXBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFXBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      // retrieve the areas and then the x-axis values
      bweights[0] = IPDFXBiasH(biasn2) - IPDFXBiasH(biasn2 - 1);
      G4double xaxisl = IPDFXBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFXBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      //G4cout << "X Bin weight " << bweights[0] << " " << rndm << G4endl;
      //G4cout << "lower and upper xaxis vals "<<xaxisl<<" "<<xaxisu<<G4endl;
      bweights[0] = NatProb/bweights[0];
      if(verbosityLevel >= 1)
	G4cout << "X bin weight " << bweights[0] << " " << rndm << G4endl;
      return(IPDFXBiasH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenRandY()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenRandY" << G4endl;
  if(YBias == false)
    {
      // Y is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // Y is biased
      if(IPDFYBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(YBiasH.GetVectorLength());
	  bins[0] = YBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = YBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = YBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = YBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + YBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFYBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFYBias = true
	  IPDFYBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      size_t numberOfBin = IPDFYBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFYBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      bweights[1] =  IPDFYBiasH(biasn2) - IPDFYBiasH(biasn2 - 1);
      G4double xaxisl = IPDFYBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFYBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      bweights[1] = NatProb/bweights[1];
      if(verbosityLevel >= 1)
	G4cout << "Y bin weight " << bweights[1] << " " << rndm << G4endl;
      return(IPDFYBiasH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenRandZ()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenRandZ" << G4endl;
  if(ZBias == false)
    {
      // Z is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // Z is biased
      if(IPDFZBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(ZBiasH.GetVectorLength());
	  bins[0] = ZBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = ZBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = ZBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = ZBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + ZBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFZBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFZBias = true
	  IPDFZBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      //      size_t weight_bin_no = IPDFZBiasH.FindValueBinLocation(rndm);
      size_t numberOfBin = IPDFZBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFZBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      bweights[2] =  IPDFZBiasH(biasn2) - IPDFZBiasH(biasn2 - 1);
      G4double xaxisl = IPDFZBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFZBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      bweights[2] = NatProb/bweights[2];
      if(verbosityLevel >= 1)
	G4cout << "Z bin weight " << bweights[2] << " " << rndm << G4endl;
      return(IPDFZBiasH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenRandTheta()
{
  if(verbosityLevel >= 1)
    {
      G4cout << "In GenRandTheta" << G4endl;
      G4cout << "Verbosity " << verbosityLevel << G4endl;
    }
  if(ThetaBias == false)
    {
      // Theta is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // Theta is biased
      if(IPDFThetaBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(ThetaBiasH.GetVectorLength());
	  bins[0] = ThetaBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = ThetaBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = ThetaBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = ThetaBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + ThetaBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFThetaBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFThetaBias = true
	  IPDFThetaBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      //      size_t weight_bin_no = IPDFThetaBiasH.FindValueBinLocation(rndm);
      size_t numberOfBin = IPDFThetaBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFThetaBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      bweights[3] =  IPDFThetaBiasH(biasn2) - IPDFThetaBiasH(biasn2 - 1);
      G4double xaxisl = IPDFThetaBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFThetaBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      bweights[3] = NatProb/bweights[3];
      if(verbosityLevel >= 1)
	G4cout << "Theta bin weight " << bweights[3] << " " << rndm << G4endl;
      return(IPDFThetaBiasH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenRandPhi()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenRandPhi" << G4endl;
  if(PhiBias == false)
    {
      // Phi is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // Phi is biased
      if(IPDFPhiBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(PhiBiasH.GetVectorLength());
	  bins[0] = PhiBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = PhiBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = PhiBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = PhiBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + PhiBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFPhiBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFPhiBias = true
	  IPDFPhiBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      //      size_t weight_bin_no = IPDFPhiBiasH.FindValueBinLocation(rndm);
      size_t numberOfBin = IPDFPhiBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFPhiBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      bweights[4] =  IPDFPhiBiasH(biasn2) - IPDFPhiBiasH(biasn2 - 1);
      G4double xaxisl = IPDFPhiBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFPhiBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      bweights[4] = NatProb/bweights[4];
      if(verbosityLevel >= 1)
	G4cout << "Phi bin weight " << bweights[4] << " " << rndm << G4endl;
      return(IPDFPhiBiasH.GetEnergy(rndm));
    }
}

G4double G4GeneralParticleSource::GenRandEnergy()
{
  if(verbosityLevel >= 1)
    G4cout << "In GenRandEnergy" << G4endl;
  if(EnergyBias == false)
    {
      // Energy is not biased
      G4double rndm = G4UniformRand();
      return(rndm);
    }
  else
    {
      // ENERGY is biased
      if(IPDFEnergyBias == false)
	{
	  // IPDF has not been created, so create it
	  G4double bins[256],vals[256], sum;
	  G4int ii;
	  G4int maxbin = G4int(EnergyBiasH.GetVectorLength());
	  bins[0] = EnergyBiasH.GetLowEdgeEnergy(size_t(0));
	  vals[0] = EnergyBiasH(size_t(0));
	  sum = vals[0];
	  for(ii=1;ii<maxbin;ii++)
	    {
	      bins[ii] = EnergyBiasH.GetLowEdgeEnergy(size_t(ii));
	      vals[ii] = EnergyBiasH(size_t(ii)) + vals[ii-1];
	      sum = sum + EnergyBiasH(size_t(ii));
	    }

	  for(ii=0;ii<maxbin;ii++)
	    {
	      vals[ii] = vals[ii]/sum;
	      IPDFEnergyBiasH.InsertValues(bins[ii], vals[ii]);
	    }
	  // Make IPDFEnergyBias = true
	  IPDFEnergyBias = true;
	}
      // IPDF has been create so carry on
      G4double rndm = G4UniformRand();
      //  size_t weight_bin_no = IPDFEnergyBiasH.FindValueBinLocation(rndm);
      size_t numberOfBin = IPDFEnergyBiasH.GetVectorLength();
      G4int biasn1 = 0;
      G4int biasn2 = numberOfBin/2;
      G4int biasn3 = numberOfBin - 1;
      while (biasn1 != biasn3 - 1) {
      if (rndm > IPDFEnergyBiasH(biasn2))
         biasn1 = biasn2;
      else
         biasn3 = biasn2;
      biasn2 = biasn1 + (biasn3 - biasn1 + 1)/2;
      }
      bweights[5] =  IPDFEnergyBiasH(biasn2) - IPDFEnergyBiasH(biasn2 - 1);
      G4double xaxisl = IPDFEnergyBiasH.GetLowEdgeEnergy(size_t(biasn2-1));
      G4double xaxisu = IPDFEnergyBiasH.GetLowEdgeEnergy(size_t(biasn2));
      G4double NatProb = xaxisu - xaxisl;
      bweights[5] = NatProb/bweights[5];
      if(verbosityLevel >= 1)
	G4cout << "Energy bin weight " << bweights[5] << " " << rndm << G4endl;
      return(IPDFEnergyBiasH.GetEnergy(rndm));
    }
}

// Verbosity
void G4GeneralParticleSource::SetVerbosity(int vL)
{
  verbosityLevel = vL;
  G4cout << "Verbosity Set to: " << verbosityLevel << G4endl;
}

void G4GeneralParticleSource::SetParticleDefinition
  (G4ParticleDefinition* aParticleDefinition)
{
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}

// SR1.3
//void G4GeneralParticleSource::SetNucleus (Nucleus theIon1)
//{
//  theIon = theIon1;

//  G4IonTable *theIonTable = (G4IonTable *)(G4ParticleTable::GetParticleTable()->GetIonTable());
//  G4ParticleDefinition *aIon = NULL;

//  G4int A = theIon.GetA();
//  G4int Z = theIon.GetZ();
//  G4double E = theIon.GetE();

//  aIon = theIonTable->GetIon (Z, A, E);

//  SetParticleDefinition(aIon);
//}


void G4GeneralParticleSource::GeneratePrimaryVertex(G4Event *evt)
{
  if(particle_definition==NULL) return;

  // Position stuff
  G4bool srcconf = false;
  G4int LoopCount = 0;
  while(srcconf == false)
    {
      if(SourcePosType == "Point")
	GeneratePointSource();
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
  // Angular stuff
  if(AngDistType == "iso")
    GenerateIsotropicFlux();
  else if(AngDistType == "cos")
    GenerateCosineLawFlux();
  else if(AngDistType == "user")
    GenerateUserDefFlux();
  else
    G4cout << "Error: AngDistType has unusual value" << G4endl;
  // Energy stuff
  if(EnergyDisType == "Mono")
    GenerateMonoEnergetic();
  else if(EnergyDisType == "Lin")
    GenerateLinearEnergies();
  else if(EnergyDisType == "Pow")
    GeneratePowEnergies();
  else if(EnergyDisType == "Exp")
    GenerateExpEnergies();
  else if(EnergyDisType == "Brem")
    GenerateBremEnergies();
  else if(EnergyDisType == "Bbody")
    GenerateBbodyEnergies();
  else if(EnergyDisType == "Cdg")
    GenerateCdgEnergies();
  else if(EnergyDisType == "User")
    GenUserHistEnergies();
  else if(EnergyDisType == "Arb")
    GenArbPointEnergies();
  else if(EnergyDisType == "Epn")
    GenEpnHistEnergies();
  else
    G4cout << "Error: EnergyDisType has unusual value" << G4endl;

  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  if(verbosityLevel == 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;
  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();

  for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ )
  {
    G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
    particle->SetMass( mass );
    particle->SetCharge( particle_charge );
    particle->SetPolarization(particle_polarization.x(),
                               particle_polarization.y(),
                               particle_polarization.z());
    vertex->SetPrimary( particle );

    // Set bweight equal to the multiple of all non-zero weights
    bweight = 0.;
    for(int bib=0; bib<6; bib++)
      {
	if(bweights[bib] > 0. && bweight == 0.)
	  bweight = bweight + bweights[bib];
	else if(bweights[bib] > 0.)
	  bweight = bweight * bweights[bib];
      }
    // bweight will now contain the events final weighting.
    // now pass it to the primary vertex
    vertex->SetWeight(bweight);
  }

  evt->AddPrimaryVertex( vertex );
}






