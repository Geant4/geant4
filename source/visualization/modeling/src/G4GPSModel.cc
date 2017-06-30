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
//
// $Id: G4GPSModel.cc 100959 2016-11-03 14:13:37Z allison $
//
//
// John Allison  26th April 2017.
// Model for a representation of the General Paricle Source.

#include "G4GPSModel.hh"

#include "G4VGraphicsScene.hh"
#include "G4GeneralParticleSourceData.hh"
#include "G4VisAttributes.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4Transform3D.hh"
#include "G4Point3D.hh"
#include "G4Circle.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4Para.hh"

#include <sstream>

G4GPSModel::G4GPSModel (const G4Colour& colour)
: fColour(colour)
{
  fType = "G4GPSModel";
  std::ostringstream oss;
  oss << "G4GPSModel for General Particle Source " << fColour;
  fGlobalTag = oss.str();
  fGlobalDescription = fGlobalTag;
}

G4GPSModel::~G4GPSModel () {}

G4String G4GPSModel::GetCurrentTag () const
{
return "";
}

G4String G4GPSModel::GetCurrentDescription () const
{
  return "G4GPSModel " + GetCurrentTag ();
}

void G4GPSModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
// The main task of a model is to describe itself to the graphics scene
// handler (a object which inherits G4VSceneHandler, which inherits
// G4VGraphicsScene).
{
  G4GeneralParticleSourceData* pGPSData = G4GeneralParticleSourceData::Instance();
  // Note: As far as I can see, if this is the first time Instance has been
  // called, it will, nevertheless, instantiate a default source, Type:"Point",
  // Shaep: "NULL", which will be drawn as a small circle at the origin of
  // coordinates whether you have set up GPS or not.  Sorry, can't think of a
  // way to avoid that.  Mostly, of course, you will only invoke this function,
  // if you have - or are about to - set up GPS, in which case all will be well.
  if (!pGPSData) return;

  G4int nSources = pGPSData->GetSourceVectorSize();
  for (G4int iSource = 0; iSource < nSources; ++iSource) {

    const G4SingleParticleSource* pCurrentSingleSource = pGPSData->GetCurrentSource(iSource);
    if (!pCurrentSingleSource) return;

    const G4SPSPosDistribution* pSPSPosDistribution = pCurrentSingleSource->GetPosDist();
    if (!pSPSPosDistribution) return;

    G4String Type = pSPSPosDistribution->GetPosDisType();
    G4String Shape = pSPSPosDistribution->GetPosDisShape();
    // Type can be: Point, Plane, Surface or Volume
    // Shape can be: Square, Circle, Ellipse, Rectangle,
    //    Sphere, Ellipsoid, Cylinder, Parallelepiped
//    G4cout
//    << "G4GPSModel::DescribeYourselfTo"
//    << ": PosDisType: " << Type
//    << ", Shape: " << Shape
//    << G4endl;

    const G4double& halfx = pSPSPosDistribution->GetHalfX();
    const G4double& halfy = pSPSPosDistribution->GetHalfY();
    const G4double& halfz = pSPSPosDistribution->GetHalfZ();
    const G4double& Radius  = pSPSPosDistribution->GetRadius();
    const G4double& Radius0 = pSPSPosDistribution->GetRadius0();
    const G4double& ParAlpha = pSPSPosDistribution->GetParAlpha();
    const G4double& ParTheta = pSPSPosDistribution->GetParTheta();
    const G4double& ParPhi   = pSPSPosDistribution->GetParPhi();

    const G4ThreeVector& Rotx = pSPSPosDistribution->GetRotx();
    const G4ThreeVector& Roty = pSPSPosDistribution->GetRoty();
    const G4ThreeVector& Rotz = pSPSPosDistribution->GetRotz();

    const G4ThreeVector& position = pSPSPosDistribution->GetCentreCoords();
    G4Transform3D transform(CLHEP::HepXHat,CLHEP::HepYHat,CLHEP::HepZHat,Rotx,Roty,Rotz);
    transform = G4Translate3D(position) * transform;

    G4double surfaceTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    G4double smallHalfThickness = 10.*surfaceTolerance;

    G4VisAttributes gpsAtts;
    gpsAtts.SetColour(fColour);
    gpsAtts.SetForceSolid();
    
    if (Type == "Point") {

      G4Circle circle;
      circle.SetPosition(position);
      circle.SetScreenDiameter(10.);
      circle.SetVisAttributes(gpsAtts);
      sceneHandler.BeginPrimitives(transform);
      sceneHandler.AddPrimitive(circle);
      sceneHandler.EndPrimitives();

    } else if (Type == "Plane") {

      // Code based on G4SPSPosDistribution::GeneratePointsInPlane.
      sceneHandler.PreAddSolid(transform,gpsAtts);
      if (Shape == "Circle") {
        sceneHandler.AddSolid
        (G4Tubs("GPS_Circle",0.,Radius,smallHalfThickness,0.,twopi));
      } else if (Shape == "Annulus") {
        sceneHandler.AddSolid
        (G4Tubs("GPS_Annulus",Radius0,Radius,smallHalfThickness,0.,twopi));
      } else if (Shape == "Ellipse") {
        sceneHandler.AddSolid
        (G4EllipticalTube("GPS_Ellipse",halfx,halfy,smallHalfThickness));
      } else if (Shape == "Square") {
        sceneHandler.AddSolid
        (G4Box("GPS_Ellipse",halfx,halfx,smallHalfThickness));
      } else if (Shape == "Rectangle") {
        sceneHandler.AddSolid
        (G4Box("GPS_Rectangle",halfx,halfy,smallHalfThickness));
      }
      sceneHandler.PostAddSolid();

    } else if (Type == "Surface" || Type == "Volume") {
      
      // Code based on G4SPSPosDistribution::GeneratePointsOnSurface.
      // and           G4SPSPosDistribution::GeneratePointsInVolume.
      sceneHandler.PreAddSolid(transform,gpsAtts);
      if (Shape == "Sphere") {
        sceneHandler.AddSolid
        (G4Orb("GPS_Sphere",Radius));
      } else if (Shape == "Ellipsoid") {
        sceneHandler.AddSolid
        (G4Ellipsoid("GPS_Ellipsoid",halfx,halfy,halfz));
      } else if (Shape == "Cylinder") {
        sceneHandler.AddSolid
        (G4Tubs("GPS_Cylinder",0.,Radius, halfz, 0., twopi));
      } else if (Shape == "Para") {
        sceneHandler.AddSolid
        (G4Para("GPS_Para",halfx,halfy,halfz,ParAlpha,ParTheta,ParPhi));
      }
      sceneHandler.PostAddSolid();

    }
  }
}
