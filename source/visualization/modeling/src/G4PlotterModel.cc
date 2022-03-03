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
// Guy Barrand 25th September 2021

#include "G4PlotterModel.hh"

#include "G4VGraphicsScene.hh"

G4PlotterModel::G4PlotterModel(G4Plotter& a_plotter,const G4String& description,const G4Transform3D& a_transform)
:fPlotter(a_plotter),fTransform(a_transform)
{
  fType = "G4PlotterModel";
  fGlobalTag = fType;
  fGlobalDescription = fType + ": " + description;
  // Have the extent so that the viewer will have an ortho camera
  // with a height of 1 in order to map the tools::sg::plots height of 1:
  //   radius = sqrt(half_x*half_x+half_y*half_y+half_z*half_z);
  //   radius = sqrt(3*half*half)
  //   half = sqrt(radius*radius/3)
  double radius = 0.5;
  double half = ::sqrt(radius*radius/3.0);
  fExtent = G4VisExtent(-half,half,-half,half,-half,half); //x,y,z min/max.
}

void G4PlotterModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  sceneHandler.BeginPrimitives(fTransform);
  sceneHandler.AddPrimitive(fPlotter);
  sceneHandler.EndPrimitives();
}
