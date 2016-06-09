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
// $Id: G4VScoreColorMap.cc,v 1.3.2.1 2009/08/11 13:47:53 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-04 $
//

#include "G4VScoreColorMap.hh"
#include <string>
#include <sstream>
#include <iomanip>

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"

G4VScoreColorMap::G4VScoreColorMap(G4String mName)
:fName(mName),ifFloat(true),fMinVal(0.),fMaxVal(DBL_MAX)
{;}

G4VScoreColorMap::~G4VScoreColorMap()
{;}

void G4VScoreColorMap::DrawColorChart(G4int _nPoint) {

  fVisManager = G4VVisManager::GetConcreteInstance();
  if(!fVisManager) {
    G4cerr << "G4VScoringMesh::DrawColorChart(): no visualization system" << G4endl;
    return;
  }

  DrawColorChartBar(_nPoint);
  DrawColorChartText(_nPoint);
}

void G4VScoreColorMap::DrawColorChartBar(G4int _nPoint) {

  G4double min = this->GetMin();
  G4double max = this->GetMax();
  G4double smin = -0.89, smax = smin + 0.05*(_nPoint)*0.83, step=0.001;
  G4double c[4];
  for(G4double y = smin; y < smax; y+=step) {
    G4double ra = (y-smin)/(smax-smin), rb = 1.-ra;
    G4Polyline line;
    line.push_back(G4Point3D(-0.96, y, 0.));
    line.push_back(G4Point3D(-0.91, y, 0.));
    this->GetMapColor((ra*max+rb*min)/(ra+rb), c);
    G4Colour col(c[0], c[1], c[2]);
    G4VisAttributes att(col);
    line.SetVisAttributes(&att);
    fVisManager->Draw2D(line);
  }

}
void G4VScoreColorMap::DrawColorChartText(G4int _nPoint) {
  G4double min = this->GetMin();
  G4double max = this->GetMax();
  G4double c[4];
  G4Colour black(0., 0., 0.);
  for(int n = 0; n < _nPoint; n++) {
    G4double a = n/(_nPoint-1.), b = 1.-a;
    G4double v = (a*max + b*min)/(a+b);
    // background color
    for(int l = 0; l < 21; l++) {
      G4Polyline line;
      line.push_back(G4Point3D(-0.908, -0.905+0.05*n+0.002*l, 0.));
      line.push_back(G4Point3D(-0.705, -0.905+0.05*n+0.002*l, 0.));
      G4VisAttributes attblack(black);
      line.SetVisAttributes(&attblack);
      fVisManager->Draw2D(line);
    }
    // text
    //char cstring[80]; 
    //std::sprintf(cstring, "%8.2e", v);
    //G4String value(cstring);
    std::ostringstream oss;
    oss << std::setw(8) << std::setprecision(1) << std::scientific << v;
    std::string str = oss.str();
    G4String value(str.c_str());

    G4Text text(value, G4Point3D(-0.9, -0.9+0.05*n, 0));
    G4double size = 12.;
    text.SetScreenSize(size);
    this->GetMapColor(v, c);
    G4Colour color(c[0], c[1], c[2]);
    G4VisAttributes att(color);
    text.SetVisAttributes(&att);

    fVisManager->Draw2D(text);
  }
}
