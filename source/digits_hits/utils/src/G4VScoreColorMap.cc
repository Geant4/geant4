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
#include "G4Polyhedron.hh"

G4VScoreColorMap::G4VScoreColorMap(G4String mName)
  : fName(mName)
{}

void G4VScoreColorMap::DrawColorChart(G4int _nPoint)
{
  fVisManager = G4VVisManager::GetConcreteInstance();
  if(fVisManager == nullptr)
  {
    G4cerr << "G4VScoringMesh::DrawColorChart(): no visualization system"
           << G4endl;
    return;
  }

  DrawColorChartBar(_nPoint);
  DrawColorChartText(_nPoint);
}

void G4VScoreColorMap::DrawColorChartBar(G4int _nPoint)
{
  G4double min  = this->GetMin();
  G4double max  = this->GetMax();
  G4double smin = -0.89, smax = smin + 0.05 * (_nPoint) *0.83, step = 0.001;
  G4double c[4];

  fVisManager->BeginDraw2D();

  for(G4double y = smin; y < smax; y += step)
  {
    G4double ra = (y - smin) / (smax - smin), rb = 1. - ra;
    G4Polyline line;
    line.push_back(G4Point3D(-0.96, y, 0.));
    line.push_back(G4Point3D(-0.91, y, 0.));
    this->GetMapColor((ra * max + rb * min) / (ra + rb), c);
    G4Colour col(c[0], c[1], c[2]);
    G4VisAttributes att(col);
    line.SetVisAttributes(&att);
    fVisManager->Draw2D(line);
  }

  fVisManager->EndDraw2D();
}

void G4VScoreColorMap::DrawColorChartText(G4int _nPoint)
{
  G4double min = this->GetMin();
  G4double max = this->GetMax();
  G4double c[4];
  G4Colour black(0.1, 0.1, 0.1, 0.8);

  fVisManager->BeginDraw2D();

  for(int n = 0; n < _nPoint; n++)
  {
    G4double a = n / (_nPoint - 1.), b = 1. - a;
    G4double v = (a * max + b * min) / (a + b);
    // background color

    for(int l = 0; l < 21; l++)
    {
      G4Polyline line;
      line.push_back(G4Point3D(-0.9, -0.905 + 0.05 * n + 0.002 * l, 0.));
      line.push_back(G4Point3D(-0.75, -0.905 + 0.05 * n + 0.002 * l, 0.));
      G4VisAttributes attblack(black);
      line.SetVisAttributes(&attblack);
      fVisManager->Draw2D(line);
    }

    // text
    std::ostringstream oss;
    oss << std::setw(8) << std::setprecision(1) << std::scientific << v;
    G4String value  = oss.str();

    G4Text text(value, G4Point3D(-0.9, -0.9 + 0.05 * n, 0.4));
    G4double size = 12.;
    text.SetScreenSize(size);
    this->GetMapColor(v, c);
    G4Colour color(c[0], c[1], c[2], 1.);
    G4VisAttributes att(color);
    text.SetVisAttributes(&att);

    fVisManager->Draw2D(text);
  }

  // draw ps name
  // background
  G4double lpsname = 2. + fPSName.size() * 0.95;
  if(lpsname > 0)
  {
    for(int l = 0; l < 22; l++)
    {
      G4Polyline line;
      line.push_back(G4Point3D(-0.92, -0.965 + 0.002 * l, 0.));
      line.push_back(
        G4Point3D(-0.92 + 0.025 * lpsname, -0.965 + 0.002 * l, 0.));
      G4VisAttributes attblack(black);
      line.SetVisAttributes(&attblack);
      fVisManager->Draw2D(line);
    }

    // ps name
    G4Text txtpsname(fPSName, G4Point3D(-0.9, -0.96, 0.1));
    G4double size = 12.;
    txtpsname.SetScreenSize(size);
    G4Colour color(1., 1., 1.);
    G4VisAttributes att(color);
    txtpsname.SetVisAttributes(&att);
    fVisManager->Draw2D(txtpsname);
  }

  // draw unit
  // background
  G4double len = 2. + fPSUnit.size();
  if(len > 0)
  {
    for(int l = 0; l < 21; l++)
    {
      G4Polyline line;
      line.push_back(G4Point3D(-0.7, -0.9 + 0.002 * l, 0.));
      line.push_back(G4Point3D(-0.7 + 0.025 * len, -0.9 + 0.002 * l, 0.));
      G4VisAttributes attblack(black);
      line.SetVisAttributes(&attblack);
      fVisManager->Draw2D(line);
    }

    // unit
    G4String psunit = "[" + fPSUnit + "]";
    G4Text txtunit(psunit, G4Point3D(-0.69, -0.9, 0.1));
    G4double size = 12.;
    txtunit.SetScreenSize(size);
    G4Colour color(1., 1., 1.);
    G4VisAttributes att(color);
    txtunit.SetVisAttributes(&att);
    fVisManager->Draw2D(txtunit);
  }

  fVisManager->EndDraw2D();
}
