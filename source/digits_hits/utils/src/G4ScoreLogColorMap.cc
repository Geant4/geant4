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
// $Id: G4ScoreLogColorMap.cc,v 1.2 2008/03/25 04:06:08 akimura Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "G4ScoreLogColorMap.hh"
#include <cmath>


#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"

G4ScoreLogColorMap::G4ScoreLogColorMap(G4String mName)
:G4VScoreColorMap(mName)
{;}

G4ScoreLogColorMap::~G4ScoreLogColorMap()
{;}

void G4ScoreLogColorMap::GetMapColor(G4double val, G4double color[4])
{
  G4bool lmin = true, lmax = true, lval = true;
  if(fMinVal <= 0.) lmin = false;
  if(fMaxVal <= 0.) lmax = false;
  if(val <= 0.) lval = false;
  G4double logmin = 0., logmax = 0., logval = 0.;
  if(lmin) logmin = std::log10(fMinVal);
  if(lmax) logmax = std::log10(fMaxVal);
  if(lval) logval = std::log10(val);
  G4double value = 0.;
  if(lmax) value = (logval-logmin)/(logmax-logmin);

  if(value > 1.) {value=1.;}
  if(value < 0.) {value=0.;}

  // color map
  const int NCOLOR = 6;
  struct ColorMap {
    G4double val;
    G4double rgb[4];
  } colormap[NCOLOR] = {{0.0, {1., 1., 1., 1.}}, // value, r, g, b, alpha
			{0.2, {0., 0., 1., 1.}},
			{0.4, {0., 1., 1., 1.}},
			{0.6, {0., 1., 0., 1.}},
			{0.8, {1., 1., 0., 1.}},
			{1.0, {1., 0., 0., 1.}}};
  
  // search
  G4int during[2] = {0, 0};
  for(int i = 1; i < NCOLOR; i++) {
    if(colormap[i].val >= value) {
      during[0] = i-1;
      during[1] = i;
      break;
    }
  }

  // interpolate
  G4double a = std::fabs(value - colormap[during[0]].val);
  G4double b = std::fabs(value - colormap[during[1]].val);
  for(int i = 0; i < 4; i++) {
    color[i] = (b*colormap[during[0]].rgb[i] + a*colormap[during[1]].rgb[i])
      /(colormap[during[1]].val - colormap[during[0]].val);
  } 

}


void G4ScoreLogColorMap::DrawColorChartBar(G4int _nPoint) {

  G4bool lmin = true, lmax = true;
  if(fMinVal <= 0.) lmin = false;
  if(fMaxVal <= 0.) lmax = false;
  G4double min = 0.;
  if(lmin) min = std::log10(fMinVal);
  G4double max = 0.;
  if(lmax) max = std::log10(fMaxVal);

  G4double smin = -0.89, smax = smin + 0.05*(_nPoint)*0.83, step=0.001;
  G4double c[4];
  for(G4double y = smin; y < smax; y+=step) {
    G4double ra = (y-smin)/(smax-smin), rb = 1.-ra;
    G4Polyline line;
    line.push_back(G4Point3D(-0.96, y, 0.));
    line.push_back(G4Point3D(-0.91, y, 0.));
    G4double val = std::pow(10., (ra*max+rb*min)/(ra+rb));
    this->GetMapColor(val, c);
    G4Colour col(c[0], c[1], c[2]);
    G4VisAttributes att(col);
    line.SetVisAttributes(&att);
    fVisManager->Draw2D(line);
  }

}
void G4ScoreLogColorMap::DrawColorChartText(G4int _nPoint) {
  G4bool lmin = true, lmax = true;
  if(fMinVal <= 0.) lmin = false;
  if(fMaxVal <= 0.) lmax = false;

  G4double min = 0.;
  if(lmin) min = std::log10(fMinVal);
  if(min > 0.) min = std::floor(min);
  else min = std::ceil(min);

  G4double max = 0.;
  if(lmax) max = std::log10(fMaxVal);
  if(max > 0.) max = std::floor(max);
  else max = std::ceil(max);

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
    char cstring[80]; 
    std::sprintf(cstring, "%8.1e", std::pow(10., v));
    G4String value(cstring);
    G4Text text(value, G4Point3D(-0.9, -0.9+0.05*n, 0));
    G4double size = 12.;
    text.SetScreenSize(size);
    this->GetMapColor(std::pow(10., v), c);
    G4Colour color(c[0], c[1], c[2]);
    G4VisAttributes att(color);
    text.SetVisAttributes(&att);

    fVisManager->Draw2D(text);
  }
}
