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

#include "G4DefaultLinearColorMap.hh"

void G4DefaultLinearColorMap::GetMapColor(G4double val, G4double color[4])
{
  G4double value;
  if(fMaxVal == fMinVal)
    value = 0.;
  else
    value = (val - fMinVal) / (fMaxVal - fMinVal);

  if(value > 1.)
  {
    value = 1.;
  }
  if(value < 0.)
  {
    value = 0.;
  }

  // color map
  const int NCOLOR = 6;
  struct ColorMap
  {
    G4double val;
    G4double rgb[4];
  } colormap[NCOLOR] = { { 0.0, { 1., 1., 1., 1. } },  // value, r, g, b, alpha
                         { 0.2, { 0., 0., 1., 1. } },
                         { 0.4, { 0., 1., 1., 1. } },
                         { 0.6, { 0., 1., 0., 1. } },
                         { 0.8, { 1., 1., 0., 1. } },
                         { 1.0, { 1., 0., 0., 1. } } };

  // search
  G4int during[2] = { 0, 0 };
  for(int i = 1; i < NCOLOR; i++)
  {
    if(colormap[i].val >= value)
    {
      during[0] = i - 1;
      during[1] = i;
      break;
    }
  }

  // interpolate
  G4double a = std::fabs(value - colormap[during[0]].val);
  G4double b = std::fabs(value - colormap[during[1]].val);
  for(int i = 0; i < 4; i++)
  {
    color[i] =
      (b * colormap[during[0]].rgb[i] + a * colormap[during[1]].rgb[i]) /
      (colormap[during[1]].val - colormap[during[0]].val);
    if(color[i] > 1.)
      color[i] = 1.;
  }
}
