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
#include <cmath>
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4XpimNTotal.hh"
#include "G4SystemOfUnits.hh"

G4XpimNTotal::
 G4XpimNTotal()
 {
   std::pair<double,double> it;
   it.first=1105.46 ;it.second= 8.8;  theLowEData.push_back(it);
   it.first=1139.41 ;it.second= 13.6;  theLowEData.push_back(it);
   it.first=1165.12 ;it.second= 26.0;  theLowEData.push_back(it);
   it.first=1198.8  ;it.second= 60.;  theLowEData.push_back(it);
   it.first=1212.51 ;it.second= 70.;  theLowEData.push_back(it);
   it.first=1222.84 ;it.second= 72.;  theLowEData.push_back(it);
   it.first=1233.18 ;it.second= 68.;  theLowEData.push_back(it);
   it.first=1240.08 ;it.second= 65.;  theLowEData.push_back(it);
   it.first=1301.86 ;it.second= 30.;  theLowEData.push_back(it);
   it.first=1335.65 ;it.second= 26.;  theLowEData.push_back(it);
   it.first=1368.93 ;it.second= 27.;  theLowEData.push_back(it);
   it.first=1433.81 ;it.second= 30.;  theLowEData.push_back(it);
   it.first=1496.42 ;it.second= 45.;  theLowEData.push_back(it);
   it.first=1508.67 ;it.second= 47.;  theLowEData.push_back(it);
   it.first=1520.84 ;it.second= 45.;  theLowEData.push_back(it);
   it.first=1568.67 ;it.second= 35.;  theLowEData.push_back(it);
   it.first=1603.69 ;it.second= 37.;  theLowEData.push_back(it);
   it.first=1660.54 ;it.second= 58.;  theLowEData.push_back(it);
   it.first=1674.47 ;it.second= 59.;  theLowEData.push_back(it);
   it.first=1688.29 ;it.second= 58.;  theLowEData.push_back(it);
   it.first=1779.57 ;it.second= 36.;  theLowEData.push_back(it);
   it.first=1881.49 ;it.second= 36.8;  theLowEData.push_back(it);
   it.first=1978.31 ;it.second= 34.5;  theLowEData.push_back(it);
   it.first=2038.82 ;it.second= 34.5;  theLowEData.push_back(it);
   it.first=2115.39 ;it.second= 36.;  theLowEData.push_back(it);
   it.first=2159.18 ;it.second= 36.3;  theLowEData.push_back(it);
   it.first=2244.22 ;it.second= 36.;  theLowEData.push_back(it);
   it.first=2424.78 ;it.second= 33.;  theLowEData.push_back(it);
   it.first=2664.2  ;it.second= 32.;  theLowEData.push_back(it);
   it.first=3487.43 ;it.second= 28.;  theLowEData.push_back(it);
   it.first=4434.76 ;it.second= 26.7;  theLowEData.push_back(it);
 } 

G4double G4XpimNTotal::
 CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
 {
   G4double sqrts = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
   if(sqrts > theLowEData.back().first/MeV) return thePDGData.CrossSection(trk1, trk2);
   G4double result = 0;
   size_t i(0), it(0);
   if(sqrts<theLowEData[0].first) return 0;
   for(i=0; i<theLowEData.size(); i++)
   {
     if(theLowEData[i].first/MeV>sqrts) break;
     it = i;
   }
   G4double x1 = G4Log(theLowEData[it].first);
   G4double x2 = G4Log(theLowEData[it+1].first);
   G4double y1 = G4Log(theLowEData[it].second);
   G4double y2 = G4Log(theLowEData[it+1].second);
   G4double x = G4Log(sqrts);
   G4double y = y1+(x-x1)*(y2-y1)/(x2-x1);
   result = G4Exp(y);
   return result*millibarn;
 }
