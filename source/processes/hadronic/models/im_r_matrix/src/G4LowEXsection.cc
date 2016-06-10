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
#include "G4LowEXsection.hh"
#include "G4SystemOfUnits.hh"

G4double G4LowEXsection::
 CrossSection(G4double aX) const
 {
   G4double result = 0;
   if(aX<front().first) return 0;
   G4LowEXsection::const_iterator i;
   G4LowEXsection::const_iterator it=end();
   for(i=begin(); i!=end(); i++)
   {
     if((*i).first/MeV>aX) break;
     it = i;
   }
   G4double x1 = G4Log((*it).first);
   G4double x2 = G4Log((*(it+1)).first);
   G4double y1 = G4Log((*it).second);
   G4double y2 = G4Log((*(it+1)).second);
   G4double x = G4Log(aX);
   G4double y = y1+(x-x1)*(y2-y1)/(x2-x1);
   result = G4Exp(y);
   return result*millibarn;
 }
