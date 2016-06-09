//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4LowEXsection.hh"
#include <cmath>

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
   G4double x1 = std::log((*it).first);
   G4double x2 = std::log((*(it+1)).first);
   G4double y1 = std::log((*it).second);
   G4double y2 = std::log((*(it+1)).second);
   G4double x = std::log(aX);
   G4double y = y1+(x-x1)*(y2-y1)/(x2-x1);
   result = std::exp(y);
   return result*millibarn;
 }
