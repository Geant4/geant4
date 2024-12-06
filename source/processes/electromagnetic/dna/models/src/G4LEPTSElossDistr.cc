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
// G4LEPTSElossDistr
//
// Author: Pedro Arce (CIEMAT), 2014
// --------------------------------------------------------------------

#include "G4LEPTSElossDistr.hh"
#include "G4LEPTSDistribution.hh"

#include <stdio.h>
#include <iostream>

G4LEPTSElossDistr::G4LEPTSElossDistr(const G4String& fname)
{
  fileName = fname;
  ReadFile();
}

void G4LEPTSElossDistr::ReadFile() 
{
  theNDistributions = 0;

  FILE * fp;

  if ((fp=fopen(fileName.c_str(), "r"))==nullptr)
  {
    NoBins = 0;
    bFileFound = false;
    return;
  } 

  bFileFound = true;

  G4int nEnergies;
  G4int nAngles;
  G4int nData;
  (void) fscanf(fp,"%i \n",&nEnergies);
  for( G4int ie = 0; ie < nEnergies; ++ie )
  {
    G4float energySep; 
    (void) fscanf(fp,"%f \n",&energySep);
    (void) fscanf(fp,"%i \n",&nAngles);
    for( G4int ia = 0; ia < nAngles; ++ia )
    {
      G4float angleSep; 
      (void) fscanf(fp,"%f \n",&angleSep);
      auto dist = new G4LEPTSDistribution();
      ++theNDistributions;
      std::map<G4double, G4LEPTSDistribution *> angleDist;
      angleDist[angleSep] = dist;
      theDistributions[energySep] = std::move(angleDist);
      
      (void) fscanf(fp,"%i \n",&nData);
      if( dist->ReadFile( fp, nData ) )
      {
        std::ostringstream message;
        message << "End of file found while reading file: "
                << fileName;
        G4Exception("G4LEPTSElossDistr::ReadFile()", "ReadError",
                    FatalException, message);	
      }
    }
  }
  
  fclose(fp);
}

G4double G4LEPTSElossDistr::Sample( G4double eMin, G4double eMax ) 
{
  // Sample Energy from Cumulative distr. G4interval [eMin, eMax]

  if( eMin > eMax) return 0.0;

  // Get the distribution to do the sampling
  G4LEPTSDistribution* distr  = nullptr;
  if( theNDistributions == 1 )
  {
    distr = (*( (*(theDistributions.cbegin())).second ).cbegin()).second;
  }
  else
  {
    for( auto itedd = theDistributions.cbegin(); itedd != theDistributions.cend(); ++itedd ){
      G4double energySep = (*itedd).first;
      if( eMax < energySep )
      {
        const auto& dist1 = (*itedd).second;
        for( auto ited = dist1.cbegin(); ited != dist1.cend(); ++ited )
        {
          G4double angleSep = (*ited).first;
	  if( 1 < angleSep )
          {
            distr = (*ited).second;
            break;
          }
        }
        break;
      }
    }
  }

  return (nullptr != distr) ? distr->Sample(eMin, eMax) : 0.0;
}
