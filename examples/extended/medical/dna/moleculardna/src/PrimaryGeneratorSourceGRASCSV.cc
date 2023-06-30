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
/// \file PrimaryGeneratorSourceGRASCSV.cc
/// \brief Primary Generator source class implementation for GRAS CSV phase space for Molecular DNA simulation

#include "PrimaryGeneratorSourceGRASCSV.hh"


PrimaryGeneratorSourceGRASCSV::PrimaryGeneratorSourceGRASCSV(const G4String& filename)
{
  // set a default buffer size
  fEndOfFile = false;

  // open file
  fInputFile.open(filename);
  if (fInputFile.is_open())
  {
    G4cout << "*** Opening particle source file " << filename << " ***" << G4endl;
    // find the starting particle source box
    G4String line = "";
    G4int totalevents = 0;
    while(std::getline(fInputFile, line))
    {
      if (line.find("TOTAL_EVENTS") != std::string::npos) {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        G4String tmp = "";
        ss >> tmp >> tmp >> totalevents;
      }
      if (line.find("TWO-STAGE PARTICLE PHASE SPACE OUTPUT FILE: PARTICLE PHASE-SPACE INFORMATION") != std::string::npos) {
        for(G4int i=0; i<18; i++)
       {
          std::getline(fInputFile, line);
          // Print the header information
          G4cout << line << G4endl;
       }
       break;
      }
    }
    // set the number of particles
    fnParticles = totalevents;
  }
  else
  {
    G4cout << "*** Warning: can't open particle source file in reading mode ***" << G4endl;
  }
}

PrimaryGeneratorSourceGRASCSV::~PrimaryGeneratorSourceGRASCSV()
{
  // close file and clear structures
  if(fInputFile.is_open())
  {
    fInputFile.close();
  }
  fPrimaryList.clear();
}

G4double PrimaryGeneratorSourceGRASCSV::RecomputeNParticles()
{
  // to be implemented
  // the variable in CSV has to be recomputed reading all the file
  return fnParticles;
}


Primary* PrimaryGeneratorSourceGRASCSV::GetPrimary()
{
  if(fInputFile.peek() == EOF || fEndOfFile == true)
    return nullptr;
  // Read in primaries if the primary list is empty
  if( fPrimaryList.size() == 0 )
  {
    // Read 100 primaries at a time
    G4int num = fBufferSize;
    if (num > fnParticles)
    {
      num = fnParticles;
    }

    for(G4int i = 0; i < num; i++ )
    {
      if(fInputFile.peek() == EOF)
      {
        // End of file
        fEndOfFile = true;
        break;
      }
      G4String line = "";
      std::getline(fInputFile, line);
      if (line.find("TWO-STAGE PARTICLE PHASE SPACE OUTPUT FILE: PARTICLE PHASE-SPACE INFORMATION") != std::string::npos) 
      {
        const G4int lines = 18;
        for(G4int j=0; j<lines; j++)
        {
          std::getline(fInputFile, line);
        }
        i--;
        continue;
      }
      if (line.find("Block") != std::string::npos || line.find("File") != std::string::npos || line.find("End") != std::string::npos || line.find("*") != std::string::npos) 
      {
        // End of Block, skip
        i--;
        continue;
      }
      std::stringstream ss(line);
      std::vector<G4String> tempVect;
      G4String tempString = "";
      while(ss>>tempString) tempVect.push_back(tempString);

      try
     {
      // Define position and momentum direction vectors
      G4ThreeVector pos = G4ThreeVector( std::stod(tempVect.at(10)), std::stod(tempVect.at(11)), std::stod(tempVect.at(12)));
      G4ThreeVector momDir = G4ThreeVector( std::stod(tempVect.at(16)), std::stod(tempVect.at(17)), std::stod(tempVect.at(18)));
      G4int PDGEncoding = std::stoi(tempVect.at(2));
      G4double KE = 0;
      KE = std::stod(tempVect.at(8));

      // Generate primary particle
      auto *primary = new Primary();
      primary->SetName( PDGEncoding );
      primary->SetPosition( pos );
      primary->SetMomentumDirection( momDir );
      primary->SetEnergy( KE );

      // Populate the primaries list
      fPrimaryList.push_back( primary );
      }
      catch (const std::invalid_argument&)
      {
        G4cerr << "Phase Space reading error: invalid_argument"<< G4endl;
      }
      catch (const std::out_of_range&)
      {
        G4cerr << "Phase Space reading error: out_of_range"<< G4endl;
      }
    }
  }

  // Get first element and delete it so it's not reused
  Primary* primary = nullptr;
  if( fPrimaryList.size() > 0 )
  {
    primary = fPrimaryList.front();
    fPrimaryList.pop_front();
  }
  return primary;
}



