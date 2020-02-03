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
///////////////////////////////////////////////////////////////////////////////
// File: CCaloOrganization.hh
// Description: Packing, unpacking and other related utilities for 
//              calorimetric numbering schema
///////////////////////////////////////////////////////////////////////////////

#ifndef CCaloOrganization_h
#define CCaloOrganization_h 1

#include <map>
#include "globals.hh"

class CCaloOrganization
{
public:
  CCaloOrganization(){}
  virtual ~CCaloOrganization(){}
	 
public:
  unsigned int packindex(G4int det, G4int z,
                         G4int eta, G4int phi) const;
  unsigned int packindex(G4int det, G4int depth, G4int z,
                         G4int eta, G4int phi) const;
  void unpackindex(const unsigned int& idx, G4int& det, G4int& z,
                   G4int& eta, G4int& phi) const;
  void unpackindex(const unsigned int& idx, G4int& det, G4int& depth, G4int& z,
                   G4int& eta, G4int& phi) const;
    
public:
      
  // additional tool:
  G4int getUnitWithMaxEnergy(std::map<G4int,G4float,std::less<G4int> >& themap);

  G4float  energyInMatrix(G4int nCellInEta, G4int nCellInPhi, 
			G4int crystalWithMaxEnergy, 
			std::map<G4int,G4float,std::less<G4int> >& themap);
};

#endif
