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
/// \file runAndEvent/RE02/include/RE02PSEnergyDeposit.hh
/// \brief Definition of the RE02PSEnergyDeposit class
//
//
//

#ifndef RE02PSEnergyDeposit_h
#define RE02PSEnergyDeposit_h 1

#include "G4PSEnergyDeposit.hh"
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring energy deposit.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura 
///////////////////////////////////////////////////////////////////////////////

class RE02PSEnergyDeposit : public G4PSEnergyDeposit
{
   public: // with description
      RE02PSEnergyDeposit(G4String name,G4int nx,G4int ny, G4int nz);
      virtual ~RE02PSEnergyDeposit();

  protected: // with description
      virtual G4int GetIndex(G4Step*);

  private:
      G4int fNx, fNy, fNz;
};
#endif

