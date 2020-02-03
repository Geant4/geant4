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
/// \file runAndEvent/RE02/include/RE02PSCellFlux.hh
/// \brief Definition of the RE02PSCellFlux class
//
//
//

#ifndef RE02PSCellFlux_h
#define RE02PSCellFlux_h 1

#include "G4PSCellFlux.hh"
///////////////////////////////////////////////////////////////////////////////
//
/// Cell flux class for a parameterized volume in a three dimentional structure
///
/// (Description)
///   This is a primitive scorer class for scoring cell flux.
///   The Cell Flux is defined by a sum of track length divided by the geometry
///   volume, where all of the tracks in the geometry are taken into account.
///   e.g. the unit of Cell Flux is mm/mm3.
///
///   If you score only tracks passing through the geometry volume, use
///   G4PSPassageCellFlux.
///
///
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

class RE02PSCellFlux : public G4PSCellFlux
{
   public: // with description
      RE02PSCellFlux(G4String name,G4int nx,G4int ny, G4int nz);
      virtual ~RE02PSCellFlux();

  protected: // with description
      virtual G4int GetIndex(G4Step*);

  private:
      G4int fNx, fNy, fNz;
};
#endif

