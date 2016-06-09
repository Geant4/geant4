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
// $Id: G4CascadeMomentum.hh,v 1.1 2008/09/22 10:06:32 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// Class G4CascadeMomentum
//
// Class description:
//
// A simple wrapper class meant to replace the widespread use of
// std::vector<double> in the cascade mode code, which causes
// problems for performance due to excess memory allocations.

// Author: Peter Elmer, Princeton University                  7-Aug-2008
// --------------------------------------------------------------------
#ifndef G4CASCADE_MOMENTUM_HH
#define G4CASCADE_MOMENTUM_HH

#include <cassert>

#include "G4Types.hh"

class G4CascadeMomentum
{
  public:

    G4CascadeMomentum() {for (int i=0; i<4; ++i) data_[i]=0.0;}

    G4double& operator[](int i)
    {
      assert(i>=0 && i<4);
      return data_[i];
    }
    const G4double& operator[](int i) const
    {
      assert(i>=0 && i<4);
      return data_[i];
    }

  private:

    G4double data_[4];
};
#endif

