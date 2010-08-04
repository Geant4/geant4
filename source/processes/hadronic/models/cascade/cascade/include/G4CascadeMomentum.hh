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
// $Id: G4CascadeMomentum.hh,v 1.6 2010-08-04 05:28:24 mkelsey Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class G4CascadeMomentum
//
// Class description:
//
// A simple wrapper class meant to replace the widespread use of
// std::vector<double> in the cascade mode code, which causes
// problems for performance due to excess memory allocations.
//
// NOTE:  The Bertini code does not pass legitimate four-vectors when
//	  creating new particles; the new getLV() function takes an
//	  optional mass argument (in Bertini units [GeV]) so that a
//	  valid G4LorentzVector can be returned.
//
// Author: Peter Elmer, Princeton University                  7-Aug-2008
// Update: Michael Kelsey, SLAC (support G4LorentzVector)     8-Jan-2010
// --------------------------------------------------------------------
#ifndef G4CASCADE_MOMENTUM_HH
#define G4CASCADE_MOMENTUM_HH

#include <cassert>

#include "G4Types.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

class G4CascadeMomentum
{
  public:
    G4CascadeMomentum() {for (int i=0; i<4; ++i) data_[i]=0.0;}

    // WARNING!  This metric is (t,x,y,z), DIFFERENT FROM HepLV!
    G4CascadeMomentum(const G4LorentzVector& lv) {
      setLV(lv);
    }

    G4double& operator[](int i) {
      assert(i>=0 && i<4);
      return data_[i];
    }
    const G4double& operator[](int i) const {
      assert(i>=0 && i<4);
      return data_[i];
    }

    operator const G4LorentzVector&() const {
      return getLV();			// Casting can't do mass repairs
    }

    const G4LorentzVector& getLV(G4double mass=-1.) const {
      if (mass>=0.) lv.setVectM(get3V(),mass);		// Force input mass!
      else lv.set(data_[1],data_[2],data_[3],data_[0]);
      
      return lv;
    }

    G4ThreeVector get3V() const {
      return getLV().vect();
    }

    void setLV(const G4LorentzVector& lv) {
      data_[0] = lv.t();		// NOTE DIFFERENT METRIC CONVENTION!
      data_[1] = lv.x();
      data_[2] = lv.y();
      data_[3] = lv.z();
    }

    G4CascadeMomentum& operator=(const G4LorentzVector& lv) {
      setLV(lv);
      return *this;
    }

  private:
    G4double data_[4];
    mutable G4LorentzVector lv;		// Buffer for conversion operations
};

#endif

