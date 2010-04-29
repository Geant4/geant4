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
// $Id: G4FinalStateSampler.hh,v 1.6 2010-04-29 00:30:02 mkelsey Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: D. H. Wright
// Date:   26 March 2009
//
// 20100405  M. Kelsey -- Pass const-ref std::vector<>, improve base interface
// 20100413  M. Kelsey -- Move subclass functionality here
// 20100428  M. Kelsey -- Use G4InuclParticleNames enums in .cc file

#ifndef G4FinalStateSampler_h
#define G4FinalStateSampler_h 1

// Class Description:
//
// Implementation base class for sampling partial cross sections and final
// states for inelastic hadron-nucleon interactions

#include "globals.hh"
#include <vector>
#include "G4FastVector.hh"
#include "G4ReactionProduct.hh"


class G4FinalStateSampler
{
public:  // with description
  G4FinalStateSampler() { }
     
  virtual ~G4FinalStateSampler() { }
    
  enum { energyBins=30 };

protected:  // with description
  typedef std::pair<G4int, G4double> Interpolation;
  Interpolation interpolateEnergy(G4double ke) const;
  
  G4double interpolateCrossSection(const Interpolation& epair,
				   const G4double xsec[energyBins]) const;

  // Optional start/stop arguments default to inclusive arrays
  void fillSigmaBuffer(const Interpolation& epair, 
		       const G4double xsec_dMult[][energyBins],
		       G4int startBin=0, G4int stopBin=8) const;

  G4int sampleFlat() const;

  // This version will disappear
  G4int sampleFlat(const std::vector<G4double>& sigma) const;

  void CheckQnums(const G4FastVector<G4ReactionProduct,256> &vec,
		  G4int &vecLen,
		  G4ReactionProduct &currentParticle,
		  G4ReactionProduct &targetParticle,
		  G4double Q, G4double B, G4double S);
  
private:
  mutable std::vector<G4double> sigmaBuf;
  static const G4double energyScale[energyBins];
};

#endif
