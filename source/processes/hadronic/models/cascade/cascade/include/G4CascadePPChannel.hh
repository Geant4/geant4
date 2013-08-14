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
// $Id: G4CascadePPChannel.hh 67796 2013-03-08 06:18:39Z mkelsey $
//
// 20120907  M. Kelsey -- Subclass and overload findCrossSection() function.

#ifndef G4_CASCADE_PP_CHANNEL_HH
#define G4_CASCADE_PP_CHANNEL_HH

#include "G4CascadeData.hh"
#include "G4CascadeFunctions.hh"
#include "G4PionNucSampler.hh"

struct G4CascadePPChannelData {
  typedef G4CascadeData<30,1,6,18,32,7,8,10,11> data_t;
  static const data_t data;
};


class G4CascadePPChannel
  : public G4CascadeFunctions<G4CascadePPChannelData,G4PionNucSampler> {
public:
  G4CascadePPChannel()
    : G4CascadeFunctions<G4CascadePPChannelData,G4PionNucSampler>() {;}
  virtual ~G4CascadePPChannel() {;}

  // Will replace interpolation of 0-10 MeV bin on total and elastic
  virtual G4double 
  findCrossSection(G4double ke, const G4double (&xsec)[30]) const;
};

#endif	/* G4_CASCADE_PP_CHANNEL_HH */
