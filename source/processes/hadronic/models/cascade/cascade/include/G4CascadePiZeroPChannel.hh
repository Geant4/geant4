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
// $Id: G4CascadePiZeroPChannel.hh 67796 2013-03-08 06:18:39Z mkelsey $

#ifndef G4_CASCADE_PIZEROP_CHANNEL_HH
#define G4_CASCADE_PIZEROP_CHANNEL_HH

#include "G4CascadeData.hh"
#include "G4CascadeFunctions.hh"
#include "G4PionNucSampler.hh"

struct G4CascadePiZeroPChannelData {
  typedef G4CascadeData<30,5,13,21,30,6,7,8,9> data_t;
  static const data_t data;
};

typedef G4CascadeFunctions<G4CascadePiZeroPChannelData,G4PionNucSampler> G4CascadePiZeroPChannel;

#endif	/* G4_CASCADE_PIZEROP_CHANNEL_HH */
