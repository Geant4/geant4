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
// $Id: G4CascadeKplusPChannel.hh 67796 2013-03-08 06:18:39Z mkelsey $
//
// 20100507  M. Kelsey -- Remove redundant total-bins template argument
// 20100510  M. Kelsey -- Add initial "31" template arg.  Add G4CascSampler
//		to template for channel typedef
// 20100514  M. Kelsey -- Replace G4CascadeSampler with G4KaonHypSampler.
// 20150713  M. Kelsey -- Replace G4KaonHypSampler with G4KaonSampler, new
//		template args up to 9-body final states.

#ifndef G4_CASCADE_KPLUSP_CHANNEL_HH
#define G4_CASCADE_KPLUSP_CHANNEL_HH

#include "G4CascadeData.hh"
#include "G4CascadeFunctions.hh"
#include "G4KaonSampler.hh"

struct G4CascadeKplusPChannelData {
  typedef G4CascadeData<30,1,6,16,29,42,54,41,47> data_t;
  static const data_t data;
};

typedef G4CascadeFunctions<G4CascadeKplusPChannelData,G4KaonSampler> G4CascadeKplusPChannel;

#endif
