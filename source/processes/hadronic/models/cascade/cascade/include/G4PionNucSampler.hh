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
// $Id: G4PionNucSampler.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100512  M. Kelsey -- Replaces G4FinalStateSampler

#ifndef G4_PION_NUC_SAMPLER_HH
#define G4_PION_NUC_SAMPLER_HH

#include "G4CascadeSampler.hh"

// Entire interface is inherited from base (struct makes default public:)
struct G4PionNucSampler : public G4CascadeSampler<30,8> { G4PionNucSampler(); };

#endif	/* G4_PION NUC_SAMPLER_HH */
