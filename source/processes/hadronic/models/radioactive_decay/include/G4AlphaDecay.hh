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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4AlphaDecay.hh                                                   //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   20 October 2014                                                   //
//  Description: performs alpha emission from radioactive nuclei, and         //
//               returns daughter particles in rest frame of parent nucleus   // 
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4AlphaDecay_h
#define G4AlphaDecay_h 1

#include "G4VDecayChannel.hh"


class G4AlphaDecay : public G4VDecayChannel
{
  public:
    G4AlphaDecay(const G4ParticleDefinition* theParentNucleus,
                 const G4double& theBR, const G4double& Qvalue,
                 const G4double& excitation);

    G4AlphaDecay(const G4AlphaDecay&);

    ~G4AlphaDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    void SetHLThreshold(G4double HLT) {halflifeThreshold = HLT;}
    G4double GetHLThreshold() {return halflifeThreshold;}

    void DumpInfo();

  private:
    const G4double transitionQ;
    const G4double daughterEx;
    G4double halflifeThreshold;  // for variance reduction mode
};
#endif

