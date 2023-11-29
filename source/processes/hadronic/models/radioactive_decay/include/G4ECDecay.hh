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
//  File:   G4ECDecay.hh                                                      //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   25 November 2014                                                  //
//  Description: performs electron capture decay of radioactive nuclei, and   //
//               returns daughter particles in rest frame of parent nucleus   // 
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4ECDecay_h
#define G4ECDecay_h 1

#include "G4NuclearDecay.hh"


class G4ECDecay : public G4NuclearDecay
{
  public:
    G4ECDecay(const G4ParticleDefinition* theParentNucleus,
              const G4double& theBR, const G4double& Qvalue,
              const G4double& excitation, const G4Ions::G4FloatLevelBase& flb,
              const G4RadioactiveDecayMode& mode);

    virtual ~G4ECDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    virtual void DumpNuclearInfo();

    void SetARM(G4bool onoff) {applyARM = onoff;}

  private:
    void DefineSubshellProbabilities(G4int Z, G4int A);


  private:
    const G4double transitionQ;
    G4double PL1,PL2,PM1,PM2,PN1,PN2;
    G4bool applyARM;

    //Ratio of subshells probability
    static const G4double PL2overPL1[100];
    static const G4double PM2overPM1[100];
    static const G4double PN2overPN1[100];

};
#endif

