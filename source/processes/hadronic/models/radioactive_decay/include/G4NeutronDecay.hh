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
//  File:   G4NeutronDecay.hh                                                 //
//  Author: L.G. Sarmiento (Lund)                                             //
//  Date:   10 October 2015                                                   //
//  Description: performs protom emission from radioactive nuclei, and        //
//               returns daughter particles in rest frame of parent nucleus   //
//                                                                            //
//               This class is created based on G4AlphaDecay                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4NeutronDecay_h
#define G4NeutronDecay_h 1

#include "G4NuclearDecay.hh"


class G4NeutronDecay : public G4NuclearDecay
{
  public:
    G4NeutronDecay(const G4ParticleDefinition* theParentNucleus,
                  const G4double& theBR, const G4double& Qvalue,
                  const G4double& excitation, const G4Ions::G4FloatLevelBase& flb);

    virtual ~G4NeutronDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    virtual void DumpNuclearInfo();

  private:
    const G4double transitionQ;
};
#endif

