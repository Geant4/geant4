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
//  File:   G4LENDGammaCrossSection.hh                                        //
//  Date:   30 March 2020                                                     //
//  Author: Dennis H. Wright                                                  //
//                                                                            //
//  Description: cross sections for inelastic scattering of gammas from       //
//               nuclei including gamma-induced fission.  This cross section  //
//               is very similar to G4LENDCombinedCrossSection except that    //
//               does not sample elastic or capture reactions since there are //
//               no such data for gammas in GND.                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4LENDGammaCrossSection_h
#define G4LENDGammaCrossSection_h 1

#include "G4LENDCrossSection.hh"

class G4LENDInelasticCrossSection;
class G4LENDFissionCrossSection;
class G4HadProjectile;

class G4LENDGammaCrossSection : public G4LENDCrossSection
{
  public:
    G4LENDGammaCrossSection(G4ParticleDefinition* pd);
    ~G4LENDGammaCrossSection(){;};

    void BuildPhysicsTable(const G4ParticleDefinition&);

    G4double GetIsoCrossSection(const G4DynamicParticle*, G4int /*Z*/,
                                G4int /*A*/, const G4Isotope*,
                                const G4Element*, const G4Material*);

    G4int SelectChannel(const G4DynamicParticle*, G4int /*Z*/,
                        G4int /*A*/, const G4Isotope*, const G4Element*,
                        const G4Material*);

  private:
    G4LENDInelasticCrossSection* inelasticXS;
    G4LENDFissionCrossSection* fissionXS;

};
#endif
