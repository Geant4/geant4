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
//
// GEANT4 physics class: G4PhotoNuclearCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02

#ifndef G4PhotoNuclearCrossSection_h
#define G4PhotoNuclearCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include <vector>

class G4PhotoNuclearCrossSection : public G4VCrossSectionDataSet
{
public:
    
    G4PhotoNuclearCrossSection();
    ~G4PhotoNuclearCrossSection() override;
    
    static const char* Default_Name() {return "PhotoNuclearXS";}
    
    void CrossSectionDescription(std::ostream&) const override;

    G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Element* elm = nullptr,
			   const G4Material* mat = nullptr) override;
    
    G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
                               const G4Material* mat = nullptr) override;

    G4double GetIsoCrossSection(const G4DynamicParticle*,
                                G4int Z, G4int A,
                                const G4Isotope* iso = nullptr,
                                const G4Element* elm = nullptr,
                                const G4Material* mat = nullptr) override;
    
    G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
                                    const G4Material*) override;

    G4double ComputeElementXSection(G4double energy, G4int Z);

    G4double ComputeIsoXSection(G4double energy, G4int Z, G4int A);

    G4PhotoNuclearCrossSection& operator=
    (const G4PhotoNuclearCrossSection& right) = delete;
    G4PhotoNuclearCrossSection(const G4PhotoNuclearCrossSection&) = delete;

private:
    
    G4int GetFunctions(G4double a, G4double* y, G4double* z);
    G4double EquLinearFit(G4double X, G4int N, const G4double X0,
                          const G4double XD, const G4double* Y);
    G4double ThresholdEnergy(G4int Z, G4int N);
    
    G4int     lastZ = 0;         // The last Z of calculated nucleus
    G4double  lastSig = 0.0;     // Last value of the Cross Section
    G4double* lastGDR = nullptr; // Pointer to the last array of GDR cross sections
    G4double* lastHEN = nullptr; // Pointer to the last array of HEn cross sections
    G4double  lastE = 0.0;       // Last used in the cross section Energy
    G4double  lastTH = 0.0;      // Last value of the Energy Threshold (A-dependent)
    G4double  lastSP = 0.0;      // Last value of the ShadowingPomeron (A-dependent)
    
    // Vector of pointers to the GDRPhotonuclearCrossSection
    std::vector <G4double*> GDR;

    // store deuteron, triton, He3 XS
    G4double* deuteron_GDR = nullptr;
    G4double* deuteron_HR = nullptr;
    G4double deuteron_TH = 0.0;
    G4double deuteron_SP = 0.0;
    G4double* triton_GDR = nullptr;
    G4double* triton_HR = nullptr;
    G4double triton_TH = 0.0;
    G4double triton_SP = 0.0;
    G4double* he3_GDR = nullptr;
    G4double* he3_HR = nullptr;
    G4double he3_TH = 0.0;
    G4double he3_SP = 0.0;
    
    // Vector of pointers to the HighEnPhotonuclearCrossSect
    std::vector <G4double*> HEN;
    
    std::vector <G4double> spA;  // shadowing coefficients (A-dependent)
    std::vector <G4double> eTH;    // energy threshold (A-dependent)
    
    G4NistManager* nistmngr;
    
    G4double mNeut;
    G4double mProt;
};

#endif
