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
// $Id: G4IonisParamElm.hh 96794 2016-05-09 10:09:30Z gcosmo $
//

// class description
//
// The class contains few (physical) quantities related to the Ionisation
// process, for the Element defined by its atomic number Z
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 09.03.01: copy constructor and assignement operator in public (mma)
// 09.07.98: data moved from G4Element (mma) 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4IonisParamElm_HH
#define G4IonisParamElm_HH

#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4IonisParamElm
{
public:  // with description

    G4IonisParamElm(G4double Z);
    ~G4IonisParamElm();

    // retrieval methods
    
    G4double  GetZ()        const {return fZ;}       // atomic number Z
    G4double  GetZ3()       const {return fZ3;}      // std::pow (Z,1/3)
    G4double  GetZZ3()      const {return fZZ3;}     // std::pow (Z(Z+1),1/3)
    G4double  GetlogZ3()    const {return flogZ3;}   // std::log(Z)/3

    G4double  GetTau0() const {return fTau0;};
                       // 0.1*std::pow(Z,1/3)*MeV/proton_mass_c2
    
    G4double  GetTaul() const {return fTaul;}        // 2*MeV/proton mass
    
    G4double  GetAlow() const {return fAlow;}
    G4double  GetBlow() const {return fBlow;}
    G4double  GetClow() const {return fClow;}
                       // parameters for the low energy ion.loss
    
    G4double  GetMeanExcitationEnergy()  const {return fMeanExcitationEnergy;}
                       // ICRU'37 report 

    G4double  GetFermiVelocity()          const {return fVFermi;};
    G4double  GetLFactor()                const {return fLFactor;};
      
    G4double* GetShellCorrectionVector() const {return fShellCorrectionVector;}
                                       // shell correction coefficients
        
    G4IonisParamElm(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

private:

    // operators
    G4IonisParamElm& operator=(const G4IonisParamElm&) = delete;
    G4int operator==(const G4IonisParamElm&) const = delete;
    G4int operator!=(const G4IonisParamElm&) const = delete;
    G4IonisParamElm(G4IonisParamElm&) = delete;
    //
    //  data members
    //
    G4double fZ;                 // effective Z
    G4double fZ3;                // std::pow (Z,1/3)
    G4double fZZ3;               // std::pow (Z(Z+1),1/3)
    G4double flogZ3;             // std::log(Z)/3

    //  ------ ionisation loss ---------------------------------
    
    G4double  fTau0 ;                 // 0.1*std::pow(Z,1/3)*MeV/proton_mass_c2
    G4double  fTaul ;                 // 2*MeV/proton mass
    G4double  fBetheBlochLow;         // Bethe-Bloch at fTaul*particle mass   
    G4double  fAlow,fBlow,fClow;      // parameters for the low energy ion.loss
    G4double  fMeanExcitationEnergy;  //     
    G4double* fShellCorrectionVector; // shell correction coefficients

    // parameters for ion corrections computations
    G4double fVFermi;
    G4double fLFactor;
};

#endif
