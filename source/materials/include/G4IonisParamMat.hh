// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonisParamMat.hh,v 1.5 2001-03-12 17:48:47 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class description
//
// The class contains few (physical) quantities related to the Ionisation
// process, for a material defined by its pointer G4Material*
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 09.03.01: copy constructor and assignement operator in public (mma)
// 09-07-98: data moved from G4Material (mma)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4IonisParamMat_HH
#define G4IonisParamMat_HH

#include "G4ios.hh"
#include "globals.hh"

class G4Material;                        // forward declaration

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4IonisParamMat  // with description
{
public:

    G4IonisParamMat(G4Material*); 
    virtual ~G4IonisParamMat();

    //
    // retrieval methods
    //
     
    // parameters for mean energy loss calculation:
    G4double  GetMeanExcitationEnergy()   const {return fMeanExcitationEnergy;};
    G4double  GetLogMeanExcEnergy()       const {return fLogMeanExcEnergy;};
    G4double* GetShellCorrectionVector()  const {return fShellCorrectionVector;};
    G4double  GetTaul()                   const {return fTaul;};
    
    // parameters of the density correction:
    G4double  GetCdensity()               const {return fCdensity;};
    G4double  GetMdensity()               const {return fMdensity;};
    G4double  GetAdensity()               const {return fAdensity;};
    G4double  GetX0density()              const {return fX0density;};
    G4double  GetX1density()              const {return fX1density;};

    // parameters of the energy loss fluctuation model:
    G4double  GetF1fluct()                const {return fF1fluct;};
    G4double  GetF2fluct()                const {return fF2fluct;};
    G4double  GetEnergy1fluct()           const {return fEnergy1fluct;};
    G4double  GetLogEnergy1fluct()        const {return fLogEnergy1fluct;};
    G4double  GetEnergy2fluct()           const {return fEnergy2fluct;};
    G4double  GetLogEnergy2fluct()        const {return fLogEnergy2fluct;};
    G4double  GetEnergy0fluct()           const {return fEnergy0fluct;};
    G4double  GetRateionexcfluct()        const {return fRateionexcfluct;};

public:  // without description

    G4IonisParamMat(const G4IonisParamMat&);
    const G4IonisParamMat& operator=(const G4IonisParamMat&);          
    G4int operator==(const G4IonisParamMat&) const;
    G4int operator!=(const G4IonisParamMat&) const;

private:
    
    // Compute mean parameters : ExcitationEnergy,Shell corretion vector ...
    void ComputeMeanParameters();

    // Compute parameters for the density effect
    void ComputeDensityEffect();

    // Compute parameters for the energy fluctuation model
    void ComputeFluctModel();

private:

//
// data members
//
    G4Material* fMaterial;                    // this material

   // parameters for mean energy loss calculation
    G4double  fMeanExcitationEnergy;          // 
    G4double  fLogMeanExcEnergy;              // 
    G4double* fShellCorrectionVector;         // shell correction coefficients
    G4double  fTaul;                          // lower limit of Bethe-Bloch formula

   // parameters of the density correction....
    G4double fCdensity;                      // mat.constant
    G4double fMdensity;                      // exponent
    G4double fAdensity;                      //
    G4double fX0density;                     //
    G4double fX1density;                     //

   // parameters of the energy loss fluctuation model
    G4double fF1fluct;                       
    G4double fF2fluct;                       
    G4double fEnergy1fluct;
    G4double fLogEnergy1fluct;
    G4double fEnergy2fluct;
    G4double fLogEnergy2fluct;
    G4double fEnergy0fluct;
    G4double fRateionexcfluct;
};

#endif
