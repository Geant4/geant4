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
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// class description
//
// An isotope is a chemical isotope defined by its name,
//                                                 Z: atomic number
//                                                 N: number of nucleons 
//                                                 A: mass of a mole (optional)
//                                                 m: isomer state (optional)
// If A is not defined it is taken from Geant4 database
//
// The class contains as a private static member the table of defined
// isotopes (an ordered vector of isotopes).
//
// Isotopes can be assembled into elements via the G4Element class.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// 20.08.11: Add flag fm for isomer level (mma)
// 15.11.05: GetIsotope(isotopeName, G4bool warning=false)
// 31.03.05: A becomes optional. Taken from Nist data base by default (mma)  
// 26.02.02: fIndexInTable renewed 
// 14.09.01: fCountUse: nb of elements which use this isotope 
// 13.09.01: stl migration. Suppression of the data member fIndexInTable
// 30.03.01: suppression of the warning message in GetIsotope
// 04.08.98: new method GetIsotope(isotopeName) (mma)
// 17.01.97: aesthetic rearrangement (mma)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#ifndef G4ISOTOPE_HH
#define G4ISOTOPE_HH

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

class G4Isotope;
using G4IsotopeTable = std::vector<G4Isotope*>;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class G4Isotope
{
 public:  // with description

    // Make an isotope
    //
    G4Isotope(const G4String& name,		//its name
                    G4int     z,		//atomic number
                    G4int     n,		//number of nucleons
                    G4double  a = 0.,		//mass of mole
                    G4int     mlevel = 0);      //isomer level		    
                    
    ~G4Isotope();

    // Retrieval methods
    //
    const G4String& GetName()  const {return fName;}

    // Atomic number 
    G4int    GetZ()     const {return fZ;}

    // Number of nucleous
    G4int    GetN()     const {return fN;}

    // Atomic mass of mole in Geant4 units with electron shell 
    G4double GetA()     const {return fA;}

    // Isomer level      
    G4int    Getm()     const {return fm;}
        
    static  
    G4Isotope* GetIsotope(const G4String& name, G4bool warning=false);
    
    static const
    G4IsotopeTable* GetIsotopeTable();
    
    static 
    size_t GetNumberOfIsotopes();
    
    size_t GetIndex() const {return fIndexInTable;}
    
    friend
    std::ostream& operator<<(std::ostream&, const G4Isotope*);
    
    friend
    std::ostream& operator<<(std::ostream&, const G4Isotope&);
    
    friend
    std::ostream& operator<<(std::ostream&, G4IsotopeTable);
     
 public:  // without description
 
    G4bool operator==(const G4Isotope&) const;
    G4bool operator!=(const G4Isotope&) const;

    G4Isotope(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.
   
    void SetName(const G4String& name) {fName=name;}

 private:
     
    G4Isotope(G4Isotope&);
    G4Isotope& operator=(const G4Isotope&);

 private:

    G4String fName;              // name of the Isotope
    G4int    fZ;                 // atomic number
    G4int    fN;                 // number of nucleons
    G4double fA;                 // atomic mass of a mole
    G4int    fm;                 // isomer level
        
    static 
    G4IsotopeTable theIsotopeTable;
    
    size_t fIndexInTable;        // index in the Isotope table
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
