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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4DNADAMAGE_HH
#define G4DNADAMAGE_HH 1

#include "G4Molecule.hh"

class G4VDNAHit
{
public :
    G4VDNAHit(){;}
    virtual ~G4VDNAHit(){;}
};

class G4DNAIndirectHit : public G4VDNAHit
{
public :
    G4DNAIndirectHit(const G4String& baseName, const G4Molecule* molecule,
                     const G4ThreeVector& position, G4double time);
    virtual ~G4DNAIndirectHit();

    inline const G4Molecule* GetMolecule() {return fpMolecule;}
    inline const G4ThreeVector& GetPosition() {return fPosition;}
    inline const G4String& GetBaseName() {return fBaseName;}
    inline double GetTime() {return fTime;}

    void Print();

protected :
    const G4Molecule* fpMolecule;
    G4ThreeVector fPosition;
    G4double fTime;
    G4String fBaseName;
};

class G4DNADamage
{
public:
    static G4DNADamage* Instance();
    static void DeleteInstance();

    virtual void Reset();

    //void AddDirectDamage();
    virtual void AddIndirectDamage(const G4String& baseName,
                                   const G4Molecule* molecule,
                                   const G4ThreeVector& position, 
                                   G4double time);

    inline const std::vector<G4DNAIndirectHit*>* GetIndirectHits();
    inline virtual G4int GetNIndirectHits() const
    {
        if(fJustCountDamage)
            return fNIndirectDamage;

        return (G4int)fIndirectHits.size();
    }

    inline virtual void SetOnlyCountDamage(G4bool flag = true)
    {
        fJustCountDamage = flag;
    }

    inline virtual G4bool OnlyCountDamage() const
    {
        return fJustCountDamage;
    }

protected :
    G4DNADamage();
    static G4ThreadLocal G4DNADamage* fpInstance;
    virtual ~G4DNADamage();

    G4bool fJustCountDamage;
    G4int fNIndirectDamage;
    std::vector<G4DNAIndirectHit*> fIndirectHits;
    std::map<G4Molecule, const G4Molecule*> fMolMap;
};

inline const std::vector<G4DNAIndirectHit*>* G4DNADamage::GetIndirectHits()
{
    return &fIndirectHits;
}

#endif // G4DNADAMAGES_HH
