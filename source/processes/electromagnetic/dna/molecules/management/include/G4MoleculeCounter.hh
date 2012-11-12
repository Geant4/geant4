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
#ifndef G4MoleculeCounter_h
#define G4MoleculeCounter_h

#include "G4Molecule.hh"
#include <map>
#include "G4ParticleDefinition.hh"
#include <memory>

using namespace std;

struct compDoubleWithPrecision
{
    bool operator() (const double& a, const double& b) const
    {
        if(fabs(a - b) < fPrecision)
        {
            return false;
        }
        else
        {
            return a < b;
        }
    }

    static double fPrecision ;
};

typedef map<G4double, G4int, compDoubleWithPrecision> NbMoleculeAgainstTime;

class G4MoleculeCounter
{
protected:
    G4MoleculeCounter();
    virtual ~G4MoleculeCounter(){;}
    static G4MoleculeCounter* fpInstance;
    typedef std::map<G4Molecule, NbMoleculeAgainstTime> CounterMapType;

    CounterMapType fCounterMap;

    std::map<const G4MoleculeDefinition*, G4bool> fDontRegister ;
    G4bool fUse;

    G4int fVerbose ;

public:
    static void DeleteInstance();

    static  G4MoleculeCounter* GetMoleculeCounter();
    inline const NbMoleculeAgainstTime& GetNbMoleculeAgainstTime(const G4Molecule &molecule);
    std::auto_ptr<vector<G4Molecule> > GetRecordedMolecules();
    virtual void AddAMoleculeAtTime(const G4Molecule&, G4double);
    virtual void RemoveAMoleculeAtTime(const G4Molecule&, G4double);

    // inline void DontRegister(G4MoleculeID);
    inline virtual void DontRegister(const G4MoleculeDefinition*);
    inline virtual void ResetDontRegister();

    void Use(G4bool flag = true)
    {
        fUse=flag;
    }
    G4bool InUse()
    {
        return fUse;
    }

    inline void SetVerbose(G4int);
    inline G4int GetVerbose();

    virtual void ResetCounter();
};

inline void G4MoleculeCounter::ResetCounter()
{
    fCounterMap.clear();
}

inline const NbMoleculeAgainstTime& G4MoleculeCounter::GetNbMoleculeAgainstTime(const G4Molecule& molecule)
{
    return fCounterMap[molecule];
}

inline void G4MoleculeCounter::SetVerbose(G4int level)
{
    fVerbose = level;
}

inline G4int G4MoleculeCounter::GetVerbose()
{
    return fVerbose ;
}

inline void G4MoleculeCounter::DontRegister(const G4MoleculeDefinition* molDef)
{
    fDontRegister[molDef] = true ;
}

inline void G4MoleculeCounter::ResetDontRegister()
{
    fDontRegister.clear();
}

#endif
