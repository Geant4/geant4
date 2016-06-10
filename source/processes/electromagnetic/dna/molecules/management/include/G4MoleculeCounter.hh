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
#include <memory>

struct compDoubleWithPrecision
{
    bool operator() (const double& a, const double& b) const
    {
        if(std::fabs(a - b) < fPrecision)
        {
            return false;
        }
        else
        {
            return a < b;
        }
    }

    static G4ThreadLocal double fPrecision ;
};

typedef std::map<G4double, G4int, compDoubleWithPrecision> NbMoleculeAgainstTime;

class G4MoleculeCounter
{
public:
	typedef std::map<G4Molecule, NbMoleculeAgainstTime> CounterMapType;

#if __cplusplus > 199711L && !defined __clang__
    typedef std::unique_ptr<std::vector<G4Molecule> > RecordedMolecules;
#else
    typedef std::auto_ptr<std::vector<G4Molecule> > RecordedMolecules;
#endif

protected:
    G4MoleculeCounter();
    virtual ~G4MoleculeCounter(){;}
    static G4ThreadLocal G4MoleculeCounter* fpInstance;

    CounterMapType fCounterMap;
    std::map<const G4MoleculeDefinition*, G4bool> fDontRegister ;
    static G4bool fUse;

    G4int fVerbose ;

    friend class G4Molecule;
    virtual void AddAMoleculeAtTime(const G4Molecule&, G4double);
    virtual void RemoveAMoleculeAtTime(const G4Molecule&, G4double);

public:
    static void DeleteInstance();

    static  G4MoleculeCounter* GetMoleculeCounter();
    inline const NbMoleculeAgainstTime& GetNbMoleculeAgainstTime(const G4Molecule &molecule);

    RecordedMolecules GetRecordedMolecules();

    /*
     * The dynamics of the given molecule won't be saved into memory.
     */
    inline virtual void DontRegister(const G4MoleculeDefinition*);
    inline virtual void RegisterAll();

    /*
     * If the molecule counter is used, it will be called
     * at every creation/deletion of a molecule to
     * to increase/decrease the number at a given time.
     */
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

    /*
     * It sets the min time difference in between two time slices.
     */
    void SetTimeSlice(double);

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

inline void G4MoleculeCounter::RegisterAll()
{
    fDontRegister.clear();
}

#endif
