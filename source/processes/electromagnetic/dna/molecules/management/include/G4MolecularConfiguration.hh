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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------


#ifndef G4MolecularConfiguration_
#define G4MolecularConfiguration_ 1
#include <G4MoleculeDefinition.hh>
#include <map>
#include <vector>
#include <CLHEP/Utility/memory.h>

struct comparator;

class G4MolecularConfiguration;
class G4MoleculeDefinition;
class G4ElectronOccupancy;
class G4MolecularDecayChannel;

/** The pointer G4MolecularConfiguration will be shared by all the
* molecules having the same molecule definition and the same
* electron occupancy
* BE CAREFUlL !!! : If you change the mass for instance of a OH^-,
* this will affect all the OH^- molecule diffusing around
*/
class G4MolecularConfiguration
{
public :

    /////////////////
    // Static methods

    // Get for a given moleculeDefinition and a given electronic configuration, the mol conf
    static G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition*,
            const G4ElectronOccupancy& electronOccupancy);

    // Get ground state electronic configuration
    static G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition*);

    // Release memory of the mol conf manager
    static void DeleteManager();
    ///////////////

    // Methods
    const G4MoleculeDefinition* GetDefinition() const;

    /** Returns the name of the molecule
    */
    const G4String& GetName() const;

    /** Returns the nomber of atoms compouning the molecule
       */
    G4int GetAtomsNumber() const;

    /** Method used in Geant4-DNA to excite water molecules
    */
    G4MolecularConfiguration* ExciteMolecule(G4int);

    /** Method used in Geant4-DNA to ionize water molecules
    */
    G4MolecularConfiguration* IonizeMolecule(G4int);

    /** Add n electrons to a given orbit.
    * Note : You can add as many electrons to a given orbit, the result
    * may be unrealist.
    */
    G4MolecularConfiguration* AddElectron(G4int orbit, G4int n =1);

    /** Remove n electrons to a given orbit.
    */
    G4MolecularConfiguration* RemoveElectron(G4int,G4int number=1);

    /** Move one electron from an orbit to another.
    */
    G4MolecularConfiguration* MoveOneElectron(G4int /*orbit*/,G4int /*orbit*/);

    /** Returns the number of electron.
    */
    G4double GetNbElectrons() const;

    /** Show the electronic state of the molecule.
    */
    void PrintState() const;

    const std::vector <const G4MolecularDecayChannel*>* GetDecayChannel() const;

    G4int GetMoleculeID() const;

    /** Sets the diffusion coefficient D of the molecule used in diffusion
       * processes to calculate the mean square jump distance between two
       * changes of direction. In three dimension : <x^2> = 6 D t where t is
       * the mean jump time between two changes of direction.
       */
    inline void SetDiffusionCoefficient(G4double);

    /** Returns the diffusion coefficient D.
    */
    inline G4double GetDiffusionCoefficient() const;

    /** Set the decay time of the molecule.
    */
    inline void SetDecayTime(G4double);

    /** Returns the decay time of the molecule.
    */
    inline G4double GetDecayTime() const;

    /** The Van Der Valls Radius of the molecule
    */
    inline void SetVanDerVaalsRadius(G4double);
    inline G4double GetVanDerVaalsRadius() const ;

    /** Returns the object ElectronOccupancy describing the electronic
    * configuration of the molecule.
    */
    inline const G4ElectronOccupancy* GetElectronOccupancy() const;

    /** Returns the charge of molecule.
    */
    inline G4int GetCharge() const;

    /** Set the total mass of the molecule.
    */
    inline void SetMass(G4double);

    /** Returns the total mass of the molecule.
    */
    inline G4double GetMass() const;

protected :
    G4MolecularConfiguration(const G4MoleculeDefinition*, const G4ElectronOccupancy&);
    G4MolecularConfiguration(const G4MolecularConfiguration&);
    G4MolecularConfiguration & operator=(G4MolecularConfiguration &right);
    ~G4MolecularConfiguration();
    G4MolecularConfiguration* ChangeConfiguration(const G4ElectronOccupancy& newElectronOccupancy);

    const G4MoleculeDefinition* fMoleculeDefinition;
    const G4ElectronOccupancy* fElectronOccupancy;

    struct G4MolecularConfigurationManager
    {
        G4MolecularConfigurationManager(){;}
        ~G4MolecularConfigurationManager();

        typedef std::map<const G4MoleculeDefinition*, std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator> > MolecularConfigurationTable;
        MolecularConfigurationTable fTable;
    };

    static G4MolecularConfigurationManager* fgManager;

    static G4MolecularConfigurationManager* GetManager();

    G4double fDynDiffusionCoefficient;
    G4double fDynVanDerVaalsRadius;
    G4double fDynDecayTime;
    G4double fDynMass;
    G4int    fDynCharge;
    mutable G4String fName; // mutable allowed this member to be changed in const methods
};

struct comparator
{
    bool operator() (const G4ElectronOccupancy& occ1, const G4ElectronOccupancy& occ2) const
    {
        // Since this method is called a lot of time,
        // we retrieve only once the totOcc
        G4int totalOcc1 = occ1.GetTotalOccupancy() ;
        G4int totalOcc2 = occ2.GetTotalOccupancy() ;
        if ( totalOcc1!= totalOcc2)
        {
            return totalOcc1<totalOcc2;
        }
        else
        {
            G4int occupancy1 = -1 ;
            G4int occupancy2 = -1 ;
            const G4int sizeOrbit = occ1.GetSizeOfOrbit() ;
            for (G4int i=0; i<occ1.GetSizeOfOrbit();)
            {
                // Since this method is called a lot of time,
                // we retrieve only once the Occ

                occupancy1 = occ1.GetOccupancy(i);
                occupancy2 = occ2.GetOccupancy(i);

                if (occupancy1 != occupancy2)
                {
                    return occupancy1 < occupancy2;
                }
                else
                {
                    i++;
                    if (i >= sizeOrbit) return false;
                }
            }
        }
        return false;
    }
};


inline const G4MoleculeDefinition* G4MolecularConfiguration::GetDefinition() const
{
    return fMoleculeDefinition;
}

inline const G4ElectronOccupancy* G4MolecularConfiguration::GetElectronOccupancy() const
{
    return fElectronOccupancy ;
}

inline void G4MolecularConfiguration::SetDiffusionCoefficient(G4double dynDiffusionCoefficient)
{
    fDynDiffusionCoefficient = dynDiffusionCoefficient ;
}

inline G4double G4MolecularConfiguration::GetDiffusionCoefficient() const
{
    return fDynDiffusionCoefficient;
}

inline void G4MolecularConfiguration::SetDecayTime(G4double dynDecayTime)
{
    fDynDecayTime = dynDecayTime;
}

inline G4double G4MolecularConfiguration::GetDecayTime() const
{
    return fDynDecayTime;
}

inline void G4MolecularConfiguration::SetVanDerVaalsRadius(G4double dynVanDerVaalsRadius)
{
    fDynVanDerVaalsRadius = dynVanDerVaalsRadius ;
}

inline G4double G4MolecularConfiguration::GetVanDerVaalsRadius() const
{
    return fDynVanDerVaalsRadius;
}

inline G4int G4MolecularConfiguration::GetCharge() const
{
    return fDynCharge ;
}

inline void G4MolecularConfiguration::SetMass(G4double aMass)
{
    fDynMass = aMass ;
}

inline G4double G4MolecularConfiguration::GetMass() const
{
    return fDynMass;
}
#endif
