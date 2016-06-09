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
// $Id: G4MolecularConfiguration.cc 65022 2012-11-12 16:43:12Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4MolecularConfiguration.hh"
#include "G4UIcommand.hh"

using namespace std;

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// G4MolecularConfigurationManager
G4MolecularConfiguration::G4MolecularConfigurationManager* G4MolecularConfiguration::fgManager = 0 ;

G4MolecularConfiguration::G4MolecularConfigurationManager*
        G4MolecularConfiguration::GetManager()
{
    if(!fgManager)
    {
        fgManager = new G4MolecularConfiguration::G4MolecularConfigurationManager;
    }

    return fgManager;
}

G4MolecularConfiguration::G4MolecularConfigurationManager::~G4MolecularConfigurationManager()
{
    G4MolecularConfigurationManager::MolecularConfigurationTable::iterator it1;
    std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it2;

    for(it1 = fTable.begin() ; it1 != fTable.end() ; it1++)
    {
        for(it2=it1->second.begin(); it2!=it1->second.end(); it2++)
        {
            if(it2->second)
            {
                delete it2->second;
            }
        }
    }
}

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Static method in G4MolecularConfiguration
G4MolecularConfiguration* G4MolecularConfiguration::GetMolecularConfiguration(const G4MoleculeDefinition* molDef)
{
    const G4ElectronOccupancy& elecOcc = *molDef->GetGroundStateElectronOccupancy();
    if(GetManager()->fTable[molDef][elecOcc])
    {
        return GetManager()->fTable[molDef][elecOcc];
    }
    else
    {
        G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef, elecOcc);
        return newConf ;
    }
}

G4MolecularConfiguration* G4MolecularConfiguration::GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
        const G4ElectronOccupancy& elecOcc )
{
    if(GetManager()->fTable[molDef][elecOcc])
    {
        return GetManager()->fTable[molDef][elecOcc];
    }
    else
    {
        G4MolecularConfiguration* newConf = new G4MolecularConfiguration(molDef, elecOcc);
        return newConf ;
    }
}

void G4MolecularConfiguration::DeleteManager()
{
    if(fgManager) delete fgManager;
    fgManager = 0;
}

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// G4MolecularConfiguration
G4MolecularConfiguration::G4MolecularConfiguration(const G4MoleculeDefinition* moleculeDef,
        const G4ElectronOccupancy& elecOcc)
{
    fMoleculeDefinition = moleculeDef ;
    fgManager->fTable[fMoleculeDefinition][elecOcc] = this;
    std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator>::iterator it ;
    it = fgManager->fTable[moleculeDef].find(elecOcc);
    fElectronOccupancy = &(it->first);

    fDynCharge = fMoleculeDefinition->GetNbElectrons()-fElectronOccupancy->GetTotalOccupancy();
    fDynMass = fMoleculeDefinition->GetMass() ;

    fDynDiffusionCoefficient = fMoleculeDefinition->GetDiffusionCoefficient() ;
    fDynVanDerVaalsRadius = fMoleculeDefinition->GetVanDerVaalsRadius() ;
    fDynDecayTime = fMoleculeDefinition->GetDecayTime() ;
}

G4MolecularConfiguration::~G4MolecularConfiguration()
{
    if(fElectronOccupancy)
    {
        delete fElectronOccupancy;
        fElectronOccupancy = 0;
    }
}

G4MolecularConfiguration* G4MolecularConfiguration::ChangeConfiguration(const G4ElectronOccupancy& newElectronOccupancy)
{
    G4MolecularConfiguration* output = fgManager->fTable[fMoleculeDefinition][newElectronOccupancy] ;
    if(! output)
    {
        output = new G4MolecularConfiguration(fMoleculeDefinition, newElectronOccupancy);
    }
    return output ;
}

G4MolecularConfiguration& G4MolecularConfiguration::operator=(G4MolecularConfiguration& right)
{
    if (&right==this) return *this;
    return *this;
}


/** Method used in Geant4-DNA to excite water molecules
 */
G4MolecularConfiguration* G4MolecularConfiguration::ExciteMolecule(G4int ExcitedLevel)
{
    G4ElectronOccupancy newElectronOccupancy (*fElectronOccupancy);

    newElectronOccupancy.RemoveElectron(ExcitedLevel,1);
    newElectronOccupancy.AddElectron(5,1);

    return ChangeConfiguration(newElectronOccupancy);
}

/** Method used in Geant4-DNA to ionize water molecules
 */
G4MolecularConfiguration* G4MolecularConfiguration::IonizeMolecule(G4int IonizedLevel)
{
    G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);

    if(newElectronOccupancy.GetOccupancy(IonizedLevel) != 0)
    {
        newElectronOccupancy.RemoveElectron(IonizedLevel,1);
    }
    else
    {
        G4String errMsg = "There is no electron on the orbit " + G4UIcommand::ConvertToString(IonizedLevel) +
                          " you want to free. The molecule's name you want to ionized is "+ GetName();
        G4Exception("G4Molecule::IonizeMolecule","",FatalErrorInArgument, errMsg);
        PrintState();
    }

    // PrintState();

    return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::AddElectron(G4int orbit, G4int number)
{
    G4ElectronOccupancy newElectronOccupancy(*fElectronOccupancy);
    newElectronOccupancy.AddElectron(orbit, number);
    return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::RemoveElectron(G4int orbit,G4int number)
{
    G4ElectronOccupancy newElectronOccupancy (*fElectronOccupancy);

    if(newElectronOccupancy.GetOccupancy(orbit) != 0)
    {
        newElectronOccupancy.RemoveElectron(orbit, number );
    }
    else
    {
        G4String errMsg = "There is already no electron into the orbit " + G4UIcommand::ConvertToString(orbit) +
                          " you want to free. The molecule's name is "+ GetName();
        G4Exception("G4Molecule::RemoveElectron","",JustWarning, errMsg);
        PrintState();
    }

    return ChangeConfiguration(newElectronOccupancy);
}

G4MolecularConfiguration* G4MolecularConfiguration::MoveOneElectron(G4int orbitToFree,G4int orbitToFill)
{
    G4ElectronOccupancy newElectronOccupancy (*fElectronOccupancy);

    if(newElectronOccupancy . GetOccupancy(orbitToFree)>=1)
    {
        newElectronOccupancy . RemoveElectron(orbitToFree,1);
        newElectronOccupancy . AddElectron(orbitToFill,1);
    }
    else
    {
        G4String errMsg = "There is no electron on the orbit " + G4UIcommand::ConvertToString(orbitToFree) +
                          " you want to free. The molecule's name is "+ GetName();
        G4Exception("G4Molecule::MoveOneElectron","",FatalErrorInArgument, errMsg);
        PrintState();
    }

    return ChangeConfiguration(newElectronOccupancy);
}

const G4String& G4MolecularConfiguration::GetName() const
{
    if(fName.isNull())
    {
        fName = fMoleculeDefinition->GetName();
        fName+= "^";
        fName+= "{";
        fName+= G4UIcommand::ConvertToString(fDynCharge);
        fName+= "}";
    }
    return fName;
}

G4int G4MolecularConfiguration::GetAtomsNumber() const
{
    return fMoleculeDefinition->GetAtomsNumber();
}

G4double G4MolecularConfiguration::GetNbElectrons() const
{
    return fElectronOccupancy->GetTotalOccupancy();
}

void G4MolecularConfiguration::PrintState() const
{
    G4cout<<"--------------Print electronic state of "<<GetName()<<"---------------"<<G4endl;
    fElectronOccupancy->DumpInfo();
    if(fElectronOccupancy==fMoleculeDefinition->GetGroundStateElectronOccupancy())
    {
        G4cout<<"At ground state"<<G4endl;
    }
    else
    {
        if(fMoleculeDefinition->GetDecayTable())
            G4cout<<"Transition :"<<(fMoleculeDefinition->GetDecayTable())->GetExcitedState(fElectronOccupancy)<<G4endl;
    }
}

// added - to be transformed in a "Decay method"
const vector <const G4MolecularDecayChannel*>* G4MolecularConfiguration::GetDecayChannel() const
{
    return fMoleculeDefinition-> GetDecayChannels(fElectronOccupancy);
}

G4int G4MolecularConfiguration::GetMoleculeID() const
{
    if(fMoleculeDefinition)
        return fMoleculeDefinition->GetPDGEncoding();
    else
        G4Exception("G4Molecule::GetMoleculeID","",FatalErrorInArgument, "You should first enter a molecule defintion");

    return INT_MAX;
}
