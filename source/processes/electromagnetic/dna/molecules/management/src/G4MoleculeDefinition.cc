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
// $Id: G4MoleculeDefinition.cc 93883 2015-11-03 08:25:04Z gcosmo $
//
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      21 Oct 2009 first implementation by A. Mantero and M.Karamitros
//                  Based on prototype of A.Mantero
// **********************************************************************

#include "G4MoleculeDefinition.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4Serialize.hh"

using namespace std;

// -----------------------------------------------------------------------------
// ###                          MoleculeDefinition                           ###
// -----------------------------------------------------------------------------

G4MoleculeDefinition::G4MoleculeDefinition(const G4String& name,
                                           G4double mass,
                                           G4double diffCoeff,
                                           G4int charge,
                                           G4int electronicLevels,
                                           G4double radius,
                                           G4int atomsNumber,
                                           G4double lifetime,
                                           G4String aType,
                                           G4FakeParticleID ID) :
    G4ParticleDefinition(name, mass, 0., charge, 0, 0, 0, 0, 0, 0, "Molecule",
                         0, 0, ID, false, lifetime, NULL, false, aType, 0, 0.0),
    fDiffusionCoefficient(diffCoeff),
    fAtomsNb(atomsNumber),
    fVanDerVaalsRadius(radius)

{
  fCharge = charge;
  if(electronicLevels)
  {
    fElectronOccupancy = new G4ElectronOccupancy(electronicLevels);
  }
  else
  {
    fElectronOccupancy = 0;
  }
  fDecayTable = NULL;
  G4MoleculeTable::Instance()->Insert(this);
}

G4MoleculeDefinition* G4MoleculeDefinition::Load(std::istream& in)
{
  G4String name;
  G4double mass;
  G4double diffCoeff;
  G4int charge;
  G4int electronicLevels;
  G4double radius;
  G4int atomsNumber;
  G4double lifetime;
  G4String aType;

  READ(in,name);
  READ(in,mass);
  READ(in,diffCoeff);
  READ(in,charge);
  READ(in,electronicLevels);
  READ(in,radius);
  READ(in,atomsNumber);
  READ(in,lifetime);
  READ(in,aType);

  return new G4MoleculeDefinition(name,
                                  mass,
                                  diffCoeff,
                                  charge,
                                  electronicLevels,
                                  radius,
                                  atomsNumber,
                                  lifetime,
                                  aType);
}

void G4MoleculeDefinition::Serialize(std::ostream& out)
{
  WRITE(out,this->G4ParticleDefinition::GetParticleName());
  WRITE(out,this->G4ParticleDefinition::GetPDGMass());
  WRITE(out,this->G4ParticleDefinition::GetPDGLifeTime());
  WRITE(out,this->G4ParticleDefinition::GetParticleType());
  WRITE(out,fDiffusionCoefficient);
  WRITE(out,fCharge);
  if(fElectronOccupancy)
  {
    WRITE(out,fElectronOccupancy->GetSizeOfOrbit());
  }
  else
  {
    WRITE(out,(int) 0);
  }
  WRITE(out,fVanDerVaalsRadius);
  WRITE(out,fAtomsNb);
}

//______________________________________________________________________________

G4MoleculeDefinition::~G4MoleculeDefinition()
{
  if (fElectronOccupancy)
  {
    delete fElectronOccupancy;
    fElectronOccupancy = 0;
  }
  if (fDecayTable)
  {
    delete fDecayTable;
    fDecayTable = 0;
  }
}

//___________________________________________________________________________

G4MolecularConfiguration*
 G4MoleculeDefinition::NewConfiguration(const G4String& molConfLabel)
{
  bool alreadyExist(false);
  return G4MolecularConfiguration::CreateMolecularConfiguration(GetName()+"_"+
                                                                molConfLabel,
                                                                this,
                                                                molConfLabel,
                                                                alreadyExist);
}

//___________________________________________________________________________

G4MolecularConfiguration*
G4MoleculeDefinition::GetConfigurationWithLabel(const G4String& molecularConfLabel)
{
  return G4MolecularConfiguration::GetMolecularConfiguration(this,
                                                             molecularConfLabel);
}

//______________________________________________________________________________

G4MolecularConfiguration*
G4MoleculeDefinition::
  NewConfigurationWithElectronOccupancy(const G4String& exStId,
                                        const G4ElectronOccupancy& elecConf,
                                        double decayTime)
{
  bool alreadyExist(false);
  G4MolecularConfiguration* conf =
      G4MolecularConfiguration::CreateMolecularConfiguration(GetName()+"_"+
                                                             exStId,
                                                             this,
                                                             exStId,
                                                             elecConf,
                                                             alreadyExist);

  conf->SetDecayTime(decayTime);

  return conf;
}

//______________________________________________________________________________

void G4MoleculeDefinition::SetLevelOccupation(G4int shell, G4int eNb)
{
  if(fElectronOccupancy)
  {
    G4int levelOccupancy = fElectronOccupancy->GetOccupancy(shell);

    if(levelOccupancy)
    {

      fElectronOccupancy->RemoveElectron(shell, levelOccupancy);
    }

    fElectronOccupancy->AddElectron(shell, eNb);
  }
}

//______________________________________________________________________________

void G4MoleculeDefinition::AddDecayChannel(const G4String& molecularConfLabel,
                                           const G4MolecularDissociationChannel* channel)
{
  if(!fDecayTable)
  {
    fDecayTable = new G4MolecularDissociationTable();
  }

  fDecayTable->AddChannel(G4MolecularConfiguration::GetMolecularConfiguration(this,
                                                                              molecularConfLabel),
                          channel);
}

//______________________________________________________________________________

void G4MoleculeDefinition::AddDecayChannel(const G4MolecularConfiguration* molConf,
                                           const G4MolecularDissociationChannel* channel)
{
  if (!fDecayTable)
  {
    fDecayTable = new G4MolecularDissociationTable();
  }
  fDecayTable->AddChannel(molConf, channel);
}

//___________________________________________________________________________

const vector<const G4MolecularDissociationChannel*>*
G4MoleculeDefinition::GetDecayChannels(const G4String& ExState) const
{
  if (fDecayTable)
  {
    const vector<const G4MolecularDissociationChannel*>* output =
        fDecayTable->GetDecayChannels(ExState);
    return output;
  }
  else
  {
    G4ExceptionDescription errMsg;
    errMsg <<  ": no Excited States and Decays for"
           << GetName()
           << " are defined.";
    G4Exception("G4MoleculeDefinition::GetDecayChannels", "",
                FatalErrorInArgument, errMsg);
  }
  return 0;
}

//___________________________________________________________________________

const vector<const G4MolecularDissociationChannel*>*
G4MoleculeDefinition::GetDecayChannels(const G4MolecularConfiguration* conf)
  const
{
  if(fDecayTable)
  {
    const vector<const G4MolecularDissociationChannel*>* output =
        fDecayTable->GetDecayChannels(conf);
    return output;
  }
//  else
//  {
//    G4ExceptionDescription errMsg;
//    errMsg << ": no Excited States and Decays for"
//           << GetName()
//           << " are defined.";
//    G4Exception("G4MoleculeDefinition::GetDecayChannels",
//                "",
//                FatalErrorInArgument,
//                errMsg);
//  }
  return 0;
}

//___________________________________________________________________________
// Protected
//___________________________________________________________________________

G4MoleculeDefinition::G4MoleculeDefinition(const G4MoleculeDefinition& right) :
    G4ParticleDefinition((const G4ParticleDefinition &) right),
    fDiffusionCoefficient(right.fDiffusionCoefficient),
    fAtomsNb(right.fAtomsNb),
    fVanDerVaalsRadius(right.fVanDerVaalsRadius)
{
  if (right.fElectronOccupancy != 0)
  {
    fElectronOccupancy = new G4ElectronOccupancy(*(right.fElectronOccupancy));
  }
  else fElectronOccupancy = 0;

  if (right.fDecayTable != 0)
  {
    fDecayTable = new G4MolecularDissociationTable(*(right.fDecayTable));
  }
  else fDecayTable = 0;

  fCharge = right.fCharge;
}

//___________________________________________________________________________

const G4MoleculeDefinition&
G4MoleculeDefinition::operator=(const G4MoleculeDefinition &right)
{
  if (this != &right)
  {
  }
  return *this;
}

//___________________________________________________________________________

void G4MoleculeDefinition::Finalize()
{
  G4MoleculeTable::Instance()->Finalize(this);
}

//___________________________________________________________________________


