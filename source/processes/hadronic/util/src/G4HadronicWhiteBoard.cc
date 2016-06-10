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

#include "G4HadronicWhiteBoard.hh"

G4ThreadLocal G4HadronicWhiteBoard * G4HadronicWhiteBoard::theInstance = 0;

G4HadronicWhiteBoard::G4HadronicWhiteBoard()
 : theProjectile(0), theDef(0), theName(0), theE(0.0), thePx(0.0),
   thePy(0.0), thePz(0.0), theA(0), theZ(0)
{}


G4HadronicWhiteBoard & G4HadronicWhiteBoard::Instance()
{
  if (!theInstance) theInstance = new G4HadronicWhiteBoard;
  return *theInstance;
}

const G4HadProjectile * G4HadronicWhiteBoard::GetProjectile() 
{
  return theProjectile;
}

const G4Nucleus & G4HadronicWhiteBoard::GetTargetNucleus() 
{ 
  return theTarget;
}

const G4ParticleDefinition * G4HadronicWhiteBoard::GetPDef() {return theDef;}
G4String G4HadronicWhiteBoard::GetParticleName() {return theName;}
G4double G4HadronicWhiteBoard::GetEnergy() {return theE;}
G4double G4HadronicWhiteBoard::GetPx(){return thePx;}
G4double G4HadronicWhiteBoard::GetPy(){return thePy;}
G4double G4HadronicWhiteBoard::GetPz(){return thePz;}
G4int G4HadronicWhiteBoard::GetA(){return theA;}
G4int G4HadronicWhiteBoard::GetZ(){return theZ;}

  
void G4HadronicWhiteBoard::SetProjectile(const G4HadProjectile & aProjectile)
{
  theProjectile = & aProjectile;
  theDef = theProjectile->GetDefinition();
  theName = theDef->GetParticleName().c_str();
  theE = theProjectile->Get4Momentum().t();
  thePx = theProjectile->Get4Momentum().vect().x();
  thePy = theProjectile->Get4Momentum().vect().y();
  thePz = theProjectile->Get4Momentum().vect().z();
}

    
void G4HadronicWhiteBoard::SetTargetNucleus(const G4Nucleus & aTarget) 
{
  theTarget = aTarget;
  theA = theTarget.GetA_asInt();
  theZ = theTarget.GetZ_asInt();
}


void G4HadronicWhiteBoard::SetProcessName(const G4String& aProcessName) 
{
  theProcessName = aProcessName;
}


void G4HadronicWhiteBoard::SetModelName(const G4String& aModelName) 
{
  theModelName = aModelName;
}


void G4HadronicWhiteBoard::Dump()
{
  std::cerr << std::endl;
  std::cerr << "*** Geant4 Hadronic Reaction Information ***" 
            << std::endl;
  std::cerr << "    Process: " << theProcessName << " , Model: " 
            << theModelName << std::endl;
  std::cerr << "    Nucleus A, Z = " << theA << " " << theZ 
            << std::endl;
  std::cerr << "    Projectile was a " << theName 
            << std::endl;
  std::cerr << "    projectile momentum (px, py, pz) = (" << thePx << ", " 
            << thePy << ", " << thePz << ")" << std::endl;
  std::cerr << "    Projectile energy = "<< theE 
            << std::endl;
  std::cerr << "*** End of Geant4 Hadronic Reaction Information ***"
            << std::endl;
  G4Exception(theModelName, "001", FatalException, "segmentation fault");
}
