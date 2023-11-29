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

#include <sstream>

#include "G3EleTable.hh"

#include "G4Types.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G3EleTable::G3EleTable()
{
  // create an array of pointers to elements
  _Ele = new G4Element*[_MaxEle];
  LoadUp();
}

G3EleTable::~G3EleTable()
{
  delete [] _Ele;
}

G4Element* 
G3EleTable::GetEle(G4double Z)
{
  G4double A = 0.;
  char name[20], sym[3];
  G4int index = (G4int) Z-1;
  if (!parse(Z, name, sym, A))
  {
    G4String na(name);
    G4String sy(sym);
    if (_Ele[index] == nullptr)
    {
      if ( A == 0. )
      {
        // Cannot create element with A = 0.
        G4String text = "Failed to get element Z = " + std::to_string(Z);
        G4Exception("G3EleTable::GetEle", "G3toG40016", FatalException, text);
      }
      // add an element to the element table here
      _Ele[index] = new G4Element(na, sy, Z, A*g/mole);
    }
  }
  return _Ele[index];
}

G4int 
G3EleTable::parse(G4double& Z, char* name, char* sym, G4double& A)
{ 
  G4int rc = 0;
  if (Z>0 && Z <=_MaxEle)
  {
    G4int z = (G4int) Z-1;
    G4String str(_EleNames[z]);
    char* cstr = new char [str.length()+1];
    std::strcpy(cstr, str.c_str());
    char* p = std::strtok(cstr," ");
    std::strcpy(name, p);
    p = std::strtok(NULL," ");
    std::strcpy(sym, p);
    p = std::strtok(NULL," ");
    std::istringstream in(p);
    in >> A;
    delete [] cstr;
  }
  else
  {
    rc = -1;
  }
  return rc;
}

void
G3EleTable::LoadUp()
{
  // initialize element pointers to 0
  for (G4int j=0; j<_MaxEle; ++j)
  {
    _Ele[j]=nullptr;
  }
}
