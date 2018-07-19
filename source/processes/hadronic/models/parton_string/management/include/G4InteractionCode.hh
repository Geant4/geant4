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
// $Id: G4InteractionCode.hh 100828 2016-11-02 15:25:59Z gcosmo $
//
#ifndef G4InteractionCode_h
#define G4InteractionCode_h 1


class G4InteractionCode
{
  public:
    G4InteractionCode(G4String & aCode);

    void SetCode(G4String & aCode);
    G4String GetCode();

  private:
    G4String theCode;
};

inline void G4InteractionCode::SetCode(G4String & aCode)
{
  theCode = aCode;
}

inline G4String G4InteractionCode::GetCode()
{
  return theCode;
}

#endif

