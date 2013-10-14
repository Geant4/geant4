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
//      LaunchG4.hh
//
//      Copyright 2010 M. Karamitros <kara@cenbg.in2p3.fr>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#ifndef LaunchG4_HH
#define LaunchG4_HH 1

#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIExecutive;
class G4VisManager;
class G4RunManager;

class LaunchG4
{
public:
    LaunchG4();
    ~LaunchG4();

    void Initialize(G4bool chemistryflag = true);
    void RunSimu(G4String macFile);
    void StartSession();

    void NewSession(int argc,char** argv, const G4String& sessionType = "");
    void BuildGUIFrame();

private :
    PrimaryGeneratorAction* fpPrimGenAct ;

    G4RunManager * fpRunManager ;

#ifdef G4UI_USE
    G4UIExecutive * fpSession ;
#endif

#ifdef G4VIS_USE
    G4VisManager* fpVisManager ;
#endif

};
#endif /* LaunchG4_HH */
