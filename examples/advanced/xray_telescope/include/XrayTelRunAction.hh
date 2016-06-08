//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelRunAction.hh                             *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              16/10/01                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 07.11.2001 M.G. Pia
// - Modified the analysis management
//
// 16.10.2001 R.Nartallo
// - Clean up code to avoid 'pedantic' and 'ANSI' compiler warnings 
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelRunAction_h
#define XrayTelRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class XrayTelRunAction : public G4UserRunAction
{
public:
  XrayTelRunAction();

  ~XrayTelRunAction();

  void BeginOfRunAction(const G4Run*);

  void EndOfRunAction(const G4Run*);

  void Update(G4double energy);

private:
 
  G4int nEnteringTracks;
  G4double totEnteringEnergy;

};

#endif

