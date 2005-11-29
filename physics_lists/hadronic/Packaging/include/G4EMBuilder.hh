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
// $Id: G4EMBuilder.hh,v 1.3 2005-11-29 16:53:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EMBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
//
//----------------------------------------------------------------------------
//

#ifndef G4EMBuilder_h
#define G4EMBuilder_h 1

#include "G4EMTailorer.hh"
#include "globals.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4ElectroNuclearBuilder.hh"

class G4EMBuilder
{
public:
  G4EMBuilder();
  virtual ~G4EMBuilder();

  void Build();

  void Synch(G4String & aState);
  void GammaNuclear(G4String & aState);

private:
  G4bool wasActivated;
  G4bool synchOn;
  G4bool gammNucOn;

  G4EMTailorer * theT;
  G4SynchrotronRadiation*  theElectronSynch;
  G4SynchrotronRadiation*  thePositronSynch;
  G4ElectroNuclearBuilder* theGNPhysics;
};

#endif





