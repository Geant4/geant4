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

// ====================================================================
//
//   XXVisManager.hh
//   $Id: XXVisManager.hh,v 1.1 2002-04-29 20:43:55 asaim Exp $
//
// ====================================================================
#ifndef XX_VIS_MANAGER_H
#define XX_VIS_MANAGER_H

#include "G4VisManager.hh"

class XXVisManager : public G4VisManager {
public:
  XXVisManager();
  ~XXVisManager();

private:
  virtual void RegisterGraphicsSystems();

};

#endif
