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
/// \file runAndEvent/RE03/include/RE03WorkerInitialization.hh
/// \brief Definition of the RE03WorkerInitialization class
//
//
// $Id: RE03WorkerInitialization.hh 66780 2013-01-12 14:56:35Z gcosmo $
//

#ifndef RE03WorkerInitialization_h
#define RE03WorkerInitialization_h 1

#include "G4UserWorkerInitialization.hh"
#include "globals.hh"

class RE03WorkerInitialization : public G4UserWorkerInitialization
{
  public:
    RE03WorkerInitialization();    
    virtual ~RE03WorkerInitialization();

  public:
    virtual void WorkerStart() const;
 
};

#endif


