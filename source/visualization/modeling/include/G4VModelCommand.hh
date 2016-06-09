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
// $Id: G4VModelCommand.hh,v 1.3 2006/06/29 21:32:27 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
// 
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description
// Templated base class for model messengers. Commands specific to a particular
// concrete model should inherit from G4VModelCommand, with the concrete model
// type as the template parameter.
// Class Description - End:

#ifndef G4VMODELCOMMAND_HH
#define G4VMODELCOMMAND_HH

#include "G4UImessenger.hh"
#include "G4String.hh"

class G4UIcommand;

template <typename T>
class G4VModelCommand : public G4UImessenger {

public: // With description

  G4VModelCommand(T* model);
  // Input model

  virtual ~G4VModelCommand();

  G4String GetCurrentValue(G4UIcommand* command);

protected:

  T* Model();
  // Access to model

private:

  T* fpModel;

};

template <typename T>
G4VModelCommand<T>::G4VModelCommand(T* model)
  :fpModel(model)
{}

template <typename T>
G4VModelCommand<T>::~G4VModelCommand() {}

template <typename T>
G4String 
G4VModelCommand<T>::GetCurrentValue(G4UIcommand*) 
{
  return "";
}

template <typename T>
T*
G4VModelCommand<T>::Model() 
{
  return fpModel;
}

#endif



