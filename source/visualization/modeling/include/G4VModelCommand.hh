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
// $Id: G4VModelCommand.hh,v 1.2 2005/11/28 20:07:11 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
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



