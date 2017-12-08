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
#include "G4RNGHelper.hh"
#include "Randomize.hh"

template<>
G4TemplateRNGHelper<G4long>* G4TemplateRNGHelper<G4long>::instance = 0;

template<>
G4TemplateRNGHelper<G4String>* G4TemplateRNGHelper<G4String>::instance = 0;

template<class T>
G4TemplateRNGHelper<T>* G4TemplateRNGHelper<T>::GetInstance()
{
  if (!instance)
  {
    instance = new G4TemplateRNGHelper<T>();
  }
  return instance;
}
template<class T>
G4TemplateRNGHelper<T>* G4TemplateRNGHelper<T>::GetInstanceIfExist()
{
  return instance;
}

template<>
G4TemplateRNGHelper<G4long>* G4TemplateRNGHelper<G4long>::GetInstance()
{
  if (!instance)
  {
    instance = new G4TemplateRNGHelper<G4long>();
  }
  return instance;
}
template<>
G4TemplateRNGHelper<G4long>* G4TemplateRNGHelper<G4long>::GetInstanceIfExist()
{
  return instance;
}

template<>
G4TemplateRNGHelper<G4String>* G4TemplateRNGHelper<G4String>::GetInstance()
{
  if (!instance)
  {
    instance = new G4TemplateRNGHelper<G4String>();
  }
  return instance;
}
template<>
G4TemplateRNGHelper<G4String>* G4TemplateRNGHelper<G4String>::GetInstanceIfExist()
{
  return instance;
}

template<class T>
G4TemplateRNGHelper<T>::~G4TemplateRNGHelper()
{
  Clear();
  instance = 0;
}
