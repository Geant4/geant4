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
// File name:     RadmonGeneratorsWithLabelFactory.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.cc,v 1.3 2006/06/29 16:16:37 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//

// Include files
#include "RadmonGeneratorsWithLabelFactory.hh"
#include "RadmonVGeneratorWithLabel.hh"
#include "RadmonGeneratorsWithLabelMessenger.hh"



                                                RadmonGeneratorsWithLabelFactory :: RadmonGeneratorsWithLabelFactory()
{
}



                                                RadmonGeneratorsWithLabelFactory :: ~RadmonGeneratorsWithLabelFactory()
{
 GeneratorsList::iterator i(generatorsList.begin());
 const GeneratorsList::iterator end(generatorsList.end());
 
 RadmonGeneratorsWithLabelMessenger * messenger(RadmonGeneratorsWithLabelMessenger::Instance());

 while (i!=end)
 {
  messenger->RemoveAvailableGenerator((*i)->GetLabel());
  
  delete (*i);

  i++;
 }
}





RadmonVGenerator *                              RadmonGeneratorsWithLabelFactory :: GetGenerator(const G4String & generatorType)
{
 GeneratorsList::const_iterator i(generatorsList.begin());
 const GeneratorsList::const_iterator end(generatorsList.end());
 
 while (i!=end)
 {
  if ((*i)->GetLabel()==generatorType)
   return (*i)->New();

  i++;
 }

 return 0; 
}





void                                            RadmonGeneratorsWithLabelFactory :: AppendGenerator(RadmonVGeneratorWithLabel * generator)
{
 generatorsList.push_back(generator);

 RadmonGeneratorsWithLabelMessenger * messenger(RadmonGeneratorsWithLabelMessenger::Instance());
 messenger->AddAvailableGenerator(generator->GetLabel());
}
