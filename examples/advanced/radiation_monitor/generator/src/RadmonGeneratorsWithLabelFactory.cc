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
// File name:     RadmonGeneratorsWithLabelFactory.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.cc,v 1.2 2006-06-28 13:54:03 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
