//
// File name:     RadmonGeneratorsWithLabelFactory.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.cc,v 1.1 2005-10-25 16:36:43 capra Exp $
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
