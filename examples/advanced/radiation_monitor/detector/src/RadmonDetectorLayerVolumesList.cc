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
// File name:     RadmonDetectorLayerVolumesList.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.cc,v 1.3 2006/06/29 16:13:57 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

// Include files
#include "RadmonDetectorLayerVolumesList.hh"
#include "RadmonDetectorLayerVolumeItem.hh"

RadmonDetectorLayerVolumeItem *                 RadmonDetectorLayerVolumesList :: AppendItem(void)
{
 RadmonDetectorLayerVolumeItem * item(new RadmonDetectorLayerVolumeItem);
 volumesVector.push_back(item);
 
 return item;
}





void                                            RadmonDetectorLayerVolumesList :: RemoveItem(G4int index)
{
 VolumesVector::iterator i(volumesVector.begin());
 i+=index;
 
 RemoveRange(i, i+1);
}



void                                            RadmonDetectorLayerVolumesList :: RemoveItemsByRange(G4int first, G4int last)
{
 VolumesVector::iterator i(volumesVector.begin());
 
 RemoveRange(i+first, i+last);
}



void                                            RadmonDetectorLayerVolumesList :: RemoveAllItems(void)
{
 RemoveRange(volumesVector.begin(), volumesVector.end());
}

void                                            RadmonDetectorLayerVolumesList :: RemoveRange(VolumesVector::iterator begin, VolumesVector::iterator end)
{
 VolumesVector::iterator i;
 VolumesVector::iterator j;
 const VolumesVector::iterator endJ(volumesVector.end());
 
 G4bool kill;
 G4bool somethingDone(true);
 G4bool stillSomethingToDo;

 do
 {
  if (!somethingDone)
   G4Exception("RadmonDetectorLayerVolumesList::RemoveAllItems: Some dependencies cannot be removed.");
 
  stillSomethingToDo=false;
  somethingDone=false;
 
  i=begin;
  while (i!=end)
  {
   if (*i)
   {
    stillSomethingToDo=true;
    
    kill=true;
    j=volumesVector.begin();
    
    while (j!=endJ)
    {
     if (*j)
      if ((*j)->GetMotherVolumeItem()==(*i))
      {
       kill=false;
       break;
      }

     j++;
    }
   }
   else
    kill=false;
 
   if (kill)
   {
    delete (*i);
    (*i)=0;
    somethingDone=true;
   }
   
   i++;
  }
 }
 while (stillSomethingToDo);  
 
 volumesVector.erase(begin, end);
}
