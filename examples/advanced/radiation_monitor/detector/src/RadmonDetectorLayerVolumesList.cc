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
// File name:     RadmonDetectorLayerVolumesList.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.cc,v 1.2 2006-06-28 13:50:36 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
