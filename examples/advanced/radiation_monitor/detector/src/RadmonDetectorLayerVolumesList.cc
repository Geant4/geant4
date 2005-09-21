//
// File name:     RadmonDetectorLayerVolumesList.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumesList.cc,v 1.1 2005-09-21 14:52:32 capra Exp $
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
