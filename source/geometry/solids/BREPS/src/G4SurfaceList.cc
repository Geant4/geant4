// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceList.cc,v 1.4 2000-08-28 08:57:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4SurfaceList.cc
//
// ----------------------------------------------------------------------

#include "G4SurfaceList.hh"

G4SurfaceList::G4SurfaceList()
{
  first = index = last = next = temp = (G4Surface*)0;
  number_of_elements=0; 
}


G4SurfaceList::~G4SurfaceList()
{
  EmptyList();
}

void G4SurfaceList::MoveToFirst(G4Surface* srf)
{
  if(number_of_elements)
  {
    RemovePointer();
    srf->next = first;  
    first = srf;
    index=first;
    number_of_elements++;  
  }
}


void G4SurfaceList::AddSurface(G4Surface* srf)
{
  if(first == (G4Surface*)0)
  {
    index = srf;
    first = srf;
    last = srf;
    first->next = (G4Surface*)0;
  }
  else
  {
    srf->next = last->next;  
    last->next = srf;
    last = last->next;
  }
  
  number_of_elements++;  
  index=first;
}


G4Surface* G4SurfaceList::GetSurface()
{
  return index;
}


G4Surface* G4SurfaceList::GetSurface(G4int number)
{
  index = first;
  for(G4int a=0;a<number;a++)
    Step();
    
  return index;
}


G4Surface* G4SurfaceList::GetLastSurface()
{
  return last;
}


void G4SurfaceList::RemoveSurface(G4Surface* srf)
{
  if(srf!=(G4Surface*)0)
  {
    number_of_elements--;  
    temp = first;
    
    if(srf == first)
    {
      first=first->next;
      index = first;
      if(number_of_elements == 0)last = first;
      delete srf;
      return;
    }
    else	
    {
      while(temp->next != srf) temp = temp->next;
      index = srf->next;
      temp->next = index;
      if(srf == last) last = temp;
      index = first;
      delete srf;
    }
  }
}


void G4SurfaceList::RemovePointer()
{
  // Remove the current pointer from the List
  // Do not delete the object itself
  if(number_of_elements)
    if(first != index)
    {
      temp = first;
      
      // Find previous
      while(temp->next != index) temp = temp->next;
      
      // Hop over the one to be removed
      temp->next = index->next;
	
      // Correct the index pointer
      index = temp->next;
    }
    else
    {
      // Hop over the first
      first = first->next;
      index = first;
    }
  
  number_of_elements--;
}


void G4SurfaceList::EmptyList()
{
  //Deletes all surfaces in List
  while (first != (G4Surface*)0)
  {
    temp  = first;
    first = first->next;
    delete temp;
    number_of_elements--;
  }
  
  last = index = first;
}


void G4SurfaceList::MoveToFirst()
{
  index = first;
}


void G4SurfaceList::Step()
{
  if(index!=(G4Surface*)0)
    index = index->next;
}


void G4SurfaceList::G4SortList()
{
  if(number_of_elements == 1) return;
  
  // First create a vector of the surface distances
  // to the ray origin
  G4Surface** distances = new G4Surface*[number_of_elements];
  G4int x = 0;
  MoveToFirst();

  // Copy surface pointers to vector
  if(number_of_elements > 1)
  {
    while(x < number_of_elements)
    {
      distances[x] = index;
      index = index->next;
      x++;
    }
    
    MoveToFirst();

    // Sort List of pointers using quick G4Sort
    QuickG4Sort( distances, 0, number_of_elements-1 );
    
    // Organize the linked List of surfaces according
    // to the quickG4Sorted List.
    x = 0;
    first = distances[x];
    last = first;
    x++;	

    while (x < number_of_elements)
    {
      last->next = distances[x];
      last = last ->next;
      x++;
    }
	
    last->next = (G4Surface*)0;
    MoveToFirst();
  }
  
  delete[] distances;
}


void G4SurfaceList::QuickG4Sort(G4Surface** Dist, G4int left, G4int right)
{
  register G4int i=left;
  register G4int j=right;
  
  G4Surface* elem1;
  G4Surface* elem2 = Dist[(left+right)/2];
  
  do
  {
    while ( (Dist[i]->GetDistance() < elem2->GetDistance())  &&  (i < right) ) 
      i++;
    
    while ( (elem2->GetDistance() < Dist[j]->GetDistance())  &&  (j > left))
      j--;

    if(i<=j)
    {
      elem1   = Dist[i];
      Dist[i] = Dist[j];
      Dist[j] = elem1;
      i++;
      j--;
    }
  } while (i<=j);
  
  if( left < j  ) 
    QuickG4Sort(Dist, left, j );

  if( i < right ) 
    QuickG4Sort(Dist, i, right);    
}
