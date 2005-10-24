//
// File name:     RadmonVLayoutSubject.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVLayoutSubject.cc,v 1.1 2005-10-24 14:51:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonVLayoutSubject.hh"

#include "RadmonVLayoutObserver.hh"


void                                            RadmonVLayoutSubject :: AttachObserver(RadmonVLayoutObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonVLayoutSubject :: DetachObserver(RadmonVLayoutObserver * observer)
{
 observersSet.erase(observer);
}





void                                            RadmonVLayoutSubject :: NotifyChange(void)
{
 ObserversSet::iterator i=observersSet.begin();
 ObserversSet::iterator end=observersSet.end();

 while (i!=end)
 {
  (*i)->OnLayoutChange();
  i++;
 }
}
