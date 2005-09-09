//
// File name:     RadmonVDetectorLayoutSubject.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayoutSubject.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonVDetectorLayoutSubject.hh"

#include "RadmonVDetectorLayoutObserver.hh"


void                                            RadmonVDetectorLayoutSubject :: AttachObserver(RadmonVDetectorLayoutObserver * observer)
{
 observersSet.insert(observer);
}



void                                            RadmonVDetectorLayoutSubject :: DetachObserver(RadmonVDetectorLayoutObserver * observer)
{
 observersSet.erase(observer);
}





void                                            RadmonVDetectorLayoutSubject :: NotifyChange(void)
{
 ObserversSet::iterator i=observersSet.begin();
 ObserversSet::iterator end=observersSet.end();

 while (i!=end)
 {
  (*i)->OnLayoutChange();
  i++;
 }
}
