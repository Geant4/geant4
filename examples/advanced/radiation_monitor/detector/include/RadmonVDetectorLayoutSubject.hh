//
// File name:     RadmonVDetectorLayoutSubject.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayoutSubject.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Subject class of the observer-subjectmodel for the detector
//                layout
//

#ifndef   RADMONVDETECTORLAYOUTSUBJECT_HH
 #define  RADMONVDETECTORLAYOUTSUBJECT_HH
 
 // Include files
 #include <set>
 
 // Forward declarations
 class RadmonVDetectorLayoutObserver;
 
 class RadmonVDetectorLayoutSubject
 {
  public:
   void                                         AttachObserver(RadmonVDetectorLayoutObserver * observer);
   void                                         DetachObserver(RadmonVDetectorLayoutObserver * observer);

  protected:
   inline                                       RadmonVDetectorLayoutSubject();
   inline                                      ~RadmonVDetectorLayoutSubject();

   void                                         NotifyChange(void);

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorLayoutSubject(const RadmonVDetectorLayoutSubject & copy);
   RadmonVDetectorLayoutSubject &               operator=(const RadmonVDetectorLayoutSubject & copy);

  // Private attributes
   typedef std::set<RadmonVDetectorLayoutObserver *> ObserversSet;
   ObserversSet                                 observersSet;
 };

 // Inline implementations
 #include "RadmonVDetectorLayoutSubject.icc"
#endif /* RADMONVDETECTORLAYOUTSUBJECT_HH */
