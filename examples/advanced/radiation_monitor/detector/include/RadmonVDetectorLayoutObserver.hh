//
// File name:     RadmonVDetectorLayoutObserver.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLayoutObserver.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class of the observer-subjectmodel for the detector
//                layout
//

#ifndef   RADMONVDETECTORLAYOUTOBSERVER_HH
 #define  RADMONVDETECTORLAYOUTOBSERVER_HH
 
 class RadmonVDetectorLayoutSubject;
 
 class RadmonVDetectorLayoutObserver
 {
  public:
   virtual void                                 OnLayoutChange(void) = 0;

  protected:
   inline                                       RadmonVDetectorLayoutObserver();
   inline                                      ~RadmonVDetectorLayoutObserver();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorLayoutObserver(const RadmonVDetectorLayoutObserver & copy);
   RadmonVDetectorLayoutObserver &              operator=(const RadmonVDetectorLayoutObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonVDetectorLayoutObserver.icc"
#endif /* RADMONVDETECTORLAYOUTOBSERVER_HH */
