//
// File name:     RadmonVLayoutObserver.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVLayoutObserver.hh,v 1.1 2005-10-24 14:51:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Observer class of the observer-subjectmodel for the detector
//                layout
//

#ifndef   RADMONVLAYOUTOBSERVER_HH
 #define  RADMONVLAYOUTOBSERVER_HH
 
 class RadmonVLayoutSubject;
 
 class RadmonVLayoutObserver
 {
  public:
   virtual void                                 OnLayoutChange(void) = 0;

  protected:
   inline                                       RadmonVLayoutObserver();
   inline                                      ~RadmonVLayoutObserver();

  private:
  // Hidden constructors and operators
                                                RadmonVLayoutObserver(const RadmonVLayoutObserver & copy);
   RadmonVLayoutObserver &                      operator=(const RadmonVLayoutObserver & copy);
 };
 
 // Inline implementations
 #include "RadmonVLayoutObserver.icc"
#endif /* RADMONVLAYOUTOBSERVER_HH */
