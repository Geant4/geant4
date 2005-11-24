//
// File name:     RadmonSensitiveDetectorDataStorer.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetectorDataStorer.hh,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Stores data step into a radmon hit
//

#ifndef   RADMONSENSITIVEDETECTORDATASTORER_HH
 #define  RADMONSENSITIVEDETECTORDATASTORER_HH
 
 // Include files
 
 class RadmonHit;
 class G4Step;
 
 class RadmonSensitiveDetectorDataStorer
 {
  public:
   virtual void                                 StoreIntoHit(G4Step * step, RadmonHit * hit) = 0;
   
  protected:
   inline                                       RadmonSensitiveDetectorDataStorer();
   inline                                      ~RadmonSensitiveDetectorDataStorer();

  private:
  // Hidden constructors and operators
                                                RadmonSensitiveDetectorDataStorer(const RadmonSensitiveDetectorDataStorer & copy);
   RadmonSensitiveDetectorDataStorer &                                  operator=(const RadmonSensitiveDetectorDataStorer & copy);
 };
 
 // Inline implementations
 #include "RadmonSensitiveDetectorDataStorer.icc"
#endif /* RADMONSENSITIVEDETECTORDATASTORER_HH */
