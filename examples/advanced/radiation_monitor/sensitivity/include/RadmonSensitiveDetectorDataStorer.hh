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
// File name:     RadmonSensitiveDetectorDataStorer.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetectorDataStorer.hh,v 1.2 2006-06-28 13:57:05 gunter Exp $
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
