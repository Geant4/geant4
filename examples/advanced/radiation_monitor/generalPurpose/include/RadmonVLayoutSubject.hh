//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonVLayoutSubject.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVLayoutSubject.hh,v 1.2.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-02 $
//
// Description:   Subject class of the observer-subjectmodel for the detector
//                layout
//

#ifndef   RADMONVLAYOUTSUBJECT_HH
 #define  RADMONVLAYOUTSUBJECT_HH
 
 // Include files
 #include <set>
 
 // Forward declarations
 class RadmonVLayoutObserver;
 
 class RadmonVLayoutSubject
 {
  public:
   void                                         AttachObserver(RadmonVLayoutObserver * observer);
   void                                         DetachObserver(RadmonVLayoutObserver * observer);

  protected:
   inline                                       RadmonVLayoutSubject();
   inline                                      ~RadmonVLayoutSubject();

   void                                         NotifyChange(void);

  private:
  // Hidden constructors and operators
                                                RadmonVLayoutSubject(const RadmonVLayoutSubject & copy);
   RadmonVLayoutSubject &                       operator=(const RadmonVLayoutSubject & copy);

  // Private attributes
  // List of pointers to RadmonVLayoutObserver
  // set does not care of the order of the elements inside the list
   typedef std::set<RadmonVLayoutObserver *>    ObserversSet;
   ObserversSet                                 observersSet;
 };

 // Inline implementations
 #include "RadmonVLayoutSubject.icc"
#endif /* RADMONVLAYOUTSUBJECT_HH */
