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
// File name:     RadmonApplicationRunNumbering.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationRunNumbering.hh,v 1.2 2006-06-28 13:45:38 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Dumps run id
//

#ifndef   RADMONAPPLICATIONRUNNUMBERING_HH
 #define  RADMONAPPLICATIONRUNNUMBERING_HH
 
 // Include files
 #include "RadmonRunActionObserver.hh"
 #include "globals.hh"
 
 class RadmonApplicationRunNumbering : public RadmonRunActionObserver
 {
  public:
   inline                                       RadmonApplicationRunNumbering();
   inline virtual                              ~RadmonApplicationRunNumbering();
   
   inline void                                  Enable(void);
   inline void                                  Disable(void);
   
   virtual void                                 OnBeginOfRun(const G4Run * run);
   inline virtual void                          OnEndOfRun(const G4Run * run);
   
  private:
                                                RadmonApplicationRunNumbering(const RadmonApplicationRunNumbering & copy);
   RadmonApplicationRunNumbering &              operator=(const RadmonApplicationRunNumbering & copy);
   
   G4bool                                       enable;
 };
 
 // Inline implementations
 #include "RadmonApplicationRunNumbering.icc"
#endif /* RADMONAPPLICATIONRUNNUMBERING_HH */
