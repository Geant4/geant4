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
// File name:     RadmonDataAnalysisLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisLayout.hh,v 1.2 2006-06-28 13:44:01 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class that describes layer informations
//

#ifndef   RADMONDATAANALYSISLAYOUT_HH
 #define  RADMONDATAANALYSISLAYOUT_HH

 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonDataAnalysisLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonDataAnalysisLayout();
   inline                                       RadmonDataAnalysisLayout(const RadmonDataAnalysisLayout & copy);
   inline                                      ~RadmonDataAnalysisLayout();

   inline RadmonDataAnalysisLayout &            operator=(const RadmonDataAnalysisLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline const G4String &                      GetType(void) const;
   inline void                                  SetType(const G4String & type);

   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     dataAnalysisLabel;
   G4String                                     dataAnalysisType;
 };

 // Inline implementations
 #include "RadmonDataAnalysisLayout.icc"
#endif /* RADMONDATAANALYSISLAYOUT_HH */
