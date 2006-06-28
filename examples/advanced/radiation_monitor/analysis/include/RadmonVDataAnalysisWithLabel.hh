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
// File name:     RadmonVDataAnalysisWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysisWithLabel.hh,v 1.2 2006-06-28 13:44:23 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of an analysis piece with label
//                and attributes
//

#ifndef   RADMONVDATAANALYSISWITHLABEL_HH
 #define  RADMONVDATAANALYSISWITHLABEL_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 #include "RadmonVDataAnalysis.hh"
 
 class RadmonVDataAnalysisWithLabel : public RadmonVDataAnalysis, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVDataAnalysisWithLabel();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetDataAnalysisAttribute(const G4String & attributeName, const G4String & value);

   virtual RadmonVDataAnalysisWithLabel *       New(void) const = 0;

  protected:
   inline                                       RadmonVDataAnalysisWithLabel(const G4String & label);
   
  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysisWithLabel();
                                                RadmonVDataAnalysisWithLabel(const RadmonVDataAnalysisWithLabel & copy);
   RadmonVDataAnalysisWithLabel &               operator=(const RadmonVDataAnalysisWithLabel & copy);

  // Private attributes
   G4String                                     physiscListLabel;
 };
 
 #include "RadmonVDataAnalysisWithLabel.icc"
#endif /* RADMONVDATAANALYSISWITHLABEL_HH */
