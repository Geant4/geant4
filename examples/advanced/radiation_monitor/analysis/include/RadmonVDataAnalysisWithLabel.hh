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
// File name:     RadmonVDataAnalysisWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysisWithLabel.hh,v 1.3 2006/06/29 16:07:23 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
