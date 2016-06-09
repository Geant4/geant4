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
// File name:     RadmonVGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorLayout.hh,v 1.3 2006/06/29 16:15:54 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Abstract class for storing primary generator data
//

#ifndef   RADMONVGENERATORLAYOUT_HH
 #define  RADMONVGENERATORLAYOUT_HH
 
 // Include files
 #include "RadmonVLayoutSubject.hh"
 #include "G4String.hh"
 #include <iostream>

 class RadmonVGeneratorLayout : public RadmonVLayoutSubject
 {
  public:
   virtual void                               InsertSource(const G4String & sourceLabel) = 0;
   virtual void                               SetRelativeSourceIntensity(const G4String & sourceLabel, G4double relativeIntensity) = 0;
   virtual G4double                           GetRelativeSourceIntensity(const G4String & sourceLabel) const = 0;
   virtual void                               RemoveSource(const G4String & sourceLabel) = 0;
   virtual G4int                              GetNSources(void) const = 0;
   virtual const G4String &                   GetSourceLabel(G4int index) const = 0;
   virtual void                               AppendSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel) = 0;
   virtual void                               SetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & typeName) = 0;
   virtual void                               RemoveSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel) = 0;
   virtual G4int                              GetNSourceAlgorithms(const G4String & sourceLabel) const = 0;
   virtual const G4String &                   GetSourceAlgorithmLabel(const G4String & sourceLabel, G4int index) const = 0;
   virtual const G4String &                   GetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel) const = 0;
   virtual void                               SetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & value) = 0;
   virtual void                               ClearSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute) = 0;
   virtual G4String                           GetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & defaultValue=G4String()) const = 0;
   virtual G4int                              GetSourceAlgorithmNAttributes(const G4String & sourceLabel, const G4String & algorithmLabel) const = 0;
   virtual const G4String &                   GetSourceAlgorithmAttributeName(const G4String & sourceLabel, const G4String & algorithmLabel, G4int index) const = 0;
   virtual G4bool                             Load(std::istream & in) = 0;
   virtual G4bool                             Save(std::ostream & out) const = 0;
   virtual void                               DumpLayout(std::ostream & out) const = 0;

  protected:
   inline                                     RadmonVGeneratorLayout();
   inline                                    ~RadmonVGeneratorLayout();

  // Hidden constructors and operators
  private:
                                              RadmonVGeneratorLayout(RadmonVGeneratorLayout & copy);
   RadmonVGeneratorLayout &                   operator=(const RadmonVGeneratorLayout & copy);
 };
 
 // Inline implementations
 #include "RadmonVGeneratorLayout.icc"
#endif /* RADMONVGENERATORLAYOUT_HH */
