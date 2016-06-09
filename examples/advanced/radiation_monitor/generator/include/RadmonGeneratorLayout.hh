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
// File name:     RadmonGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorLayout.hh,v 1.3 2006/06/29 16:15:09 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Implementation of the particles source data
//

#ifndef   RADMONGENERATORLAYOUT_HH
 #define  RADMONGENERATORLAYOUT_HH

 // Include files
 #include "RadmonVGeneratorLayout.hh"
 #include "RadmonTLabelledCollection.hh"
 #include "RadmonGeneratorSourceLayout.hh"
 
 class RadmonGeneratorLayout : public RadmonVGeneratorLayout
 {
  public:
   inline                                     RadmonGeneratorLayout();
   inline                                    ~RadmonGeneratorLayout();

   virtual void                               InsertSource(const G4String & sourceLabel);
   virtual void                               SetRelativeSourceIntensity(const G4String & sourceLabel, G4double relativeIntensity);
   virtual G4double                           GetRelativeSourceIntensity(const G4String & sourceLabel) const;
   virtual void                               RemoveSource(const G4String & sourceLabel);
   virtual G4int                              GetNSources(void) const;
   virtual const G4String &                   GetSourceLabel(G4int index) const;
   virtual void                               AppendSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel);
   virtual void                               SetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & typeName);
   virtual void                               RemoveSourceAlgorithm(const G4String & sourceLabel, const G4String & algorithmLabel);
   virtual G4int                              GetNSourceAlgorithms(const G4String & sourceLabel) const;
   virtual const G4String &                   GetSourceAlgorithmLabel(const G4String & sourceLabel, G4int index) const;
   virtual const G4String &                   GetSourceAlgorithmType(const G4String & sourceLabel, const G4String & algorithmLabel) const;
   virtual void                               SetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & value);
   virtual void                               ClearSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute);
   virtual G4String                           GetSourceAlgorithmAttribute(const G4String & sourceLabel, const G4String & algorithmLabel, const G4String & attribute, const G4String & defaultValue=G4String()) const;
   virtual G4int                              GetSourceAlgorithmNAttributes(const G4String & sourceLabel, const G4String & algorithmLabel) const;
   virtual const G4String &                   GetSourceAlgorithmAttributeName(const G4String & sourceLabel, const G4String & algorithmLabel, G4int index) const;
   virtual G4bool                             Load(std::istream & in);
   virtual G4bool                             Save(std::ostream & out) const;
   virtual void                               DumpLayout(std::ostream & out) const;

  private:
   inline G4String &                          GetNullStr() const;

  // Hidden constructors and operators
                                              RadmonGeneratorLayout(const RadmonGeneratorLayout & copy);
   RadmonGeneratorLayout &                    operator=(const RadmonGeneratorLayout & copy);
   
  // Private attributes
   RadmonTLabelledCollection<RadmonGeneratorSourceLayout> labelledSourcesCollection;
 };

 // Inline implementations
 #include "RadmonGeneratorLayout.icc"
#endif /* RADMONGENERATORLAYOUT_HH */
