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
// File name:     RadmonGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorLayout.hh,v 1.2 2006-06-28 13:52:59 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
