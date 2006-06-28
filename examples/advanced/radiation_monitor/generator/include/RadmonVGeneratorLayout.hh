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
// File name:     RadmonVGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorLayout.hh,v 1.2 2006-06-28 13:53:31 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
