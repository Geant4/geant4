//
// File name:     RadmonGeneratorSourceAlgorithmLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceAlgorithmLayout.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Container for the data of an algorithm of a particles source
//

#ifndef   RADMONGENERATORSOURCEALGORITHMLAYOUT_HH
 #define  RADMONGENERATORSOURCEALGORITHMLAYOUT_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonGeneratorSourceAlgorithmLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonGeneratorSourceAlgorithmLayout();
                                                RadmonGeneratorSourceAlgorithmLayout(const RadmonGeneratorSourceAlgorithmLayout & copy);
   inline                                      ~RadmonGeneratorSourceAlgorithmLayout();

   RadmonGeneratorSourceAlgorithmLayout &       operator=(const RadmonGeneratorSourceAlgorithmLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline const G4String &                      GetType(void) const;
   inline void                                  SetType(const G4String & type);

   void                                         DumpLayout(std::ostream & out, const G4String & indent = G4String()) const;

  // Private attributes
  private:
   G4String                                     algorithmLabel;
   G4String                                     algorithmType;
 };
 
 // Inline implementations
 #include "RadmonGeneratorSourceAlgorithmLayout.icc"
#endif /* RADMONGENERATORSOURCEALGORITHMLAYOUT_HH */
