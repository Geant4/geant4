//
// File name:     RadmonDetectorPadsData.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorPadsData.hh,v 1.1 2005-09-30 08:33:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class for pads data managemenet
//

#ifndef   RADMONDETECTORPADSDATA_HH
 #define  RADMONDETECTORPADSDATA_HH
 
 // Include files
 #include "G4ThreeVector.hh" 
 #include "G4RotationMatrix.hh"
 #include <vector>
 
 // Forward declarations
 class G4Material;
 
 class RadmonDetectorPadsData
 {
  public:
   inline                                       RadmonDetectorPadsData();
                                                RadmonDetectorPadsData(const RadmonDetectorPadsData & copy);
   inline                                      ~RadmonDetectorPadsData();
   
   inline G4double                              GetWidth(void) const;
   inline G4double                              GetHeight(void) const;
   inline G4Material *                          GetMaterial(void) const;
   G4int                                        GetNPads(void) const;
   const G4ThreeVector &                        GetPosition(G4int index) const;
   const G4RotationMatrix &                     GetRotation(G4int index) const;

   inline void                                  SetWidth(G4double width);
   inline void                                  SetHeight(G4double height);
   inline void                                  SetMaterial(G4Material * material);
   inline void                                  AppendPosition(const G4ThreeVector & position);
   void                                         AppendPositionAndRotation(const G4ThreeVector & position, const G4RotationMatrix & rotation);
   bool                                         ReadPositionsAndRotationsFromString(const G4String & positionsStr);

   RadmonDetectorPadsData &                     operator=(const RadmonDetectorPadsData & copy);

  private:
   bool                                         ProcessElement(G4double & x, G4double & y, G4double & delta, const G4String & content);
   bool                                         ReadUmis(G4double & value, const G4String & text, const char * category);

   
  // Private data types
   typedef std::vector<G4ThreeVector>           PadPositions;
   typedef std::vector<G4RotationMatrix>        PadRotations;
  
  // Private attributes
   G4double                                     padWidth;
   G4double                                     padHeight;
   G4Material *                                 padMaterial;
   PadPositions                                 padPositions;
   PadRotations                                 padRotations;
 };
 
 // Inline implementations
 #include "RadmonDetectorPadsData.icc"
#endif /* RADMONDETECTORPADSDATA_HH */
