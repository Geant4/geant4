//
// File name:     RadmonPrimaryGeneratorAction.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPrimaryGeneratorAction.hh,v 1.2 2006-03-08 17:55:37 sfonseca Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Primary particles generator user action
//

#ifndef   RADMONPRIMARYGENERATORACTION_HH
 #define  RADMONPRIMARYGENERATORACTION_HH
 
 // Include files
 #include "RadmonVLayoutObserver.hh"
 #include "G4VUserPrimaryGeneratorAction.hh"
 #include "G4ParticleGun.hh"
 #include "G4String.hh"
 #include "globals.hh"
 #include <map>
 #include "G4Event.hh"
 #include "G4GeneralParticleSource.hh"

//**
class G4GeneralParticleSource;
class G4Event;

 // Forward declarations
 class RadmonVGeneratorsFactory;
 class RadmonVGeneratorLayout;
 class RadmonVGenerator;

 class RadmonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public RadmonVLayoutObserver
 {
  public:
                                                RadmonPrimaryGeneratorAction(RadmonVGeneratorLayout * layout, RadmonVGeneratorsFactory * factory);
                                               ~RadmonPrimaryGeneratorAction();

   virtual void                                 GeneratePrimaries(G4Event * anEvent);

   virtual void                                 OnLayoutChange();

  private:
   inline void                                  Update(void);
   inline void                                  CleanUp(void);
   inline static G4String                       MakeKey(const G4String & source, const G4String & algorithm);

  // Hidden constructors and operators
                                                RadmonPrimaryGeneratorAction();
                                                RadmonPrimaryGeneratorAction(const RadmonPrimaryGeneratorAction & copy);
   RadmonPrimaryGeneratorAction &               operator=(const RadmonPrimaryGeneratorAction & copy);

  // Private data types
   typedef std::map<G4String, RadmonVGenerator *> GeneratorsMap;

  // Private attributes

   G4GeneralParticleSource * ParticleGunSource;

   G4ParticleGun                                particlesGun;
   RadmonVGeneratorLayout *                     generatorLayout;
   RadmonVGeneratorsFactory *                   generatorsFactory;
   G4bool                                       needUpdate;
   GeneratorsMap                                generatorsMap;
 };
#endif /* RADMONPRIMARYGENERATORACTION_HH */
