#ifdef G4VIS_USE

 #ifndef RandomCaloVisManager_h
 #define RandomCaloVisManager_h 1

 #include "G4VisManager.hh"

 class RandomCaloVisManager: public G4VisManager {

 public:

   RandomCaloVisManager();

 private:

   void RegisterGraphicsSystems ();

 };

 #endif

#endif
