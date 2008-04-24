#ifdef G4VIS_USE

 #ifndef DetectorSliceVisManager_h
 #define DetectorSliceVisManager_h 1

 #include "G4VisManager.hh"

 class DetectorSliceVisManager: public G4VisManager {

 public:

   DetectorSliceVisManager();

 private:

   void RegisterGraphicsSystems ();

 };

 #endif

#endif
