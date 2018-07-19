// $Id: Doxymodules_biasing.h 104750 2017-06-15 08:47:26Z gcosmo $
// The example class categories definitions for Doxygen

/// \file Doxymodules_biasing.h
/// \brief The page that defines the extended/biasing examples modules 


/** @defgroup extended_biasing biasing
 *  Extended examples biasing classes
 *  @{
 */

/** @defgroup extended_biasing_B01 B01
 *  Biasing example B01
 *  @ingroup extended_biasing
 *  @{
 */

  class B01ActionInitialization {};
  class B01DetectorConstruction {};
  class B01PrimaryGeneratorAction {};
  class B01RunAction {};
  class B01Run {};


/** @} */

/** @defgroup extended_biasing_B02 B02
 *  Biasing example B02
 *  @ingroup extended_biasing
 *  @{
 */

  class B02ActionInitialization {};
  class B02DetectorConstruction {};
  class B02ImportanceDetectorConstruction {};
  class B02PrimaryGeneratorAction {};
  class B02PVolumeStore {};
  class B02RunAction {};
  class B02Run {};

/** @} */

/** @defgroup extended_biasing_B03 B03
 *  Biasing example B03
 *  @ingroup extended_biasing
 *  @{
 */

  class B03ActionInitialization {};
  class B03DetectorConstruction {};
  class B03ImportanceDetectorConstruction {};
  class B03PhysicsList {};
  class B03PrimaryGeneratorAction {};
  class B03PVolumeStore {};
  class B03RunAction {};
  class B03Run {};

/** @} */

/** @defgroup extended_biasing_GB01 GB01
 *  Biasing example GB01
 *  @ingroup extended_biasing
 *  @{
 */

  class GB01ActionInitialization {};
  class GB01BOptrChangeCrossSection {};
  class GB01BOptrMultiParticleChangeCrossSection {};
  class GB01DetectorConstruction {};
  class GB01PrimaryGeneratorAction {};

/** @} */

/** @defgroup extended_biasing_GB02 GB02
 *  Biasing example GB02
 *  @ingroup extended_biasing
 *  @{
 */

  class GB02ActionInitialization {};
  class GB02BOptrMultiParticleForceCollision {};
  class GB02DetectorConstruction {};
  class GB02PrimaryGeneratorAction {};

/** @} */

/** @defgroup extended_biasing_GB03 GB03
 *  Biasing example GB03
 *  @ingroup extended_biasing
 *  @{
 */

  class GB03ActionInitialization {};
  class GB03BOptnSplitOrKillOnBoundary{};
  class GB03BOptrGeometryBasedBiasing{};
  class GB03DetectorConstruction {};
  class GB03DetectorMessenger{};
  class GB03PrimaryGeneratorAction {};

/** @} */

/** @defgroup extended_biasing_GB04 GB04
 *  Biasing example GB04
 *  @ingroup extended_biasing
 *  @{
 */

  class GB04ActionInitialization {};
  class GB04BOptnBremSplitting {};
  class GB04BOptrBremSplitting {};
  class GB04DetectorConstruction {};
  class GB04PrimaryGeneratorAction {};

/** @} */

/** @defgroup extended_biasing_GB05 GB05
 *  Biasing example GB05
 *  @ingroup extended_biasing
 *  @{
 */

  class GB05ActionInitialization {};
  class GB05BOptnSplitAndKillByCrossSection {};
  class GB05BOptrSplitAndKillByCrossSection {};
  class GB05DetectorConstruction {};
  class GB05PrimaryGeneratorAction {};
  class GB05SD {};

/** @} */

/** @defgroup extended_biasing_GB06 GB06
 *  Biasing example GB06
 *  @ingroup extended_biasing
 *  @{
 */

  class GB06ActionInitialization {};
  class GB06BOptnSplitAndKillByImportance {};
  class GB06BOptrSplitAndKillByImportance {};
  class GB06DetectorConstruction {};
  class GB06ParallelWorldForSlices {};
  class GB06PrimaryGeneratorAction {};
  class GB06SD {};

/** @} */

/** @defgroup extended_biasing_ReverseMC01 ReverseMC01
 *  Biasing example ReverseMC01
 *  @ingroup extended_biasing
 *  @{
 */

  class G4AdjointPhysicsList {};
  class G4AdjointPhysicsMessenger {};
  class RMC01AdjointEventAction {};
  class RMC01AnalysisManager {};
  class RMC01AnalysisManagerMessenger {};
  class RMC01DetectorConstruction {};
  class RMC01DetectorMessenger {};
  class RMC01DoubleWithWeightHit {};
  class RMC01EventAction {};
  class RMC01PrimaryGeneratorAction {};
  class RMC01RunAction {};
  class RMC01SD {};

/** @} */

/** @} */
