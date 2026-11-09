// ********************************************************************
// * GPU Acceleration Prototype: Structure of Arrays (SoA) for G4Track  *
// ********************************************************************
#ifndef GPU_TRACK_SOA_HH
#define GPU_TRACK_SOA_HH

#include "globals.hh"
#include <vector>

/**
 * GpuTrackSoA replaces the array of G4Track objects with a Structure of Arrays.
 * This is essential for GPU SIMD processing to avoid memory divergence.
 */
class GpuTrackSoA {
  public:
    GpuTrackSoA(size_t maxTracks) : m_maxTracks(maxTracks) {
      posX.resize(maxTracks, 0.0);
      posY.resize(maxTracks, 0.0);
      posZ.resize(maxTracks, 0.0);
      momX.resize(maxTracks, 0.0);
      momY.resize(maxTracks, 0.0);
      momZ.resize(maxTracks, 0.0);
      energy.resize(maxTracks, 0.0);
      trackID.resize(maxTracks, 0);
      active.resize(maxTracks, false);
    }

    // Data arrays for GPU-friendly access
    std::vector<G4double> posX, posY, posZ;
    std::vector<G4double> momX, momY, momZ;
    std::vector<G4double> energy;
    std::vector<G4int> trackID;
    std::vector<G4bool> active;

  private:
    size_t m_maxTracks;
};

#endif
