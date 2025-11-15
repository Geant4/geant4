# Source Modules Reference

## Overview

The Geant4 source code is organized into specialized modules located in the `source/` directory. Each module handles a specific aspect of particle simulation.

## Module Categories

### Physics & Particles

#### particles/
Particle definitions and properties for all supported particle types.

- Defines particle characteristics (mass, charge, lifetime, etc.)
- Includes fundamental particles, nuclei, and ions
- Provides particle property databases

#### processes/
Implementation of physical processes governing particle interactions.

- Electromagnetic processes (ionization, bremsstrahlung, Compton scattering, etc.)
- Hadronic processes (nuclear interactions, fission, etc.)
- Optical processes (reflection, refraction, scintillation, etc.)
- Transportation and decay processes

#### physics_lists/
Pre-configured collections of physics processes for specific use cases.

- Standard physics lists for different energy ranges
- Modular physics constructors
- Reference physics lists validated against experimental data

#### materials/
Material definitions and properties.

- Element and material databases
- Material property tables
- Cross-section calculations

### Geometry & Visualization

#### geometry/
Geometric modeling and navigation capabilities.

- Solid definitions (box, sphere, cylinder, etc.)
- Boolean operations for complex shapes
- Navigation algorithms for particle tracking through geometry
- Volume hierarchy management

#### graphics_reps/
Graphical representations for visualization.

- Primitive graphical objects
- Visualization attributes
- Scene graph representations

#### visualization/
Visualization drivers and tools.

- Multiple visualization backends (OpenGL, Qt, VRML, etc.)
- Interactive and batch visualization modes
- Trajectory and hit visualization

### Event Processing

#### event/
Event generation and management.

- Primary particle generation
- Event structure and storage
- Stacking mechanism for secondary particles

#### run/
Run management and control.

- Run initialization and termination
- Worker/master thread management (for multithreading)
- Random number engine management

#### track/
Track representation and management.

- Track information (position, momentum, energy, etc.)
- Step information
- Track status handling

#### tracking/
Particle tracking algorithms and management.

- Track propagation through detector
- Step management
- User action hooks

### Detector Response

#### digits_hits/
Detector response simulation and sensitive detector handling.

- Sensitive detector base classes
- Hit collections
- Digitization framework

#### readout/
Readout geometry and detector response.

- Parallel readout geometries
- Detector segmentation
- Signal processing

### User Interface & Control

#### interfaces/
User interface implementations.

- Command-line interface
- Graphical user interfaces (Qt, etc.)
- Batch mode support

#### intercoms/
Command system and messenger infrastructure.

- Command directory structure
- Parameter validation
- User action messaging

### Analysis & Output

#### analysis/
Analysis tools and histogramming.

- Histogram creation and management
- N-tuple support
- Multiple output formats (ROOT, XML, CSV)

#### persistency/
Data persistency mechanisms.

- Event data output
- Geometry persistence
- Multiple formats support

### Utilities & Advanced Features

#### global/
Global utilities and definitions.

- System of units
- Physical constants
- Random number generators
- Exception handling
- Memory management utilities

#### externals/
External dependencies and third-party libraries.

- Integrated external code
- Compatibility layers

#### parameterisations/
Fast simulation parameterizations.

- Shower parameterizations
- Simplified physics models for speed
- Geometry parameterizations

#### error_propagation/
Error propagation utilities for systematic studies.

- Track error propagation
- Uncertainty quantification

#### g3tog4/
Geant3 to Geant4 geometry conversion tools.

- Legacy geometry conversion
- Geant3 compatibility layer

## Module Dependencies

The modules are interconnected with clear dependency hierarchies:

```
global
  ├── geometry
  ├── materials
  ├── particles
  │   └── processes
  │       └── physics_lists
  ├── track
  │   └── tracking
  ├── event
  │   └── run
  ├── digits_hits
  │   └── readout
  └── visualization
```

## Build Configuration

Each module has its own `CMakeLists.txt` for granular build control.

Main source build configuration: `source/CMakeLists.txt`

## Module Usage

When building applications with Geant4:

1. **Required modules**: Always linked (global, geometry, materials, particles, track, event, run)
2. **Optional modules**: Can be enabled/disabled based on needs (visualization, analysis, etc.)
3. **Physics lists**: Choose appropriate physics list for your application domain

## Further Information

For detailed documentation on using these modules:
- [Official Documentation](http://cern.ch/geant4/support/user_documentation)
- [Architecture Overview](../architecture)
- [Getting Started Guide](../getting-started)

## Related Files

| Directory | Purpose |
|-----------|---------|
| `source/` | All source modules |
| `source/CMakeLists.txt` | Source build configuration |
| `examples/` | Example usage of modules |
