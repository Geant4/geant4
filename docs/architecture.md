# Geant4 Architecture

## Overview

Geant4 is designed as a modular, object-oriented toolkit for particle physics simulations. The architecture is organized into distinct modules, each handling specific aspects of the simulation process.

## Core Components

The Geant4 source code is organized into the following main categories (located in `source/`):

### Physics & Particle Management

- **`particles/`**: Particle definitions and properties
- **`processes/`**: Physical processes (electromagnetic, hadronic, optical, etc.)
- **`physics_lists/`**: Pre-configured physics models for different use cases
- **`materials/`**: Material definitions and properties

### Geometry & Tracking

- **`geometry/`**: Geometric modeling and navigation
- **`graphics_reps/`**: Graphical representations
- **`track/`**: Track management and representation
- **`tracking/`**: Particle tracking algorithms

### Event & Run Management

- **`event/`**: Event generation and management
- **`run/`**: Run management and control
- **`digits_hits/`**: Detector response simulation
- **`readout/`**: Sensitive detector readout

### Visualization & Interfaces

- **`visualization/`**: Visualization drivers and tools
- **`interfaces/`**: User interface implementations
- **`intercoms/`**: Command system and messengers

### Analysis & Data

- **`analysis/`**: Analysis tools and histogramming
- **`persistency/`**: Data persistency mechanisms

### Utilities & Extensions

- **`global/`**: Global utilities and definitions
- **`externals/`**: External dependencies
- **`parameterisations/`**: Fast simulation parameterizations
- **`error_propagation/`**: Error propagation utilities
- **`g3tog4/`**: Geant3 to Geant4 conversion tools

## Design Principles

### Modularity

Geant4's modular design allows users to:
- Include only the components they need
- Extend functionality through inheritance
- Customize physics processes for specific applications

### Object-Oriented Design

The toolkit leverages C++ object-oriented features:
- **Inheritance**: Base classes define interfaces, derived classes provide implementations
- **Polymorphism**: Common interfaces for different physics models
- **Encapsulation**: Clear separation of concerns between modules

### Flexibility

Users can customize:
- Physics models and processes
- Detector geometry
- Event generation
- Output and analysis

## Build System

Geant4 uses **CMake** as its build system (version 3.16-3.27).

Key build files:
- `CMakeLists.txt`: Top-level build configuration (CMakeLists.txt:1-50)
- `cmake/Modules/`: Custom CMake modules
- `source/CMakeLists.txt`: Source-level build configuration

### Build Features

The build system:
- Enforces out-of-source builds (CMakeLists.txt:6-14)
- Manages inter-module dependencies
- Supports multiple platforms (Linux, macOS, Windows)
- Allows optional component selection

## Application Areas

Geant4 is used in:

1. **High Energy Physics**: Particle detector simulation, accelerator physics
2. **Nuclear Physics**: Nuclear reactions and radiation transport
3. **Medical Physics**: Radiotherapy, medical imaging, radiation protection
4. **Space Science**: Radiation effects in spacecraft and instruments
5. **Industrial Applications**: Non-destructive testing, radiation shielding

## Reference Papers

For detailed information on Geant4's architecture and physics models, see:

1. Nuclear Instruments and Methods in Physics Research A 835 (2016) 186-225
2. IEEE Transactions on Nuclear Science 53 No. 1 (2006) 270-278
3. Nuclear Instruments and Methods in Physics Research A 506 (2003) 250-303

## Further Reading

- [Source Modules Reference](./reference/source-modules): Detailed module descriptions
- [Build System Reference](./reference/build-system): In-depth build configuration
- [Official Documentation](http://cern.ch/geant4/support/user_documentation): Complete user guides
