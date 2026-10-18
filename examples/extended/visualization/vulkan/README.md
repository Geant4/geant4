\page Examplevulkan Example vulkan

This example demonstrates the Geant4 Vulkan visualization driver (`G4Vulkan`)
alongside the Geant4 OpenGL visualization driver (`G4OpenGL`), showing how both
drivers coexist independently and providing a built-in head-to-head benchmark
that outputs structured JSON metrics (`benchmark_opengl_vs_vulkan.json`).

The detector geometry consists of a world volume, an envelope volume, and a
concrete dump with a cylindrical cutout (`G4SubtractionSolid`).

 A- VISUALISATION DRIVERS
 ========================

The example registers all Vulkan visualization drivers alongside Geant4's
standard visualization drivers:

- `VulkanStoredQt` (`VKSQt`, `VKS`, `VK`)
- `VulkanImmediateQt` (`VKIQt`, `VKI`)
- `VulkanStoredX` (`VKSX`)
- `VulkanImmediateX` (`VKIX`)
- `OpenGLStoredQt` (`OGLSQt`, `OGLS`, `OGL`)
- `OpenGLImmediateQt` (`OGLIQt`, `OGLI`)
- `OpenGLStoredX` (`OGLSX`)
- `OpenGLImmediateX` (`OGLIX`)

 B- HOW TO RUN
 =============

1. Interactive mode (opens `VKSQt` via `init_vis.mac` / `vis.mac`):

```
% ./exampleVulkan
```

2. Batch macro mode with Vulkan commands (`vulkan.mac`) or OpenGL commands (`opengl.mac`):

```
% ./exampleVulkan vulkan.mac
% ./exampleVulkan opengl.mac
```

3. Head-to-head OpenGL vs Vulkan Benchmark:

```
% ./exampleVulkan --benchmark benchmark_opengl_vs_vulkan_linux.json
```

Benchmark files included in this directory:
- `benchmark_opengl_vs_vulkan.json` — Combined multi-platform benchmark report (macOS MoltenVK + Linux Native Vulkan)
- `benchmark_opengl_vs_vulkan_moltenvk.json` — macOS benchmark (MoltenVK / Metal vs Apple CGL OpenGL)
- `benchmark_opengl_vs_vulkan_linux.json` — Linux benchmark (Native Vulkan 1.4 Mesa RADV vs OpenGL 4.6 Mesa RadeonSI on AMD Radeon RX 570 Series)

| Metric | macOS (MoltenVK vs CGL OpenGL) | Linux (Mesa RADV Vulkan vs RadeonSI OpenGL) |
|---|---|---|
| Display list compile | 2364.6 us vs 1671.0 us (0.71x) | 1589.0 us vs 4769.1 us (**3.00x**) |
| Stored 256x256 | 2717 FPS vs 2067 FPS (**1.32x**) | 11278 FPS vs 11738 FPS (0.96x) |
| Stored 512x512 | 8910 FPS vs 2690 FPS (**3.31x**) | 22071 FPS vs 11245 FPS (**1.96x**) |
| Stored 1024x1024 | 9764 FPS vs 1278 FPS (**7.64x**) | 15360 FPS vs 6538 FPS (**2.35x**) |
| Immediate 512x512 | 1113 FPS vs 1446 FPS (0.77x) | 1093 FPS vs 867 FPS (**1.26x**) |
| 3D Picking (`SELECT`) | 0.688 us vs 128.6 us (**187.0x**) | 6.807 us vs 704.1 us (**103.4x**) |
| Vector EPS export | 3.46 ms vs 303.7 ms (**87.9x**) | 6.22 ms vs 279.9 ms (**45.0x**) |

 C- VULKAN COMMANDS (`/vis/vk/...`)
 ==================================

The Vulkan driver provides `/vis/vk/` messenger commands parallel to `/vis/ogl/`:

- `/vis/vk/export <filename>`
- `/vis/vk/printEPS`
- `/vis/vk/set/exportFormat (eps|ps|pdf|svg|png|jpg|bmp|tiff)`
- `/vis/vk/set/displayListLimit <limit>`
- `/vis/vk/set/transparency (true|false)`
- `/vis/vk/set/eventsDrawInterval <interval>`
- `/vis/vk/flushAt (endOfEvent|endOfRun|eachPrimitive|never)`

Alex Vergara Gil
18th October 2026


