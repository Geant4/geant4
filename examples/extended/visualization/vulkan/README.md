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

3. Head-to-head OpenGL vs Vulkan Benchmark (generates `benchmark_opengl_vs_vulkan.json`):

```
% ./exampleVulkan --benchmark benchmark_opengl_vs_vulkan.json
```

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
17th October 2026

