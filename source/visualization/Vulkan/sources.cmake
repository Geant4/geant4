# - G4Vulkan module build definition

# Define the Geant4 Module.
geant4_add_module(G4Vulkan
  PUBLIC_HEADERS
    G4VulkanQt.hh
    G4VulkanXm.hh
  PRIVATE_HEADERS
    G4Vulkan.hh
    G4VulkanFontBaseStore.hh
    G4VulkanImmediateSceneHandler.hh
    G4VulkanImmediateViewer.hh
    G4VulkanSceneHandler.hh
    G4VulkanSceneHandler.icc
    G4VulkanStoredSceneHandler.hh
    G4VulkanStoredViewer.hh
    G4VulkanTransform3D.hh
    G4VulkanViewer.hh
    G4VulkanViewerMessenger.hh
    G4VisFeaturesOfVulkan.hh
    G4vk2ps.hh
  SOURCES
    G4VulkanFontBaseStore.cc
    G4VulkanImmediateSceneHandler.cc
    G4VulkanImmediateViewer.cc
    G4VulkanSceneHandler.cc
    G4VulkanStoredSceneHandler.cc
    G4VulkanStoredViewer.cc
    G4VulkanTransform3D.cc
    G4VulkanViewer.cc
    G4VulkanViewerMessenger.cc
    G4VisFeaturesOfVulkan.cc
    G4vk2ps.cc)

geant4_module_compile_definitions(G4Vulkan
  PRIVATE G4VIS_BUILD_VULKAN_DRIVER)

geant4_module_link_libraries(G4Vulkan
  PUBLIC
    G4modeling
    G4vis_management
  PRIVATE
    G4UImanager
    G4UIbasic
    G4UIcommon
    G4clhep
    G4globman
    G4Geometry
    G4graphics_reps
    G4hepnumerics
    G4intercoms
    G4materials
    G4tools
    Vulkan::Vulkan)

# Add specific options directly cleanly here
# Qt
if(GEANT4_USE_VULKAN_QT OR GEANT4_USE_QT)
  set(G4VULKAN_QT_MOC_HEADERS
    include/private/G4VulkanQtExportDialog.hh
    include/private/G4VulkanQtMovieDialog.hh
    include/private/G4VulkanQtViewer.hh)

  if(QT_VERSION_MAJOR GREATER 5)
    qt6_wrap_cpp(G4VULKAN_QT_MOC_FILES ${G4VULKAN_QT_MOC_HEADERS} TARGET G4Vulkan)
  else()
    qt5_wrap_cpp(G4VULKAN_QT_MOC_FILES ${G4VULKAN_QT_MOC_HEADERS} TARGET G4Vulkan)
  endif()

  geant4_module_sources(G4Vulkan
    PUBLIC_HEADERS
      G4VulkanImmediateQt.hh
      G4VulkanStoredQt.hh
    PRIVATE_HEADERS
      G4VulkanImmediateQtViewer.hh
      G4VulkanQtExportDialog.hh
      G4VulkanQtMovieDialog.hh
      G4VulkanQtViewer.hh
      G4VulkanStoredQtSceneHandler.hh
      G4VulkanStoredQtViewer.hh
      ${G4VULKAN_QT_MOC_FILES}
    SOURCES
      G4VulkanImmediateQt.cc
      G4VulkanImmediateQtViewer.cc
      G4VulkanQt.cc
      G4VulkanQtExportDialog.cc
      G4VulkanQtMovieDialog.cc
      G4VulkanQtViewer.cc
      G4VulkanStoredQt.cc
      G4VulkanStoredQtSceneHandler.cc
      G4VulkanStoredQtViewer.cc)

  geant4_module_compile_definitions(G4Vulkan
    PUBLIC G4VIS_USE_VULKANQT
    PRIVATE G4VIS_BUILD_VULKANQT_DRIVER)

  # Workaround cleanly at finest scope - Very unclear why changing to target links caused this to fail.
  geant4_module_include_directories(G4Vulkan
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

  geant4_module_link_libraries(G4Vulkan
    PRIVATE
      Qt${QT_VERSION_MAJOR}::Gui
      Qt${QT_VERSION_MAJOR}::Widgets)
endif()

# X11
if(GEANT4_USE_VULKAN_X11)
  geant4_module_sources(G4Vulkan
    PUBLIC_HEADERS
      G4VulkanImmediateX.hh
      G4VulkanStoredX.hh
    PRIVATE_HEADERS
      G4VulkanImmediateXViewer.hh
      G4VulkanStoredXViewer.hh
      G4VulkanXViewer.hh
    SOURCES
      G4VulkanImmediateX.cc
      G4VulkanImmediateXViewer.cc
      G4VulkanStoredX.cc
      G4VulkanStoredXViewer.cc
      G4VulkanXViewer.cc)

  geant4_module_compile_definitions(G4Vulkan
    PUBLIC G4VIS_USE_VULKANX
    PRIVATE G4VIS_BUILD_VULKANX_DRIVER)

  geant4_module_link_libraries(G4Vulkan
    PRIVATE
      X11::SM
      X11::ICE
      X11::X11
      X11::Xext
      X11::Xmu)
endif()

# Win32
if(GEANT4_USE_VULKAN_WIN32)
  geant4_module_sources(G4Vulkan
    PUBLIC_HEADERS
      G4VulkanImmediateWin32.hh
      G4VulkanStoredWin32.hh
    PRIVATE_HEADERS
      G4VulkanImmediateWin32Viewer.hh
      G4VulkanStoredWin32Viewer.hh
      G4VulkanWin32Viewer.hh
    SOURCES
      G4VulkanImmediateWin32.cc
      G4VulkanImmediateWin32Viewer.cc
      G4VulkanStoredWin32.cc
      G4VulkanStoredWin32Viewer.cc
      G4VulkanWin32Viewer.cc)

  geant4_module_compile_definitions(G4Vulkan
    PUBLIC G4VIS_USE_VULKANWIN32
    PRIVATE G4VIS_BUILD_VULKANWIN32_DRIVER)
endif()

# Xm
if(GEANT4_USE_VULKAN_XM)
  geant4_module_sources(G4Vulkan
    PUBLIC_HEADERS
      G4VulkanImmediateXm.hh
      G4VulkanStoredXm.hh
    PRIVATE_HEADERS
      G4VulkanImmediateXmViewer.hh
      G4VulkanStoredXmViewer.hh
      G4VulkanXmBox.hh
      G4VulkanXmFourArrowButtons.hh
      G4VulkanXmFramedBox.hh
      G4VulkanXmPushButton.hh
      G4VulkanXmRadioButton.hh
      G4VulkanXmResources.hh
      G4VulkanXmSeparator.hh
      G4VulkanXmSliderBar.hh
      G4VulkanXmTextField.hh
      G4VulkanXmTopLevelShell.hh
      G4VulkanXmVWidgetComponent.hh
      G4VulkanXmVWidgetContainer.hh
      G4VulkanXmVWidgetObject.hh
      G4VulkanXmVWidgetShell.hh
      G4VulkanXmViewer.hh
      G4VulkanXmViewerMessenger.hh
      G4VulkanXViewer.hh
    SOURCES
      G4VulkanImmediateXm.cc
      G4VulkanImmediateXmViewer.cc
      G4VulkanStoredXm.cc
      G4VulkanStoredXmViewer.cc
      G4VulkanXm.cc
      G4VulkanXmBox.cc
      G4VulkanXmFourArrowButtons.cc
      G4VulkanXmFramedBox.cc
      G4VulkanXmMiscCallbacks.cc
      G4VulkanXmPanningCallbacks.cc
      G4VulkanXmPushButton.cc
      G4VulkanXmRadioButton.cc
      G4VulkanXmRotationCallbacks.cc
      G4VulkanXmSeparator.cc
      G4VulkanXmSettingCallbacks.cc
      G4VulkanXmSliderBar.cc
      G4VulkanXmStyleCallbacks.cc
      G4VulkanXmTextField.cc
      G4VulkanXmTopLevelShell.cc
      G4VulkanXmVWidgetComponent.cc
      G4VulkanXmVWidgetContainer.cc
      G4VulkanXmVWidgetObject.cc
      G4VulkanXmVWidgetShell.cc
      G4VulkanXmViewer.cc
      G4VulkanXmViewerMessenger.cc
      G4VulkanXmWindowHandlingCallbacks.cc
      G4VulkanXViewer.cc)

  geant4_module_compile_definitions(G4Vulkan
    PUBLIC G4VIS_USE_VULKANXM
    PRIVATE G4VIS_BUILD_VULKANXM_DRIVER)

  geant4_module_link_libraries(G4Vulkan
    PRIVATE
      X11::X11
      X11::Xt
      Motif::Xm)
endif()

