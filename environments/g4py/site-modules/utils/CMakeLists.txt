# - add libs components

if (${PYTHON_VERSION_MAJOR} MATCHES "2")
  add_subdirectory(MCScore)
elseif (${PYTHON_VERSION_MAJOR} MATCHES "3")
  add_subdirectory(MCScore/python3)
endif()
