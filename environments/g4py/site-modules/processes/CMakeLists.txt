# - add libs components

if (${PYTHON_VERSION_MAJOR} MATCHES "2")
  add_subdirectory(emcalculator)
elseif (${PYTHON_VERSION_MAJOR} MATCHES "3")
  add_subdirectory(emcalculator/python3)
endif()
