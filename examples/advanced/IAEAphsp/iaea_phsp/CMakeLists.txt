add_library(iaea_phsp STATIC iaea_header.cpp iaea_phsp.cpp iaea_record.cpp utilities.cpp)
target_include_directories(iaea_phsp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
