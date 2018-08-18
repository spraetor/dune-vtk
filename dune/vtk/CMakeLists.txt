dune_add_library("vtktypes" OBJECT
  vtktypes.cc)

#install headers
install(FILES
  filereader.hh
  filewriter.hh
  vtkfunction.hh
  vtkreader.hh
  vtkreader.impl.hh
  vtktypes.hh
  vtkwriter.hh
  vtkwriter.impl.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter)

add_subdirectory(utility)
