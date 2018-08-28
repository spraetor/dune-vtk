dune_add_library("vtktypes" OBJECT
  vtktypes.cc)

#install headers
install(FILES
  defaultvtkfunction.hh
  filereader.hh
  filewriter.hh
  gridcreator.hh
  legacyvtkfunction.hh
  vtkfunction.hh
  vtklocalfunction.hh
  vtklocalfunctioninterface.hh
  vtkreader.hh
  vtkreader.impl.hh
  vtktypes.hh
  vtkwriter.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter)

add_subdirectory(datacollectors)
add_subdirectory(utility)
add_subdirectory(writers)
