dune_add_library("vtktypes" OBJECT
  vtktypes.cc)

#install headers
install(FILES
  datacollectorinterface.hh
  datacollectorinterface.impl.hh
  defaultvtkfunction.hh
  filereader.hh
  filewriter.hh
  legacyvtkfunction.hh
  pvdwriter.hh
  pvdwriter.impl.hh
  vtkfunction.hh
  vtklocalfunction.hh
  vtklocalfunctioninterface.hh
  vtkreader.hh
  vtkreader.impl.hh
  vtktimeserieswriter.hh
  vtktimeserieswriter.impl.hh
  vtktypes.hh
  vtkwriter.hh
  vtkwriterinterface.hh
  vtkwriterinterface.impl.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter)

add_subdirectory(datacollectors)
add_subdirectory(gridcreators)
add_subdirectory(utility)
add_subdirectory(writers)
