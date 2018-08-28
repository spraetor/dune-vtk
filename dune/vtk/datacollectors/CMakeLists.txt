#install headers
install(FILES
  continuousdatacollector.hh
  datacollectorinterface.hh
  discontinuousdatacollector.hh
  quadraticdatacollector.hh
  spdatacollector.hh
  structureddatacollector.hh
  unstructureddatacollector.hh
  yaspdatacollector.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter/datacollectors)
