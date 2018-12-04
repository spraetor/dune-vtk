#install headers
install(FILES
  common.hh
  continuousgridcreator.hh
  derivedgridcreator.hh
  discontinuousgridcreator.hh
  parallelgridcreator.hh
  serialgridcreator.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter/gridcreators)
