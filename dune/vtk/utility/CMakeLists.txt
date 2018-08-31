dune_add_library("filesystem" OBJECT
  filesystem.cc)

#install headers
install(FILES
  enum.hh
  filesystem.hh
  string.hh
  uid.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtkwriter/utility)
