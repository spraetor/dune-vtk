add_executable("vtkwriter" vtkwriter.cc)
target_link_dune_default_libraries("vtkwriter")
target_link_libraries("vtkwriter" dunevtk)

add_executable("vtkreader" vtkreader.cc)
target_link_dune_default_libraries("vtkreader")
target_link_libraries("vtkreader" dunevtk)

add_executable("benchmark" benchmark.cc)
target_link_dune_default_libraries("benchmark")
target_link_libraries("benchmark" dunevtk)

add_executable("datacollector" datacollector.cc)
target_link_dune_default_libraries("datacollector")
target_link_libraries("datacollector" dunevtk)

add_executable("structuredgridwriter" structuredgridwriter.cc)
target_link_dune_default_libraries("structuredgridwriter")
target_link_libraries("structuredgridwriter" dunevtk)

add_executable("legacyvtkwriter" legacyvtkwriter.cc)
target_link_dune_default_libraries("legacyvtkwriter")
target_link_libraries("legacyvtkwriter" dunevtk)

if (dune-polygongrid_FOUND)
  add_executable("polygongrid" polygongrid.cc)
  target_link_dune_default_libraries("polygongrid")
  target_link_libraries("polygongrid" dunevtk dunepolygongrid)
endif ()

add_subdirectory(test)