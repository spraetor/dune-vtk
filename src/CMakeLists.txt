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

add_executable("quadratic" quadratic.cc)
target_link_dune_default_libraries("quadratic")
target_link_libraries("quadratic" dunevtk)

add_subdirectory(test)