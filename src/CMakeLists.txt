add_executable("vtkwriter" vtkwriter.cc)
target_link_dune_default_libraries("vtkwriter")
target_link_libraries("vtkwriter" dunevtk)

add_executable("vtkreader" vtkreader.cc)
target_link_dune_default_libraries("vtkreader")
target_link_libraries("vtkreader" dunevtk)