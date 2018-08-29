dune_add_test(SOURCES reader_writer_test.cc
              LINK_LIBRARIES dunevtk
              CMAKE_GUARD)

dune_add_test(SOURCES mixed_element_test.cc
              LINK_LIBRARIES dunevtk
              CMAKE_GUARD HAVE_UG)