find_package(ZLIB)
set(HAVE_ZLIB ${ZLIB_FOUND})
if (${HAVE_ZLIB})
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_ZLIB=1"
                              LIBRARIES "${ZLIB_LIBRARIES}"
                              INCLUDE_DIRS "${ZLIB_INCLUDE_DIRS}")
endif (${HAVE_ZLIB})
