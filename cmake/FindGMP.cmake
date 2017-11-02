find_path(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h)
find_library(GMP_LIBRARY NAMES gmp libgmp)
find_library(GMPXX_LIBRARY NAMES gmpxx libgmpxx)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARY GMPXX_LIBRARY GMP_INCLUDE_DIR)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

set(GMP_LIBRARIES ${GMP_LIBRARY} ${GMPXX_LIBRARY})
set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
