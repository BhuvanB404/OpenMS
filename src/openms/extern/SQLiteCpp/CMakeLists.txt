# Main CMake file for compiling the library itself, examples and tests.
#
# Copyright (c) 2012-2022 Sebastien Rombauts (sebastien.rombauts@gmail.com)
#
# Distributed under the MIT License (MIT) (See accompanying file LICENSE.txt
# or copy at http://opensource.org/licenses/MIT)
cmake_minimum_required(VERSION 3.10) # for "CMAKE_CXX_STANDARD" version
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake") # custom CMake modules like FindSQLiteCpp
project(SQLiteCpp VERSION 3.2.0)

option(SQLITECPP_BUILD_TESTS "Build and run tests." OFF)

# Define useful variables to handle OS differences:
if (WIN32)
    set(DEV_NULL                "NUL")
else (WIN32) # UNIX
    set(DEV_NULL                "/dev/null")
endif (WIN32)

## Build the C++ Wrapper ##

# adding a new file require explicitly modifying the CMakeLists.txt
# so that CMake knows that it should rebuild the project (it is best practice)

# list of sources files of the library
set(SQLITECPP_SRC
 ${PROJECT_SOURCE_DIR}/src/Backup.cpp
 ${PROJECT_SOURCE_DIR}/src/Column.cpp
 ${PROJECT_SOURCE_DIR}/src/Database.cpp
 ${PROJECT_SOURCE_DIR}/src/Exception.cpp
 ${PROJECT_SOURCE_DIR}/src/Savepoint.cpp
 ${PROJECT_SOURCE_DIR}/src/Statement.cpp
 ${PROJECT_SOURCE_DIR}/src/Transaction.cpp
)
source_group(src FILES ${SQLITECPP_SRC})

# list of header files of the library
set(SQLITECPP_INC
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/SQLiteCpp.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Assertion.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Backup.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Column.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Database.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Exception.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Savepoint.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Statement.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/Transaction.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/VariadicBind.h
 ${PROJECT_SOURCE_DIR}/include/SQLiteCpp/ExecuteMany.h
)
source_group(include FILES ${SQLITECPP_INC})

# list of test files of the library
set(SQLITECPP_TESTS
 tests/Column_test.cpp
 tests/Database_test.cpp
 tests/Savepoint_test.cpp
 tests/Statement_test.cpp
 tests/Backup_test.cpp
 tests/Transaction_test.cpp
 tests/VariadicBind_test.cpp
 tests/Exception_test.cpp
 tests/ExecuteMany_test.cpp
)
source_group(tests FILES ${SQLITECPP_TESTS})

# list of example files of the library
set(SQLITECPP_EXAMPLES
 examples/example1/main.cpp
)
source_group(example1 FILES ${SQLITECPP_EXAMPLES})

# list of doc files of the library
set(SQLITECPP_DOC
 README.md
 LICENSE.txt
 CHANGELOG.md
 TODO.txt
)
source_group(doc FILES ${SQLITECPP_DOC})

option(SQLITECPP_INCLUDE_SCRIPT "Include config & script files." OFF)
if (SQLITECPP_INCLUDE_SCRIPT)
    # list of config & script files of the library
    set(SQLITECPP_SCRIPT
     .editorconfig
     .gitbugtraq
     .github/workflows/build.yml
     .github/workflows/subdir_example.yml
     .gitignore
     .gitmodules
     .travis.yml
     appveyor.yml
     build.bat
     build.sh
     cpplint.py
     Doxyfile
     cmake/FindSQLite3.cmake
     cmake/SQLiteCppConfig.cmake.in
    )
    source_group(scripts FILES ${SQLITECPP_SCRIPT})
endif()

# add sources of the wrapper as a "SQLiteCpp" static library
add_library(SQLiteCpp STATIC ${SQLITECPP_SRC} ${SQLITECPP_INC} ${SQLITECPP_DOC} ${SQLITECPP_SCRIPT})
target_compile_features(SQLiteCpp PUBLIC cxx_std_11)

# Options relative to SQLite and SQLiteC++ functions

option(SQLITE_ENABLE_COLUMN_METADATA "Enable Column::getColumnOriginName(). Require support from sqlite3 library." ON)
if (SQLITE_ENABLE_COLUMN_METADATA)
    # Enable the use of SQLite column metadata and Column::getColumnOriginName() method,
    # Require that the sqlite3 library is also compiled with this flag (default under Debian/Ubuntu, but not on Mac OS X).
    target_compile_definitions(SQLiteCpp PUBLIC SQLITE_ENABLE_COLUMN_METADATA)
endif (SQLITE_ENABLE_COLUMN_METADATA)

option(SQLITE_ENABLE_ASSERT_HANDLER "Enable the user definition of a assertion_failed() handler." OFF)
if (SQLITE_ENABLE_ASSERT_HANDLER)
    # Enable the user definition of a assertion_failed() handler (default to false, easier to handler for beginners).
    target_compile_definitions(SQLiteCpp PUBLIC SQLITECPP_ENABLE_ASSERT_HANDLER)
endif (SQLITE_ENABLE_ASSERT_HANDLER)

option(SQLITE_HAS_CODEC "Enable database encryption API. Not available in the public release of SQLite." OFF)
if (SQLITE_HAS_CODEC)
    # Enable database encryption API. Requires implementations of sqlite3_key & sqlite3_key_v2.
    # Eg. SQLCipher (libsqlcipher-dev) is an SQLite extension that provides 256 bit AES encryption of database files. 
    target_compile_definitions(SQLiteCpp PUBLIC SQLITE_HAS_CODEC)
endif (SQLITE_HAS_CODEC)

option(SQLITE_USE_LEGACY_STRUCT "Fallback to forward declaration of legacy struct sqlite3_value (pre SQLite 3.19)" OFF)
if (SQLITE_USE_LEGACY_STRUCT)
    # Force forward declaration of legacy struct sqlite3_value (pre SQLite 3.19)
    target_compile_definitions(SQLiteCpp PUBLIC SQLITE_USE_LEGACY_STRUCT)
endif (SQLITE_USE_LEGACY_STRUCT)

option(SQLITE_OMIT_LOAD_EXTENSION "Enable omit load extension" OFF)
if (SQLITE_OMIT_LOAD_EXTENSION)
    # Enable the user definition of load_extension().
    target_compile_definitions(SQLiteCpp PUBLIC SQLITE_OMIT_LOAD_EXTENSION)
endif (SQLITE_OMIT_LOAD_EXTENSION)

if (UNIX AND (CMAKE_COMPILER_IS_GNUCXX OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    set_target_properties(SQLiteCpp PROPERTIES COMPILE_FLAGS "-fPIC")
endif (UNIX AND (CMAKE_COMPILER_IS_GNUCXX OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))

option(SQLITECPP_USE_ASAN "Use Address Sanitizer." OFF)
if (SQLITECPP_USE_ASAN)
    if ((CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 6) OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
        message (STATUS "Using Address Sanitizer")
        set_target_properties(SQLiteCpp PROPERTIES COMPILE_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
        if (CMAKE_COMPILER_IS_GNUCXX)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
        endif ()
    endif ()
endif (SQLITECPP_USE_ASAN)

if (SQLITECPP_USE_GCOV)
    # Prevent the compiler from removing the unused inline functions so that they get tracked as "non-covered"
    set_target_properties(SQLiteCpp PROPERTIES COMPILE_FLAGS "-fkeep-inline-functions -fkeep-static-functions")
endif ()

## Build provided copy of SQLite3 C library ##

option(SQLITECPP_INTERNAL_SQLITE "Add the internal SQLite3 source to the project." ON)
if (SQLITECPP_INTERNAL_SQLITE)
    message(STATUS "Compile sqlite3 from source in subdirectory")
    option(SQLITE_ENABLE_JSON1 "Enable JSON1 extension when building internal sqlite3 library." ON)
    # build the SQLite3 C library (for ease of use/compatibility) versus Linux sqlite3-dev package
    add_subdirectory(sqlite3)
    target_link_libraries(SQLiteCpp PUBLIC SQLite::SQLite3)
else (SQLITECPP_INTERNAL_SQLITE)
    # When using the SQLite codec, we need to link against the sqlcipher lib & include <sqlcipher/sqlite3.h>
    # So this gets the lib & header, and links/includes everything
    if(SQLITE_HAS_CODEC)
        # Make PkgConfig optional since Windows doesn't usually have it installed.
        find_package(PkgConfig QUIET)
        if(PKG_CONFIG_FOUND)
            # IMPORTED_TARGET was added in 3.6.3
            if(CMAKE_VERSION VERSION_LESS 3.6.3)
                pkg_check_modules(sqlcipher REQUIRED sqlcipher)
                # Only used in Database.cpp so PRIVATE to hide from end-user
                # Since we can't use IMPORTED_TARGET on this older Cmake version, manually link libs & includes
                target_link_libraries(SQLiteCpp PRIVATE ${sqlcipher_LIBRARIES})
                target_include_directories(SQLiteCpp PRIVATE ${sqlcipher_INCLUDE_DIRS})
            else()
                pkg_check_modules(sqlcipher REQUIRED IMPORTED_TARGET sqlcipher)
                # Only used in Database.cpp so PRIVATE to hide from end-user
                target_link_libraries(SQLiteCpp PRIVATE PkgConfig::sqlcipher)
            endif()
        else()
            # Since we aren't using pkgconf here, find it manually
            find_library(sqlcipher_LIBRARY "sqlcipher")
            find_path(sqlcipher_INCLUDE_DIR "sqlcipher/sqlite3.h"
                PATH_SUFFIXES
                    "include"
                    "includes"
            )
            # Hides it from the GUI
            mark_as_advanced(sqlcipher_LIBRARY sqlcipher_INCLUDE_DIR)
            if(NOT sqlcipher_INCLUDE_DIR)
                message(FATAL_ERROR "${PROJECT_NAME} requires the \"<sqlcipher/sqlite3.h>\" header to use the codec functionality but it wasn't found.")
            elseif(NOT sqlcipher_LIBRARY)
                message(FATAL_ERROR "${PROJECT_NAME} requires the sqlcipher library to use the codec functionality but it wasn't found.")
            endif()
            # Only used in Database.cpp so PRIVATE to hide from end-user
            target_include_directories(SQLiteCpp PRIVATE "${sqlcipher_INCLUDE_DIR}/sqlcipher")
            target_link_libraries(SQLiteCpp PRIVATE ${sqlcipher_LIBRARY})
        endif()
    else()
        find_package (SQLite3 REQUIRED)
        message(STATUS "Link to sqlite3 system library")
        target_link_libraries(SQLiteCpp PUBLIC SQLite::SQLite3)
        if(SQLite3_VERSION VERSION_LESS "3.19")
            set_target_properties(SQLiteCpp PROPERTIES COMPILE_FLAGS "-DSQLITECPP_HAS_MEM_STRUCT")
        endif()
    endif()
endif (SQLITECPP_INTERNAL_SQLITE)

# Link target with pthread and dl for Unix
if (UNIX)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(SQLiteCpp PUBLIC Threads::Threads ${CMAKE_DL_LIBS})
endif (UNIX)

# Set includes for target and transitive downstream targets

target_include_directories(SQLiteCpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/>)

# Optional additional targets:

# Uses ancient python find module!! DO NOT USE!
#option(SQLITECPP_RUN_CPPLINT "Run cpplint.py tool for Google C++ StyleGuide." OFF)
#if (SQLITECPP_RUN_CPPLINT)
#    find_package(PythonInterp)
#    if (PYTHONINTERP_FOUND)
#        # add a cpplint target to the "all" target
#        add_custom_target(SQLiteCpp_cpplint
#         ALL
#         COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cpplint.py ${CPPLINT_ARG_OUTPUT} ${CPPLINT_ARG_VERBOSE} ${CPPLINT_ARG_LINELENGTH} ${SQLITECPP_SRC} ${SQLITECPP_INC}
#        )
#    endif (PYTHONINTERP_FOUND)
#else (SQLITECPP_RUN_CPPLINT)
#    message(STATUS "SQLITECPP_RUN_CPPLINT OFF")
#endif (SQLITECPP_RUN_CPPLINT)

option(SQLITECPP_RUN_CPPCHECK "Run cppcheck C++ static analysis tool." OFF)
if (SQLITECPP_RUN_CPPCHECK)
    find_program(CPPCHECK_EXECUTABLE NAMES cppcheck)
    if (CPPCHECK_EXECUTABLE)
        # add a cppcheck target to the "all" target
        add_custom_target(SQLiteCpp_cppcheck
         ALL
         COMMAND ${CPPCHECK_EXECUTABLE} -j 8 cppcheck --enable=style --quiet ${CPPCHECK_ARG_TEMPLATE} ${PROJECT_SOURCE_DIR}/src
        )
        execute_process(COMMAND "${CPPCHECK_EXECUTABLE}" --version OUTPUT_VARIABLE CPPCHECK_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Found Cppcheck: ${CPPCHECK_EXECUTABLE} ${CPPCHECK_VERSION}")
    else (CPPCHECK_EXECUTABLE)
        message(STATUS "Could NOT find cppcheck")
    endif (CPPCHECK_EXECUTABLE)
else (SQLITECPP_RUN_CPPCHECK)
    message(STATUS "SQLITECPP_RUN_CPPCHECK OFF")
endif (SQLITECPP_RUN_CPPCHECK)

option(SQLITECPP_RUN_DOXYGEN "Run Doxygen C++ documentation tool." OFF)
if (SQLITECPP_RUN_DOXYGEN)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # add a Doxygen target to the "all" target
        add_custom_target(SQLiteCpp_doxygen
         ALL
         COMMAND doxygen Doxyfile > ${DEV_NULL}
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
    endif (DOXYGEN_FOUND)
else (SQLITECPP_RUN_DOXYGEN)
    message(STATUS "SQLITECPP_RUN_DOXYGEN OFF")
endif (SQLITECPP_RUN_DOXYGEN)

option(SQLITECPP_BUILD_EXAMPLES "Build examples." OFF)
if (SQLITECPP_BUILD_EXAMPLES)
    # add the basic example executable
    add_executable(SQLiteCpp_example1 ${SQLITECPP_EXAMPLES})
    target_link_libraries(SQLiteCpp_example1 SQLiteCpp)
    if (MSYS OR MINGW)
        target_link_libraries(SQLiteCpp_example1 ssp)
    endif ()
else (SQLITECPP_BUILD_EXAMPLES)
    message(STATUS "SQLITECPP_BUILD_EXAMPLES OFF")
endif (SQLITECPP_BUILD_EXAMPLES)

if (SQLITECPP_BUILD_TESTS)
    # add the unit test executable
    add_executable(SQLiteCpp_tests ${SQLITECPP_TESTS})
    target_link_libraries(SQLiteCpp_tests SQLiteCpp)

    find_package(GTest)
    if (GTEST_FOUND)
        message(STATUS "Link to GTest system library")
        target_link_libraries(SQLiteCpp_tests GTest::GTest GTest::Main)
    else (GTEST_FOUND)
        message(STATUS "Compile googletest from source in submodule")
        # deactivate some warnings for compiling the googletest library
        if (NOT MSVC)
            add_compile_options(-Wno-switch-enum)
        endif (NOT MSVC)

        # add the subdirectory containing the CMakeLists.txt for the googletest library
        if (NOT EXISTS "${PROJECT_SOURCE_DIR}/googletest/CMakeLists.txt")
            message(FATAL_ERROR "Missing 'googletest' submodule! Either use 'git submodule init' and 'git submodule update' to get googletest according to the README, or deactivate unit tests with -DSQLITECPP_BUILD_TESTS=OFF")
        endif ()
        add_subdirectory(googletest)
        include_directories("${PROJECT_SOURCE_DIR}/googletest/googletest/include")

        # Add definitions to keep googletest from making the compilation fail
        if (MSVC)
            if (MSVC_VERSION GREATER_EQUAL 1910 AND MSVC_VERSION LESS_EQUAL 1919) # OR MSVC_TOOLSET_VERSION EQUAL 141)
                target_compile_definitions(gtest PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
                target_compile_definitions(gtest_main PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
                target_compile_definitions(gmock PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
                target_compile_definitions(gmock_main PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
            endif (MSVC_VERSION GREATER_EQUAL 1910 AND MSVC_VERSION LESS_EQUAL 1919)
        endif (MSVC)

        target_link_libraries(SQLiteCpp_tests gtest_main)
    endif (GTEST_FOUND)

    # add a "test" target:
    enable_testing()

    # does the tests pass?
    add_test(UnitTests SQLiteCpp_tests)

    if (SQLITECPP_BUILD_EXAMPLES)
        # does the example1 runs successfully?
        add_test(Example1Run SQLiteCpp_example1)
    endif (SQLITECPP_BUILD_EXAMPLES)
else (SQLITECPP_BUILD_TESTS)
    message(STATUS "SQLITECPP_BUILD_TESTS OFF")
endif (SQLITECPP_BUILD_TESTS)

#export(TARGETS sqlite3 NAMESPACE sqlite3:: FILE sqlite3Targets.cmake)
#export(TARGETS SQLiteCpp NAMESPACE SQLiteCpp:: FILE SQLiteCppTargets.cmake)

# API version for SQLiteCpp shared library.
set_property(TARGET SQLiteCpp PROPERTY SOVERSION 0)
