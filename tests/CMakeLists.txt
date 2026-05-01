add_executable(test_atmosphere test_atmosphere.cpp)
target_link_libraries(test_atmosphere PRIVATE argus_core)
add_test(NAME atmosphere COMMAND test_atmosphere)

add_executable(test_radiative_transfer test_radiative_transfer.cpp)
target_link_libraries(test_radiative_transfer PRIVATE argus_core)
add_test(NAME radiative_transfer COMMAND test_radiative_transfer)

add_executable(test_ir test_ir.cpp)
target_link_libraries(test_ir PRIVATE argus_core)
add_test(NAME ir COMMAND test_ir)

add_executable(test_geometry test_geometry.cpp)
target_link_libraries(test_geometry PRIVATE argus_core)
add_test(NAME geometry COMMAND test_geometry)

add_executable(test_voigt test_voigt.cpp)
target_link_libraries(test_voigt PRIVATE argus_core)
add_test(NAME voigt COMMAND test_voigt)

add_executable(test_line_list test_line_list.cpp)
target_link_libraries(test_line_list PRIVATE argus_core)
add_test(NAME line_list COMMAND test_line_list)

add_executable(test_dual test_dual.cpp)
target_link_libraries(test_dual PRIVATE argus_core)
add_test(NAME dual COMMAND test_dual)
