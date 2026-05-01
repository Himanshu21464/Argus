add_executable(test_atmosphere test_atmosphere.cpp)
target_link_libraries(test_atmosphere PRIVATE argus_core)
add_test(NAME atmosphere COMMAND test_atmosphere)

add_executable(test_radiative_transfer test_radiative_transfer.cpp)
target_link_libraries(test_radiative_transfer PRIVATE argus_core)
add_test(NAME radiative_transfer COMMAND test_radiative_transfer)

add_executable(test_ir test_ir.cpp)
target_link_libraries(test_ir PRIVATE argus_core)
add_test(NAME ir COMMAND test_ir)
