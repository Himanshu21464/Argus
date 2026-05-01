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

add_executable(test_partition test_partition.cpp)
target_link_libraries(test_partition PRIVATE argus_core)
add_test(NAME partition COMMAND test_partition)

add_executable(test_hitran test_hitran.cpp)
target_link_libraries(test_hitran PRIVATE argus_core)
add_test(NAME hitran COMMAND test_hitran)

add_executable(test_real_h2o test_real_h2o.cpp)
target_link_libraries(test_real_h2o PRIVATE argus_core)
add_test(NAME real_h2o COMMAND test_real_h2o)

add_executable(test_finite_diff test_finite_diff.cpp)
target_link_libraries(test_finite_diff PRIVATE argus_core)
add_test(NAME finite_diff COMMAND test_finite_diff)

add_executable(test_stress test_stress.cpp)
target_link_libraries(test_stress PRIVATE argus_core)
add_test(NAME stress COMMAND test_stress)

add_executable(test_file_io test_file_io.cpp)
target_link_libraries(test_file_io PRIVATE argus_core)
add_test(NAME file_io COMMAND test_file_io)

add_executable(test_reproducibility test_reproducibility.cpp)
target_link_libraries(test_reproducibility PRIVATE argus_core)
add_test(NAME reproducibility COMMAND test_reproducibility)

add_executable(test_multi_molecule test_multi_molecule.cpp)
target_link_libraries(test_multi_molecule PRIVATE argus_core)
add_test(NAME multi_molecule COMMAND test_multi_molecule)

add_executable(test_self_broadening test_self_broadening.cpp)
target_link_libraries(test_self_broadening PRIVATE argus_core)
add_test(NAME self_broadening COMMAND test_self_broadening)

add_executable(test_voigt_reference test_voigt_reference.cpp)
target_link_libraries(test_voigt_reference PRIVATE argus_core)
add_test(NAME voigt_reference COMMAND test_voigt_reference)

add_executable(test_benchmark test_benchmark.cpp)
target_link_libraries(test_benchmark PRIVATE argus_core)
add_test(NAME benchmark COMMAND test_benchmark)

add_executable(test_guillot test_guillot.cpp)
target_link_libraries(test_guillot PRIVATE argus_core)
add_test(NAME guillot COMMAND test_guillot)

add_executable(test_clouds test_clouds.cpp)
target_link_libraries(test_clouds PRIVATE argus_core)
add_test(NAME clouds COMMAND test_clouds)
