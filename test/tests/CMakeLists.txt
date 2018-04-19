file(GLOB Test_SRC *.cpp)
set(MAIN_TEST "main.cpp")
list(REMOVE_ITEM Test_SRC "${MAIN_TEST}")

set(TEST_BED "${CMAKE_CURRENT_SOURCE_DIR}/data/test")
set(TOY_BED "${CMAKE_CURRENT_SOURCE_DIR}/data/toy")

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/log_out")
configure_file(test_config.h.in ${CMAKE_BINARY_DIR}/generated/test_config.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

function(addTestItem test_name source lib_name run_param)
    ADD_EXECUTABLE(${test_name} ${MAIN_TEST} ${source})
    TARGET_LINK_LIBRARIES(${test_name}
        ${lib_name} 
        libgtest
        libgmock
        -Wl,-Bstatic z -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp -Wl,-Bdynamic pthread m dl
        
    )
    add_test(${test_name} ${test_name} "${run_param}")
endfunction(addTestItem)

#addTestItem(hello_test "hello_test.cpp" "" "")
#addTestItem(table_test table_test.cpp tables "")
#addTestItem(logger_test test_logger.cpp logger "")
#addTestItem(thread_test test_thread.cpp "threadpool;logger" "")
#addTestItem(pheno_test test_pheno.cpp "pheno;logger" "")
#addTestItem(marker_test test_marker.cpp "marker;logger" "")
#addTestItem(buffer_test test_buffer.cpp "logger" "")
#addTestItem(geno_test test_geno.cpp "logger;geno;marker;pheno;tables" "")
#addTestItem(grm_test test_grm.cpp "logger;grm;geno;marker;pheno;tables;threadpool" "")
addTestItem(chisq_test test_chisq.cpp "statlib" "")
addTestItem(covar_test test_covar.cpp "covar" "")
