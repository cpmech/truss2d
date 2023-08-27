set(TESTS
    z_test_truss2d
)

foreach(var ${TESTS})
    add_executable(${var} "${var}.cpp")
    add_test(NAME ${var} COMMAND ${var})
    target_link_libraries(${var} PUBLIC ${LIBS} ${LIB_NAME})
endforeach(var)
