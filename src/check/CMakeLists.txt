set(TESTS
    z_test_equal_scalars
    z_test_equal_vectors
)

foreach(var ${TESTS})
    add_executable(${var} "${var}.cpp")
    add_test(NAME ${var} COMMAND ${var})
endforeach(var)
