set(EXAMPLES
    ex_three_member_truss
)

foreach(var ${EXAMPLES})
    add_executable(${var} "${var}.cpp")
    target_link_libraries(${var} PUBLIC ${LIB_NAME})

    # this will run the example after build
    add_custom_target(run_${var} ALL ${var})
endforeach(var)
