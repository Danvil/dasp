include_directories(
)

add_executable(pds_cmd main.cpp)

target_link_libraries(pds_cmd
	pds
	boost_program_options
	boost_timer
	boost_system
)

