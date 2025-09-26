cmake_minimum_required(VERSION 3.16)
project(eic_TDR_Q2 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Graf Gpad)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

add_library(eicplot STATIC
    plotting/Plotting.cpp
    utility/Utility.cpp
)
target_link_libraries(eicplot PUBLIC ${ROOT_LIBRARIES})
target_include_directories(eicplot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(ddis_skim_q2 analysis/DDIS_Skim_Q2.cpp)
target_link_libraries(ddis_skim_q2 PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(ddis_plots_q2 analysis/DDIS_Plots_Q2_OOP.cpp)
target_link_libraries(ddis_plots_q2 PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(ddis_skim_q2_xy analysis/DDIS_Skim_Q2_with_xy.cpp)
target_link_libraries(ddis_skim_q2_xy PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(ddis_plots_q2_xy analysis/DDIS_Plots_Q2_OOP_with_xy.cpp)
target_link_libraries(ddis_plots_q2_xy PRIVATE eicplot ${ROOT_LIBRARIES})

#add_executable(ddis_plots_all analysis/DDIS_Plots_Q2_OOP_allvars.cpp)
#target_link_libraries(ddis_plots_all PRIVATE eicplot ${ROOT_LIBRARIES})

#add_executable(ddis_skim_t_xp analysis/DDIS_Skim_t_xP.cpp)
#target_link_libraries(ddis_skim_t_xp PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(ddis_plots_t analysis/DDIS_Plot_t.cpp)
target_link_libraries(ddis_plots_t PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(skim_t analysis/DDIS_Skim_t.cpp)
target_link_libraries(skim_t PRIVATE eicplot ${ROOT_LIBRARIES})

#add_executable(inspect_tree analysis/inspect_tree.cpp)
#target_link_libraries(inspect_tree PRIVATE eicplot ${ROOT_LIBRARIES})

install(TARGETS ddis_plots_q2 ddis_skim_q2 RUNTIME DESTINATION bin)
