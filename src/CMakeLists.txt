cmake_minimum_required(VERSION 3.11)


set(project "main")

add_executable(${project} main.cpp)

# target_include_directories(
#     ${project} PRIVATE 
#       "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen"
#   )


# target_include_directories(
#       ${project} PUBLIC
#       "${PROJECT_SOURCE_DIR}/deps/eigen"
#   )

target_include_directories(
    ${project} PRIVATE 
      "${PROJECT_SOURCE_DIR}/deps/eigen"
  )