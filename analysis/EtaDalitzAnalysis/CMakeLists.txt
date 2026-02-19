cmake_minimum_required(VERSION 3.14)
project(plotting)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# get compile_commands.json for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree Graf Gpad MathCore)


# Add executable
add_executable(FullEtaDalitzAnalysis
  FullEtaDalitzAnalysis.cpp
)

# Optional: warnings and debug symbols
target_compile_options(FullEtaDalitzAnalysis PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wimplicit-fallthrough -g)

target_include_directories(FullEtaDalitzAnalysis
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  SYSTEM
    ${ROOT_INCLUDE_DIRS}
)



# Link against modern ROOT targets
target_link_libraries(FullEtaDalitzAnalysis
  PRIVATE
  ROOT::Core
  ROOT::Hist
  ROOT::Tree
  ROOT::Graf
  ROOT::Gpad
  ROOT::MathCore
)

