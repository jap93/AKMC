cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(search.x)

#if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
#  message(FATAL_ERROR "In-source builds are not allowed!!!!!!!!")
#endif()

if (METPOT)
    set(MLPOT OFF)
endif()

if (MLPOT)
    message("Compiling Nudged Elastic Band with ML potentials")
    set(GCC_COVERAGE_COMPILE_FLAGS "-DMLPOT -Wall -pedantic -lm -g")
    #set(GCC_COVERAGE_COMPILE_FLAGS "-DMLPOT -Wall -pedantic -lm -O3 -funroll-loops ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

    include_directories(../../rdforce/schnet ../../common/source)

    add_executable(search.x 
                   neb.cpp 
                   NudgedElasticBand.cpp 
                   ../../rdforce/schnet/Energy.cpp 
                   ../../rdforce/schnet/Field.cpp 
                   ../../rdforce/schnet/NbrListPBC.cpp 
                   ../../common/source/Species.cpp 
                   ../../common/source/MPICommunicator.cpp 
                   ../../common/source/Utility.cpp 
                   ../../common/source/Basis.cpp 
                   ../../common/source/Relax.cpp 
                   ../../common/source/JobControl.cpp 
    )
    
    set_property(TARGET search.x PROPERTY CXX_STANDARD 17)

#    add_subdirectory(/mnt/Progs/workspace/common/source build)

    find_package(Torch REQUIRED)
    find_package(MPI REQUIRED)


    target_link_libraries(search.x PUBLIC MPI::MPI_CXX "${TORCH_LIBRARIES}")

elseif(METPOT)
    message("Compiling Kinetic Monte Carlo:q with metal potentials")
    set(GCC_COVERAGE_COMPILE_FLAGS " -Wall -pedantic -lm -Ofast -fopenmp")
    #set(GCC_COVERAGE_COMPILE_FLAGS " -Wall -pedantic -DDEBUG -lm -O0 -g -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    #set(GCC_COVERAGE_COMPILE_FLAGS " -Wall -pedantic -lm -O0 -g -fsanitize=undefined -fno-omit-frame-pointer")
    #set(GCC_COVERAGE_COMPILE_FLAGS "-O0 -g -fsanitize=memory -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

    include_directories(Common kmc optim potentials potentials/metal)

    add_executable(search.x
                   kmc/kmc.cpp
                   kmc/KineticMonteCarlo.cpp  
                   kmc/Event.cpp  
                   Common/Basis.cpp  
                   Common/MPICommunicator.cpp 
                   Common/JobControl.cpp  
                   Common/RelaxCntrl.cpp
                   Common/ArtCntrl.cpp
                   Common/DimerCntrl.cpp
                   Common/KMCCntrl.cpp 
                   Common/Utility.cpp
                   optim/ART.cpp  
                   optim/Dimer.cpp  
                   optim/Relax.cpp  
                   optim/SaddleSearch.cpp
                   potentials/Species.cpp 
                   potentials/metal/Energy.cpp
                   potentials/metal/Field.cpp
                   potentials/metal/ManyBody.cpp 
                   potentials/metal/MetalPotential.cpp 
                   potentials/metal/NbrListPBC.cpp 
                   )

    #target_include_directories(field_test PUBLIC "/home/jap93/workspace/cuda_mc/common/")

    set_property(TARGET search.x PROPERTY CXX_STANDARD 17)

    find_package(MPI REQUIRED)
    #lapack
    find_package(LAPACK REQUIRED)

    target_link_libraries(search.x PUBLIC MPI::MPI_CXX lapack)
else()

    message("Compiling Nudged Elastic Band with pair potentials")
    set(GCC_COVERAGE_COMPILE_FLAGS " -Wall -pedantic -lm -Ofast")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

    include_directories(Common kmc optim potentials potentials/pair)

    add_executable(search.x
                   kmc/kmc.cpp
                   kmc/KineticMonteCarlo.cpp
                   Common/Basis.cpp
                   Common/MPICommunicator.cpp
                   Common/JobControl.cpp
                   Common/RelaxCntrl.cpp
                   Common/ArtCntrl.cpp
                   Common/DimerCntrl.cpp
                   Common/KMCCntrl.cpp
                   Common/Utility.cpp
                   kmc/Event.cpp  
                   optim/ART.cpp
                   optim/Dimer.cpp
                   optim/Relax.cpp
                   optim/SaddleSearch.cpp
                   potentials/Species.cpp
                   potentials/pair/Energy.cpp
                   potentials/pair/Field.cpp
                   potentials/pair/TwoBody.cpp
                   potentials/pair/Ewald.cpp
                   potentials/pair/NbrListPBC.cpp
                   )

    #target_include_directories(field_test PUBLIC "/home/jap93/workspace/cuda_mc/common/")

    set_property(TARGET search.x PROPERTY CXX_STANDARD 17)

    find_package(MPI REQUIRED)
    #lapack
    find_package(LAPACK REQUIRED)

    target_link_libraries(search.x PUBLIC MPI::MPI_CXX lapack)
endif()

message("Compiling complete")
