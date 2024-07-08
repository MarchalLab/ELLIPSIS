add_executable(ELLIPSIS main.cpp settings.cpp settings.h util.cpp util.h spliceGraph.cpp
        spliceGraph.h sgnode.cpp sgnode.h sgedge.cpp sgedge.h geneInfo.cpp geneInfo.h
        pool.cpp pool.h GTFReader.cpp GTFReader.h SJCounter.cpp SJCounter.h observedDepth.cpp observedDepth.h
        PSIsolver.cpp PSIsolver.h var.cpp var.h NNLSSolver.cpp NNLSSolver.h)

target_link_libraries(ELLIPSIS)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(ELLIPSIS Threads::Threads)

target_link_libraries(ELLIPSIS stdc++fs)

#Eigen is a template library => no libraries to find
find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
        PATH_SUFFIXES eigen3)

include_directories(${EIGEN3_INCLUDE_DIR})

install(TARGETS ELLIPSIS RUNTIME DESTINATION bin)