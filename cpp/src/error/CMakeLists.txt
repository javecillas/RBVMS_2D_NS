file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(error_files ${headers} ${sources} PARENT_SCOPE)
