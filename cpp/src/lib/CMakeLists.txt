# Collect all these file groups into one combined list
set(all_files
    ${analysis_files}
    ${error_files}
    ${material_files}
    ${mesh_files}
    ${pre_files}
)

# Build library
add_library(fem_lib ${all_files})
