# Helper macro to reduce boilerplate for your many utility apps
macro(add_hdgeant4_util UTIL_NAME)
    add_executable(${UTIL_NAME} ${UTIL_NAME}.cc)
    target_link_libraries(${UTIL_NAME} 
        hdgeant4 
        ${ROOT_LIBRARIES} 
        # Add DANA libs manually or via a FindDANA.cmake if available
        -L$ENV{HALLD_RECON_HOME}/$ENV{BMS_OSNAME}/lib -lHDGEOMETRY -lDANA
    )
endmacro()

add_hdgeant4_util(beamtree)
add_hdgeant4_util(genBH)
add_hdgeant4_util(geneBH)
add_hdgeant4_util(adapt)
add_hdgeant4_util(samplesep)
