target_sources(ablateLibrary
        PRIVATE
        perfectGas.cpp
        stiffenedGas.cpp
        tChem.cpp
        chemTab.cpp
        chemistryModel.cpp

        PUBLIC
        eos.hpp
        perfectGas.hpp
        tChem.hpp
        stiffenedGas.hpp
        chemTab.hpp
        chemistryModel.hpp
        )

add_subdirectory(transport)
add_subdirectory(tChem)
add_subdirectory(radiationProperties)
