#!/bin/bash
echo "! RELEASED ON "$(date +%d_%b_%Y)" AT "$(date +%H:%M) > DATx34.data
cat DATx34.data Module_SETTINGS.f90 Module_Angles.f90 Module_Molecular.f90 Module_Debug.f90 Module_Autocorrelation.f90 Module_Diffusion.f90 Module_Recognition.f90 Module_Distribution.f90 Module_Distances.f90 Module_Speciation.f90 Module_Cluster.f90 Module_Main.f90 > pre-alpha.f03
sed -i -e 's/\t/ /g' pre-alpha.f03
sed -i -e 's/!RELEASEDATE!/'$(date +%d_%b_%Y)'/g' pre-alpha.f03
sed -i -e 's/!RELEASEYEAR!/'$(date +%Y)'/g' pre-alpha.f03
sed -i -e 's/DEVELOPERS_VERSION_DEFAULT=.TRUE./DEVELOPERS_VERSION_DEFAULT=.FALSE./g' pre-alpha.f03
rm DATx34.data
zip modules_separate.zip Module_SETTINGS.f90 Module_Angles.f90 Module_Molecular.f90 Module_Debug.f90 Module_Autocorrelation.f90 Module_Diffusion.f90 Module_Recognition.f90 Module_Distribution.f90 Module_Distances.f90 Module_Speciation.f90 Module_Cluster.f90 Module_Main.f90
echo "program currently has "$(cat pre-alpha.f03 | wc -l)" lines."
mv pre-alpha.f03 ../
