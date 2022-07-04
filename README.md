# CFD_02_Heat-Convection-and-diffusion

1. Installation and Configuration of MPI in Visual Studio:  
  a) Download: https://www.microsoft.com/en-us/download/details.aspx?id=100593 (perhaps)  
  b) cd ./CFD_02_Heat-Convection-and-diffusion-> Double clicking to open up "_02_Heat-Convection-and-diffusion.sln"  
  c) Project-> Property(P):  
    i. C/C++-> Additional Include Directory-> Edit-> New Line-> C:\Program Files (x86)\Microsoft SDKs\MPI\Include  
    ii. Linker-> Additional Library Directories-> C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86  
    iii. Linker-> Input-> Additional Dependenies-> Edit-> "msmpi.lib"  
2. Define parameters in "header.h"  
3. "Build" the Source.cpp so as to generate .exe file  
4. Execute the .exe file in the CLI: mpiexec -n 4 CFD_02_Heat-Convection-and-diffusion.exe  
5. Post process the "result.csv" via "CFD_Postprocess.m".
