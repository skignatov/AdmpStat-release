# AdmpStat-release
Program for analysis of ADMP trajectories

Content:
AdmpStat.exe - the program itself, compiled with Intel Fortran for Win7/x64/Release ;
admpstat.inp - sample input file for mg7 trajectory processing; 
mg7-admp-6.xyz - trajectory of ADMP created by g16 (coordinates of trajectory points were extracted with admp2xyz.exe);
mg07-classes023-f.xyz - file with the coordinates of reference structures;
admpstat.out - sample output file (results of trajectory processing for mg7);

admp2xyz.exe - program for conversion Gaussian log files to *.xyz ;
admp2xyz.inp - its input file compiled for conversion mg7-admp-6.log_ --> mg7-admp-6.xyz ;
mg7-log.log - sample g16 log-file with admp calculation. (really this is the beginning of trajectory only. The file was shortened due to the github restrictions for file lengths);

src - directory with source code


Workflow:
1. Download the files
2. Run admp2xyz.exe (this converts *.log -> *.xyz)
3. Run AdmpStat.exe (this will process the trajectory. Warning! Execution can take long time, up to several  minutes. A bunch of new files will appear. The main results are in AdmpStat.out)
Edit *.inp files if you wish to change data files or the work options



