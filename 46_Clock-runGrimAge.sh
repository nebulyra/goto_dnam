#! /bin/bash

for i in {1..29} ;
do
printf "../GOTO_Data/Clocks/grimAge/Python_DNAmGrimAge_HorvathCellCounts input=../GOTO_Data/Clocks/grimAge/GrimAgeFile${i}.csv output=../GOTO_Data/Clocks/grimAge/GrimAge_output${i}.csv" | "sh" 
done
