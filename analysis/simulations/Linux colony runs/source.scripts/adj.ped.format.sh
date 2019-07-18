#!/usr/bin/env bash

#originally created on: August 21, 2017
#by: Nick Sard

#ABOUT: This script was created to make a slight modification to the BestConfig that comes out of colony
#basically I wanted to remove the * and # used by default in colony and use par1 and par2, so that
#I could work in R with the files with less reading-in issues.

cat *.BestConfig | awk '{print $1,"\t",$2,"\t",$3,"\t",$4}'| sed s/'\*'/'par1_'/g | sed s/'\#'/'par2_'/g > best.config.txt
