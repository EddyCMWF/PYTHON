#!/bin/env bash


for site in AU-Fog    CN-Ha2    CZ-wet    DE-Akm    DE-SfN    DE-Spw    DK-NuF    DK-ZaF   FI-Lom    NO-Adv    RU-Che    US-Los    US-Myb    US-ORv    US-Tw1    US-Tw4
do
    echo python2.7 extract_met_data.py $site 
    python2.7 extract_met_data.py $site
 
done


