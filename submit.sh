##******************##
## ASE longread run ##
##******************##
## Project: Longreads
## Lison Lemoine
## Last modification 10.2023
## Creation: 10.2023
## Submit samples for ASElongread run to LRZ-cluster


######------ Set environment ------######
pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/ge92cem2/00_scripts/longread_ASE/workflow_script.sh"

######------ submit to LRZ ------######

sbatch $pipeline -c "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/ge92cem2/00_scripts/longread_ASE/config.txt"
