#!/bin/bash

mkdir -p FYSC12-Labs

#intro meeting

#for some reason the pdflatex did not work smoothly.
#cd intro_meeting/reports_presentation_tex/
#pdflatex presentation.tex
#cd ../../
cp intro_meeting/reports_presentation_tex/presentation.pdf FYSC12-Labs/Howto-Reports-Slides.pdf
cp intro_meeting/fysc12_labguidelines.pdf FYSC12-Labs/FYSC12-LabGuidelines.pdf

#error_analysis
cp error_analysis/error_analysis.pdf FYSC12-Labs/

#KF6-Gamma
mkdir -p FYSC12-Labs/KF6-Gamma
cp KF6-Gamma/intro_tex/KF6_intro.pdf FYSC12-Labs/KF6-Gamma/KF6-intro.pdf 
cp KF6-Gamma/procedure_tex/procedure.pdf FYSC12-Labs/KF6-Gamma/KF6-procedure.pdf 
cp KF6-Gamma/KF6-RadionuclideTable-Gamma.pdf FYSC12-Labs/KF6-Gamma/KF6-RadionuclideTable-Gamma.pdf
cp KF6-Gamma/KF6-Attachments.pdf FYSC12-Labs/KF6-Gamma/KF6-Attachments.pdf

#KF7-Beta

#analysis_code
#mkdir -p FYSC12-Labs/analysis_code
cp -r analysis_code FYSC12-Labs/

#cleaning from unnecessary files
rm -rf FYSC12-Labs/analysis_code/__pycache__
rm -rf FYSC12-Labs/analysis_code/.ipynb_checkpoints
rm -rf FYSC12-Labs/analysis_code/KF6-Gamma/__pycache__
rm -rf FYSC12-Labs/analysis_code/KF6-Gamma/.ipynb_checkpoints
rm -rf FYSC12-Labs/analysis_code/FYSC12-Intro/__pycache__
rm -rf FYSC12-Labs/analysis_code/FYSC12-Intro/.ipynb_checkpoints
rm FYSC12-Labs/analysis_code/install_manual.aux
rm FYSC12-Labs/analysis_code/install_manual.log
rm FYSC12-Labs/analysis_code/install_manual.tex

#removing KF7-Beta
rm -rf FYSC12-Labs/analysis_code/KF7-Beta

zip -r FYSC12-Labs FYSC12-Labs

rm -rf FYSC12-Labs
