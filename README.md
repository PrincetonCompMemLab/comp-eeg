# comp-eeg
Code for the Competitive Retrieval EEG project

Contains all code required to run and analyze data from the paper:
    
    Reductions in Retrieval Competition Predict the Benefit of Repeated Testing
        Nicole S. Rafidi, Justin C. Hulbert, Paula P. Brooks, and Kenneth A. Norman

Paper can be found here:  https://rdcu.be/35V3

Code was written and run in MATLAB 2014b

PLEASE NOTE: due to the fact that random seeds were assigned in shuffle mode for nested cross-validation during classifier training, permutation testing, and boostrapping, it may not be possible to perfectly replicate results in paper.

NIMH Data Archive: https://ndar.nih.gov/study.html?id=604
Additional required reaction time data: https://drive.google.com/open?id=17vvkGS5ZUpUvWCDIZFDiOKmaF9rQh2E2

*Contents*

Experiment: Directory containing files needed to run the experiment (i.e. present the stimuli to participants). Code was written for Psychtoolbox Version 3.0.12 - Flavor: beta - Corresponds to SVN Revision 6094.
    
    make_*_stim.m will create stimuli for the various parts of the experiment. 
    Run_exp.m will run the selected experiment (Session 1 or 2, study or retrieval phase), chosen via dialog box.

Preprocessing: Directory containing files for preprocessing EEG data. Code was written for EEGLAB Version 13.4.4b. 
    
    preprocPipeline*.m will run the full pipeline for a given subject, producing exactly what is needed for classification
    createAnswerTraj.m will collect the in-experiment responses to Session 2 stimuli (necessary for classification)

Classification: Directory containing files to generate analyses and figures from paper.
    
    plotBehaveResults.m: Creates figures 2, 3, and 4
    runCompClass_PermTest.m, aggPermResultsComp.m, plotPermResultsComp.m: Create figures 5, S4a
    runCompClass_Models.m, plotModelRDM.m, weight2map.m: Create figures 6 and S1
    createTGM.m, plotTGM.m, plotTGM_subZ.m: Create figures 7, S2, S3, S4b-f
