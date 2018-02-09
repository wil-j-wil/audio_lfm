## LATENT FORCE MODELLING OF AUDIO SUBBAND AMPLITUDE ENVELOPES ##

### Matlab source code to accompany the research paper ‘A Generative Model for Natural Sounds Based on Latent Force Modelling’. [Research explainer web page and paper](http://c4dm.eecs.qmul.ac.uk/audioengineering/wil_j_wil/). ###

Main working script with LFM applied to amplitude envelopes is *audio_lfm_main.m*

### Matlab data structures containing all the data from previously run instances can be found here. Save as a folder called ‘/audio_lfm_data’ in the repo. ###

See the ‘/experiments’ folder for previously run instances (many optimisation settings must be set manually and vary significantly by sound, so use these as a reference).

‘/gen_model’ contains the code to generate novel instances of natural sounds based on the LFM.

The ‘/lfm_augmented’ folder contains the implementation of the theory presented in the research paper - an augmented latent force model that takes into account correlations over many discrete time steps.

Many of these scripts are extended versions of code from the [LFM toolbox](http://becs.aalto.fi/en/research/bayes/lfm/) by Jouni Hartikainen and Simo Sarkka. Unedited scripts from this toolbox are stored in the ‘/lfm_toolbox’ folder.

The ‘/gppad’ folder contains the [probabilistic amplitude demodulation](http://learning.eng.cam.ac.uk/Public/Turner/PAD) algorithm developed by Richard Turner.

Filter bank implementation is from code by [Josh McDermott](http://mcdermottlab.mit.edu/downloads.html).