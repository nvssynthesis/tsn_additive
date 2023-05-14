TSARA Additive:
an additive synthesizer using the Timbre Space Analysis Resynthesis format output from TimbreSpaceAnalysisResynthesisApp.
This allows for nonlinear playback of an analyzed sound.
That is, rather than stepping through the analysis of a sound from the start frame to the end frame and resynthesizing, navigation can be accomplished based on timbral similarity between frames.
TSARA Additive uses a graph structure built from the analysis file to accomplish this.


Dependencies:
Essentia
JUCE
nvs_libraries
Boost Graph Library
Eigen
Armadillo
MLPack

