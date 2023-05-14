/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "TsaraSynth.h"
#include "Synthesis.h"

//==============================================================================
/**
*/
class Tsara_additiveAudioProcessor  : public foleys::MagicProcessor
,	                                  private juce::AudioProcessorValueTreeState::Listener

                            #if JucePlugin_Enable_ARA
                             , public juce::AudioProcessorARAExtension
                            #endif
{
public:
    //==============================================================================
    Tsara_additiveAudioProcessor();
    ~Tsara_additiveAudioProcessor() override;

	void loadButtonClicked();
	void exportDotButtonClicked();
private:
	void loadFile (const juce::File& file, unsigned int voiceIdx);
public:
    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

	void parameterChanged (const juce::String& param, float value) override;

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
	
    //==============================================================================
	bool hasEditor() const override;
	juce::AudioProcessorEditor * createEditor() override;
    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
	//getStateInformation
	//setStateInformation
private:
	juce::AudioProcessorValueTreeState treeState;
	
    TsaraSynth   synthesiser;
    juce::ValueTree  presetNode;
	
	unsigned int nVoices = 1;
	std::atomic<bool> voicesCreated {false};
	
	juce::File lastFolder = juce::File::getSpecialLocation (juce::File::userDocumentsDirectory);

//	foleys::MagicProcessorState magicProcState { *this };
	
	// really a foleys::MagicPluginEditor
	juce::AudioProcessorEditor 	*ed;
	foleys::MagicGUIState       magicGuiState;
	foleys::MagicGUIBuilder     magicBuilder { magicGuiState };
	std::unique_ptr<juce::FileChooser> chooser = nullptr;
	
//	inline static const juce::String paramLoad		{ "load" };
//	juce::TextButton loadButton { "load" };
//	std::unique_ptr<juce::ButtonParameterAttachment> loadButtonAtchmnt = nullptr;
	
	// GUI MAGIC: define that as last member of your AudioProcessor
    foleys::MagicLevelSource*   outputMeter  = nullptr;
    foleys::MagicPlotSource*    oscilloscope = nullptr;
    foleys::MagicPlotSource*    analyser     = nullptr;
	
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Tsara_additiveAudioProcessor)
};
