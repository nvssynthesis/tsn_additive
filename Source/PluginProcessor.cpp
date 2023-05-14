/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout()
{
    juce::AudioProcessorValueTreeState::ParameterLayout layout;
    TsaraSynth::addADSRParameters (layout);
    TsaraSynth::addAdditiveParameters(layout);
	TsaraSynth::addNavigationParameters(layout);
    TsaraSynth::addGainParameters (layout);
    return layout;
}

//==============================================================================
Tsara_additiveAudioProcessor::Tsara_additiveAudioProcessor()
: foleys::MagicProcessor (juce::AudioProcessor::BusesProperties()
                            .withOutput ("Output", juce::AudioChannelSet::stereo(), true))
,     treeState (*this, nullptr, ProjectInfo::projectName, createParameterLayout())
{
	FOLEYS_SET_SOURCE_PATH(__FILE__);
	
	treeState.replaceState(this->MagicProcessor::createGuiValueTree());
	magicState.setGuiValueTree (treeState.state);
	
	TsaraSynth::TsaraSound::Ptr sound (new TsaraSynth::TsaraSound (treeState));
    synthesiser.addSound (sound);

	for (int i=0; i < nVoices; ++i)	{// 1 voice for now
        synthesiser.addVoice (new TsaraSynth::TsaraVoice (treeState));
		voicesCreated = true;
	}
	jassert(nVoices == synthesiser.getNumVoices());
	
	// MAGIC GUI: add a meter at the output
//	std::unique_ptr<juce::WildcardFileFilter> wcf = std::make_unique<juce::WildcardFileFilter>("*.yml;*.yaml;*.json", "", "Proper analysis file (json or yaml)");
	//    magicBuilder.registerJUCELookAndFeels();
	//    magicBuilder.registerJUCEFactories();
	
//	if (ed)
//		magicBuilder.createGUI (*ed); // need a component

	magicState.addTrigger ("open", [&]
    {
		loadButtonClicked();
    });
	
	outputMeter  = magicState.createAndAddObject<foleys::MagicLevelSource>("output");
    oscilloscope = magicState.createAndAddObject<foleys::MagicOscilloscope>("waveform");
    analyser     = magicState.createAndAddObject<foleys::MagicAnalyser>("analyser");
    magicState.addBackgroundProcessing (analyser);
}
Tsara_additiveAudioProcessor::~Tsara_additiveAudioProcessor()
{
	for (int i=0; i < nVoices; ++i)	// 1 voice for now
		synthesiser.removeVoice (i);
}
void Tsara_additiveAudioProcessor::loadButtonClicked(){
	const unsigned int voiceIdx = 0;	// reminder that we could load to other voices if there were any
	chooser = std::make_unique<juce::FileChooser>
								("Select a proper .json, .yaml, or .yml file...",// const String& dialogBoxTitle
									   juce::File{},// const File& initialFileOrDirectory = File(),
									   "*.json;*.yml;*.yaml",	// const String& filePatternsAllowed = String()
														false,	//  bool useOSNativeDialogBox = true
														false,	// bool treatFilePackagesAsDirectories = false
														 ed);	// Component* parentComponent = nullptr
	auto chooserFlags = juce::FileBrowserComponent::openMode
					  | juce::FileBrowserComponent::canSelectFiles;
	 // Pop up the FileChooser object
	 chooser->launchAsync (chooserFlags, [this] (const juce::FileChooser& fc)
	 {
		 juce::File file = fc.getResult();
		 // if() will succeed if the user actually selects a file (rather than cancelling)
		 if (file != juce::File{})
		 {   // The AudioFormatManager::createReaderFor() function is used attempt to create a reader for this particular file. This will return the nullptr value if it fails (for example the file is not an audio format the AudioFormatManager object can handle).
			 loadFile(file, voiceIdx);
		 }
	 }); // end chooser->launchAsync
}
void Tsara_additiveAudioProcessor::loadFile (const juce::File& file, unsigned int voiceIdx)
{
    lastFolder = file.getParentDirectory();
	juce::String j_fullpath = file.getFullPathName();
	std::string std_fn = j_fullpath.toStdString();

	juce::String fn = file.getFileName();
	juce::String ext = file.getFileExtension();
	if ((ext == ".json") | (ext == ".yaml") | (ext == ".yml")){
		ext = ext.removeCharacters(".");
		std::string std_ext = ext.toStdString();
		
		
		voicesCreated   = false;
		
		
		if(TsaraSynth::TsaraVoice *const voice =
		   dynamic_cast<TsaraSynth::TsaraVoice *>(synthesiser.getVoice(voiceIdx)) )
		{
			voice->loadPool( voice->file2pool(std_fn, std_ext) );
		}
		
		voicesCreated = true;
	}
}
//==============================================================================

bool Tsara_additiveAudioProcessor::hasEditor() const
{
	return true;
}
juce::AudioProcessorEditor* Tsara_additiveAudioProcessor::createEditor()
{
	ed = static_cast<juce::AudioProcessorEditor*>( new foleys::MagicPluginEditor(magicState) );
	return ed;
}
//==============================================================================
const juce::String Tsara_additiveAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool Tsara_additiveAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool Tsara_additiveAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool Tsara_additiveAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double Tsara_additiveAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int Tsara_additiveAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int Tsara_additiveAudioProcessor::getCurrentProgram()
{
    return 0;
}

void Tsara_additiveAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String Tsara_additiveAudioProcessor::getProgramName (int index)
{
    return {};
}

void Tsara_additiveAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void Tsara_additiveAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
	synthesiser.setCurrentPlaybackSampleRate (sampleRate);
	for (int i = 0; i < synthesiser.getNumVoices(); ++i){
		if(TsaraSynth::TsaraVoice *const voice = dynamic_cast<TsaraSynth::TsaraVoice *>(synthesiser.getVoice(i)) ){
			voice->setBlockSize(samplesPerBlock);
			voice->SynthesisContainer::setPlaybackSampleRate(sampleRate);
		}
	}
	
	// MAGIC GUI: setup the output meter
    outputMeter->setupSource (getTotalNumOutputChannels(), sampleRate, 500);
    oscilloscope->prepareToPlay (sampleRate, samplesPerBlock);
    analyser->prepareToPlay (sampleRate, samplesPerBlock);
}

void Tsara_additiveAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool Tsara_additiveAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
	auto outputChanSet = layouts.getMainOutputChannelSet();
	
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void Tsara_additiveAudioProcessor::parameterChanged (const juce::String& param, float value) {
	if (param == IDs::paramNavigationStyle){
		
	}
}


void Tsara_additiveAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
	if (!voicesCreated)
		return;
	
    juce::ScopedNoDenormals noDenormals;
	auto blockSize = buffer.getNumSamples();
	magicState.processMidiBuffer(midiMessages, blockSize, true);
	magicState.updatePlayheadInformation(getPlayHead());
	
	// around here i can figure out how to delay midi buffer for juce synth
	// the point is that essentia's spectral stuff will require a framesize (?) of delay
	// while the adsr etc. does not, so the attack of a new note will carry
	// a frame of the previous note.
	synthesiser.blockwiseCallback();
	synthesiser.renderNextBlock(buffer, midiMessages, 0, blockSize);
	
	buffer.applyGain(2.f);
	
	// MAGIC GUI: send the finished buffer to the level meter
    outputMeter->pushSamples (buffer);
    oscilloscope->pushSamples (buffer);
    analyser->pushSamples (buffer);
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Tsara_additiveAudioProcessor();
}
