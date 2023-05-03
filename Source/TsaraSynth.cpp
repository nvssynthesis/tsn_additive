/*
  ==============================================================================

    TsaraSynth.cpp
    Created: 9 Mar 2023 10:15:23am
    Author:  Nicholas Solem

  ==============================================================================
*/

#include "TsaraSynth.h"
#define USING_JUCE_ADSR 1



//==============================================================================

int TsaraSynth::numOscillators = 8;

void TsaraSynth::addADSRParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout)
{
    auto attack  = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID (IDs::paramAttack, 1),  "Attack",  juce::NormalisableRange<float> (0.0001f, 1.f, 0.01f), 0.10f);
    auto decay   = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID (IDs::paramDecay, 1),   "Decay",   juce::NormalisableRange<float> (0.001f, 1.f, 0.01f), 0.10f);
    auto sustain = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID (IDs::paramSustain, 1), "Sustain", juce::NormalisableRange<float> (0.0f,   1.0f, 0.01f), 1.0f);
    auto release = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID (IDs::paramRelease, 1), "Release", juce::NormalisableRange<float> (0.0001f, 1.f, 0.01f), 0.10f);

    auto group = std::make_unique<juce::AudioProcessorParameterGroup>("adsr", "ADSR", "|",
                                                                      std::move (attack),
                                                                      std::move (decay),
                                                                      std::move (sustain),
                                                                      std::move (release));
    layout.add (std::move (group));
}

void TsaraSynth::addAdditiveParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout)
{
    auto group = std::make_unique<juce::AudioProcessorParameterGroup>("additive", "Additive", "|");

	using namespace nvs;
	
	auto smtParamMap = tsaraCommon::sineModelTimbre::pMap;

	auto grab = [&](nvs::tsaraCommon::sineModelTimbre::smTimbreEnum en, tsaraCommon::baseTimbre::pMapIdx idx){
		return smtParamMap.at(en)[idx];
	};
	auto grabTranspose = [&](tsaraCommon::baseTimbre::pMapIdx idx){
		return grab(tsaraCommon::sineModelTimbre::smTimbreEnum::transpose_e, idx);
	};
	auto grabTilt = [&](tsaraCommon::baseTimbre::pMapIdx idx){
		return grab(tsaraCommon::sineModelTimbre::smTimbreEnum::tilt_e, idx);
	};
	auto grabStretch = [&](tsaraCommon::baseTimbre::pMapIdx idx){
		return grab(tsaraCommon::sineModelTimbre::smTimbreEnum::stretch_e, idx);
	};
	auto grabShift = [&](tsaraCommon::baseTimbre::pMapIdx idx){
		return grab(tsaraCommon::sineModelTimbre::smTimbreEnum::shift_e, idx);
	};
	
	auto semitAttr = juce::AudioParameterFloatAttributes().withStringFromValueFunction ([] (auto x, auto) { return juce::String (x); }).withLabel ("semitones");
	auto transpose = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramTranspose, 1), "Transpose", juce::NormalisableRange<float> (	grabTranspose(tsaraCommon::baseTimbre::pMapIdx::min_e),
											grabTranspose(tsaraCommon::baseTimbre::pMapIdx::max_e), 0.0f),
											grabTranspose(tsaraCommon::baseTimbre::pMapIdx::def_e), semitAttr);	// in semitones
	auto tilt = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramTilt, 1), "Spectral Tilt", juce::NormalisableRange<float> ( 	grabTilt(tsaraCommon::baseTimbre::pMapIdx::min_e),
											grabTilt(tsaraCommon::baseTimbre::pMapIdx::max_e), 0.0f),
											grabTilt(tsaraCommon::baseTimbre::pMapIdx::def_e));
	auto stretch = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramSpecStretch, 1), "Spectral Stretch",
	   juce::NormalisableRange<float> (		grabStretch(tsaraCommon::baseTimbre::pMapIdx::min_e),
											grabStretch(tsaraCommon::baseTimbre::pMapIdx::max_e), 0.0f),
											grabStretch(tsaraCommon::baseTimbre::pMapIdx::def_e));
	auto hzAttr = juce::AudioParameterFloatAttributes().withStringFromValueFunction ([] (auto x, auto) { return juce::String (x); }).withLabel ("Hz");
	auto shift = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramShift, 1), "Spectral Shift", juce::NormalisableRange<float> (	grabShift(tsaraCommon::baseTimbre::pMapIdx::min_e),
											grabShift(tsaraCommon::baseTimbre::pMapIdx::max_e), 0.0f),
											grabShift(tsaraCommon::baseTimbre::pMapIdx::def_e), hzAttr);	// in hertz
	group->addChild(std::move(transpose));
	group->addChild(std::move(tilt));
	group->addChild(std::move(stretch));
	group->addChild(std::move(shift));

	layout.add (std::move (group));
}
void TsaraSynth::addNavigationParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout)
{
    auto group = std::make_unique<juce::AudioProcessorParameterGroup>("navigation", "Navigation", "|");

	auto location = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramLocation, 1), "Location", juce::NormalisableRange<float>(0.f,		// min
										1.f,	// max
										0.f),	// spacing
										0.f);	// default
	auto traversalSpeed = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramTraversalSpeed, 1), "Traversal Speed", juce::NormalisableRange<float>(0.f,		// min
												20.f,	// max
												0.f,	// spacing
												0.5f,	//skewFactor
												false), // useSymmetricSkew
											0.f);	// default
	
	auto gaussRange = juce::NormalisableRange<float> (0.01f, // min
												100.f, // max
												0.f,	//spacing
												0.1f,	// skewFactor
													  false);// useSymmetricSkew
	gaussRange.setSkewForCentre(1.f);
	auto gaussianKernel = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::param_C_kernel, 1),
		  "Gaussian Kernel Scaling", gaussRange,
											  1.f);	// default
	
	auto probPowerRange = juce::NormalisableRange<float> (-100.f, // min
												100.f, // max
												0.f,	//spacing
												0.1f,	// skewFactor
													  true);// useSymmetricSkew
	probPowerRange.setSkewForCentre(0.f);
	auto probPower = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID(IDs::paramProbPower, 1),
																 "Probability Power", probPowerRange,
																 1.f);	// default
	
	auto wander = std::make_unique<juce::AudioParameterChoice>(juce::ParameterID(IDs::paramNavigationStyle, 1), //param ID
															   "Navigation Style",
															   juce::StringArray (navTypeToString.at((navigationTypes_e)0),
																				  navTypeToString.at((navigationTypes_e)1)),
															   0);	// default choice
	group->addChild(std::move(location));
	group->addChild(std::move(traversalSpeed));
	group->addChild(std::move(wander));
	group->addChild(std::move(probPower));
	group->addChild(std::move(gaussianKernel));
	layout.add(std::move(group));
}
void TsaraSynth::addGainParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout)
{
    auto gain  = std::make_unique<juce::AudioParameterFloat>(juce::ParameterID (IDs::paramGain, 1), "Gain", juce::NormalisableRange<float> (0.0f, 8.0f, 0.001f), 0.70f);

    layout.add (std::make_unique<juce::AudioProcessorParameterGroup>("output", "Output", "|", std::move (gain)));
}
void TsaraSynth::blockwiseCallback(){
	for (auto i = 0; i < voices.size(); ++i){
		if (auto v = dynamic_cast<TsaraVoice*>(voices[i])){
			v->prepareBlock();
		}
	}
}

//==============================================================================

TsaraSynth::TsaraSound::TsaraSound (juce::AudioProcessorValueTreeState& stateToUse)
  : state (stateToUse)
{
    attack = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramAttack));
    jassert (attack);
    decay = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramDecay));
    jassert (decay);
    sustain = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramSustain));
    jassert (sustain);
    release = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramRelease));
    jassert (release);
	
    soundGain = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramGain));
    jassert (soundGain);
}

juce::ADSR::Parameters TsaraSynth::TsaraSound::getADSR() const
{
    juce::ADSR::Parameters parameters;
    parameters.attack  = attack->get();
    parameters.decay   = decay->get();
    parameters.sustain = sustain->get();
    parameters.release = release->get();
    return parameters;
}


//==============================================================================

TsaraSynth::TsaraVoice::TsaraVoice (juce::AudioProcessorValueTreeState& state)
:	navPhasor(*this)
{
//	const std::string fn =  "/Users/nicholassolem/development/audio for analysis/blood is life.yaml";
	const std::string fn = "/Users/nicholassolem/development/audio for analysis/Nannou_clip_mono.yaml";
	// "/Users/nicholassolem/development/audio for analysis/sine, saw, noise.yaml"
	essentia::Pool const p = file2pool(fn, "yaml");
	loadPool(p);

	voiceGain = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramGain));
    jassert (voiceGain);
	voiceTranspose = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramTranspose));
	jassert(voiceTranspose);
	voiceTilt = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramTilt));
	jassert(voiceTilt);
	voiceShift = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramShift));
	jassert(voiceShift);
	voiceStretch = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramSpecStretch));
	jassert(voiceStretch);
	voiceNavigation = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramLocation));
	jassert(voiceNavigation);
	voiceTraversalSpeed = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter (IDs::paramTraversalSpeed));
    jassert (voiceTraversalSpeed);
	voiceNavigationStyle = dynamic_cast<juce::AudioParameterChoice*>(state.getParameter(IDs::paramNavigationStyle));
	jassert (voiceNavigationStyle);
	voice_C_kernelScaling = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::param_C_kernel));
	jassert(voice_C_kernelScaling);
	voiceProbPower = dynamic_cast<juce::AudioParameterFloat*>(state.getParameter(IDs::paramProbPower));
	jassert(voiceProbPower);
	
    voiceBuffer.setSize (1, internalBufferSize);
}
void TsaraSynth::TsaraVoice::loadPool(const essentia::Pool &p){
	std::cout << "calling base loadPool\n";
	nvs::SynthesisContainer::loadPool(p);
	std::cout << "also loading into graph\n";
	
	nvs::sound_representation soundRep = nvs::getSoundRepresentationFromPool(p);
//	auto PCAmat = soundRep.getPCAmat();
//	float maxDistanceForConnection = nvs::getMaxTimbralDistanceForConnection(PCAmat, 0.15f, 100UL);
	dg = std::make_unique<nvs::DirectedGraph_t>( nvs::createGraphFromSoundRepresentation(soundRep) );
	
	std::cout << "num vert before remove: " << boost::num_vertices(*dg) << '\n';
	std::cout << "num edge before remove: " << boost::num_edges(*dg) << '\n';
//	nvs::removeProportionOfVertices(*dg, 0.1);
	std::cout << "num vert after remove: " << boost::num_vertices(*dg) << '\n';
	std::cout << "num edge after remove: " << boost::num_edges(*dg) << '\n';
	exportGraphAsDot(*dg, "myGraph.dot");
}

bool TsaraSynth::TsaraVoice::canPlaySound (juce::SynthesiserSound* sound)
{
    return dynamic_cast<TsaraSound*>(sound) != nullptr;
}

void TsaraSynth::TsaraVoice::startNote (int midiNoteNumber,
                                          float velocity,
                                          juce::SynthesiserSound* sound,
                                          int currentPitchWheelPosition)
{
    juce::ignoreUnused (midiNoteNumber, velocity);

	float newF0 = midi2hz(midiNoteNumber);
	this->SynthesisContainer::setPitch(newF0);
	this->SynthesisContainer::resetPhase(true);
	
	const float initialLoc = voiceNavigation->get();
	initialTargetFrame = round(initialLoc * (getNumFrames() - 1));
	initialTargetVit = nvs::vertexDescriptorToIterator(*dg, nvs::DirectedGraph_t::vertex_descriptor(initialTargetFrame));

	if (auto* tsaraSound = dynamic_cast<TsaraSound*>(sound)){
		adsr.setParameters(tsaraSound->getADSR());
	}
	
    pitchWheelValue = getDetuneFromPitchWheel (currentPitchWheelPosition);
#if USING_JUCE_ADSR
    adsr.noteOn();
#endif
}

void TsaraSynth::TsaraVoice::stopNote (float velocity,
                                         bool allowTailOff)
{
    juce::ignoreUnused (velocity);

	
#if USING_JUCE_ADSR
    if (allowTailOff)
    {
        adsr.noteOff();
		std::cout << "no tail off, setting note off\n";
    }
    else
    {
//        adsr.reset();		// !!! sets env val to 0
        clearCurrentNote();
		std::cout << "yes tail off, resetting and clearing current note\n";
    }
#endif
}

void TsaraSynth::TsaraVoice::pitchWheelMoved (int newPitchWheelValue)
{
    pitchWheelValue = getDetuneFromPitchWheel (newPitchWheelValue);
}

void TsaraSynth::TsaraVoice::controllerMoved (int controllerNumber, int newControllerValue)
{
    juce::ignoreUnused (controllerNumber, newControllerValue);
}

nvs::tsaraCommon::sineModelTimbre TsaraSynth::TsaraVoice::getSineModelTimbre(){
	nvs::tsaraCommon::sineModelTimbre timbre;
	timbre.tranposeFreqs = voiceTranspose->get();
	timbre.stretchFreqs = voiceStretch->get();
	timbre.shiftFreqs = voiceShift->get();
	timbre.tiltFreqs = voiceTilt->get();
	
	return timbre;
}

void TsaraSynth::TsaraVoice::renderNextBlock (juce::AudioBuffer<float>& outputBuffer,
                                                int startSample,
                                                int numSamples) {
	if (debugCounter.incr(1) == 0){
		
	}
	if (! adsr.isActive()){
        return;
	}
	careAboutPlotting = true;

	auto jChan = outputBuffer.getNumChannels();
	auto sChan = getNumChans();

	setTimbre( getSineModelTimbre() );
	navPhasor.setFrequency(getTraversalSpeed());

	float C_kernelScaling = voice_C_kernelScaling->get();
	double probPower = static_cast<double>(voiceProbPower->get());

	if (shouldSelectNewFrame){
		nvs::DirectedGraph_t::vertex_descriptor vd = nvs::traverseToRandomVertex(*dg, *initialTargetVit, C_kernelScaling, probPower);
		initialTargetVit = nvs::vertexDescriptorToIterator(*dg, nvs::DirectedGraph_t::vertex_descriptor(vd));

		jassert(vd != getNumFrames() );
		initialTargetFrame = vd;
		requestNextTargetFrame(initialTargetFrame); // sets isCaughtUpToCurrentFrame(FALSE)

//		std::cout << "targetFrame: "<< initialTargetFrame << '\n';
		shouldSelectNewFrame = false;
	}
	
	nvs::mcSamples samps;
	for (auto i = 0; i < numSamples; ++i){
		if (navPhasor.increment()){
			shouldSelectNewFrame  = true;
			// only set true here. it is set false non-samplewise (above this loop)
		}

		produceSamples(samps);
		float val = *samps.L;
		auto destIdx = startSample + i;
		outputBuffer.setSample(0, destIdx, val);
		if (jChan == 2){
			if (sChan == 2)
				outputBuffer.setSample(1, destIdx, val);//*samps.R);
			else if (getNumChans() == 1)
				outputBuffer.setSample(1, destIdx, val);//*samps.L);
		}
	}
	adsr.applyEnvelopeToBuffer (outputBuffer, startSample, numSamples);

	const auto gain = voiceGain->get();
	outputBuffer.applyGainRamp (startSample,// int startSample,
			numSamples,					// int numSamples,
			lastVoiceGain,    			// Type startGain,
			gain);					// Type endGain)
	lastVoiceGain = gain;

	if (! adsr.isActive()){
		std::cout << "note cleared at end of render\n";
		clearCurrentNote();
	}

}

void TsaraSynth::TsaraVoice::setCurrentPlaybackSampleRate (double newRate)
{
    juce::SynthesiserVoice::setCurrentPlaybackSampleRate (newRate);
	
    juce::dsp::ProcessSpec spec;
    spec.sampleRate = newRate;
    spec.maximumBlockSize = juce::uint32 (internalBufferSize);
    spec.numChannels = 1;
    /*for (auto& osc : oscillators)
        osc->osc.prepare (spec);*/
	
	this->SynthesisContainer::setPlaybackSampleRate(newRate);
}

double TsaraSynth::TsaraVoice::getFrequencyForNote (int noteNumber, double detune, double concertPitch) const
{
    return concertPitch * std::pow (2.0, (noteNumber + detune - 69.0) / 12.0);
}

double TsaraSynth::TsaraVoice::getDetuneFromPitchWheel (int wheelValue) const
{
    return (wheelValue / 8192.0) - 1.0;
}
float TsaraSynth::TsaraVoice::getTraversalSpeed() const {
	return voiceTraversalSpeed->get();
}
navigationTypes_e TsaraSynth::TsaraVoice::getNavigationType() const {
	auto choice = voiceNavigationStyle->getIndex();
	jassert((int(choice) == 0) || (int(choice) == 1));
	return (navigationTypes_e)choice;
}
