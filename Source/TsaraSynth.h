/*
  ==============================================================================

    TsaraSynth.h
    Created: 9 Mar 2023 10:15:23am
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Synthesis.h"
#include "graph.h"

namespace IDs
{

    static juce::String paramAttack  	{ "attack" };
    static juce::String paramDecay   	{ "decay" };
    static juce::String paramSustain 	{ "sustain" };
    static juce::String paramRelease 	{ "release" };
    static juce::String paramGain    	{ "gain" };

    static juce::String paramTranspose	{ "transpose" };
    static juce::String paramTilt    	{ "tilt" };
    static juce::String paramSpecStretch{ "spectral stretch" };
    static juce::String paramShift    	{ "frequency shift" };

	static juce::String paramLocation	{ "location" };
	static juce::String paramTraversalSpeed	{ "traversal speed" };
	static juce::String param_C_kernel { "gaussian width" };
	static juce::String paramProbPower { "probability power" };
	static juce::String paramNavigationStyle 	{ "wander" };
}
enum class navigationTypes_e {
	adjacent = 0,
	wander
};
static const std::map<navigationTypes_e, juce::String> navTypeToString{
	{navigationTypes_e::adjacent, 	"adjacent"},
	{navigationTypes_e::wander, 	"wander"}
};

class TsaraSynth
: 	public juce::Synthesiser
{
public:
    static int  numOscillators;

    static void addADSRParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout);
    static void addAdditiveParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout);
    static void addGainParameters (juce::AudioProcessorValueTreeState::ParameterLayout& layout);
	static void addNavigationParameters(juce::AudioProcessorValueTreeState::ParameterLayout& layout);
	
    TsaraSynth() = default;
	
	
	
	void blockwiseCallback();
	
    class TsaraSound : public juce::SynthesiserSound
    {
    public:
        TsaraSound (juce::AudioProcessorValueTreeState& state);
        bool appliesToNote (int) override { return true; /* all notes */}
        bool appliesToChannel (int) override { return true; /* all channels */}

        juce::ADSR::Parameters getADSR() const;
		
    private:
        juce::AudioProcessorValueTreeState& state;
        juce::AudioParameterFloat* attack  = nullptr;
        juce::AudioParameterFloat* decay   = nullptr;
        juce::AudioParameterFloat* sustain = nullptr;
        juce::AudioParameterFloat* release = nullptr;
		
        juce::AudioParameterFloat* soundGain    = nullptr;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TsaraSound)
    };
	unsigned int getTotalFramesSamplewiseIdx(){
		if (auto v = dynamic_cast<TsaraVoice*>(voices[0])){
			return v->getTotalFramesSamplewiseIdx();
		}
		return 0;
	}
	class TsaraVoice : public juce::SynthesiserVoice, public nvs::SynthesisContainer
    {
    public:
        TsaraVoice (juce::AudioProcessorValueTreeState& state);
		
        bool canPlaySound (juce::SynthesiserSound *) override;
        void startNote (int midiNoteNumber,
                        float velocity,
                        juce::SynthesiserSound* sound,
                        int currentPitchWheelPosition) override;
        void stopNote (float velocity, bool allowTailOff) override;
        void pitchWheelMoved (int newPitchWheelValue) override;
        void controllerMoved (int controllerNumber, int newControllerValue) override;
        void renderNextBlock (juce::AudioBuffer<float>& outputBuffer,
                              int startSample,
                              int numSamples) override;
        void setCurrentPlaybackSampleRate (double newRate) override;
		
		void loadPool(const essentia::Pool &p) override;
    private:
		nvs::tsaraCommon::sineModelTimbre getSineModelTimbre();
		inline float getTraversalSpeed() const;
		inline navigationTypes_e getNavigationType() const;

		
		struct DebugCounter{
			unsigned int dbg_idx {0};
			inline unsigned int incr(const unsigned int modulo){
				++dbg_idx;
				dbg_idx %= modulo;
				return dbg_idx;
			}
		}  debugCounter;
		
		class NavigationPhasor{
		// when its phasor reaches peak triggers new frame
		private:
			float frequency {3.7f};
			float phase {0.0f};
			nvs::SynthesisContainer &owner;
		public:
			NavigationPhasor(nvs::SynthesisContainer &newOwner):	owner(newOwner) {}
			[[nodiscard]]
			bool increment(){
				bool trigger { false };
				phase += (frequency / owner.getPlaybackSampleRate());
				while (phase >= 1.0){
					trigger = true;
					phase -= 1.0;
				}
				return trigger;
			}
			void setFrequency(float newFreq){
				float nyq = owner.getPlaybackSampleRate() * 0.5f;
				if (newFreq < 0.f){
					newFreq = 0.f;
				} else if (newFreq > nyq) {
					newFreq = nyq;
				}
				frequency = newFreq;
			}
			inline float getFrequency() const{
				return frequency;
			}
			inline float getPhase() const{
				return phase;
			}
		};
		NavigationPhasor navPhasor;
		bool shouldSelectNewFrame {false};
		unsigned long initialTargetFrame { 0 };
		
		std::unique_ptr<nvs::sgt::DirectedGraph_t> dg = nullptr;
		nvs::sgt::DirectedGraph_t::vertex_iterator initialTargetVit ;
		
        double getDetuneFromPitchWheel (int wheelValue) const;
        double getFrequencyForNote (int noteNumber, double detune, double concertPitch = 440.0) const;

        /*void updateFrequency (BaseOscillator& oscillator, bool noteStart = false);
        std::vector<std::unique_ptr<BaseOscillator>> oscillators;*/

        double                      pitchWheelValue = 0.0;
        int                         maxPitchWheelSemitones = 12;
        const int                   internalBufferSize = 64;
//        juce::AudioBuffer<float>    oscillatorBuffer;
        juce::AudioBuffer<float>    voiceBuffer;
        juce::ADSR                  adsr;
		
		juce::AudioParameterFloat* voiceTranspose 	= nullptr;
		juce::AudioParameterFloat* voiceTilt 		= nullptr;
		juce::AudioParameterFloat* voiceStretch 	= nullptr;
		juce::AudioParameterFloat* voiceShift 		= nullptr;
		
		juce::AudioParameterFloat* voiceNavigation	= nullptr;
		juce::AudioParameterFloat* voiceTraversalSpeed 	= nullptr;
		juce::AudioParameterChoice* voiceNavigationStyle= nullptr;
		juce::AudioParameterFloat* voiceProbPower	= nullptr;
		juce::AudioParameterFloat* voice_C_kernelScaling= nullptr;

		juce::AudioParameterFloat*  voiceGain = nullptr;
        float                       lastVoiceGain = 0.0;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TsaraVoice)
    };

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TsaraSynth)
};
