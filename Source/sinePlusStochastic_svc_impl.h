//
//  sinePlusStochastic_svc_impl.h
//  tsara_additive
//
//  Created by Nicholas Solem on 5/22/23.
//  Copyright © 2023 nvssynthesis. All rights reserved.
//

#ifndef sinePlusStochastic_svc_impl_h
#define sinePlusStochastic_svc_impl_h
#include "Synthesis.h"


SynthesisContainer::sinePlusStochastic_svc::sinePlusStochastic_svc(const SynthesisContainer &owningSynth, unsigned int i)
:   									abstractSynthVoiceChannel(owningSynth, i)
{
	_voice_timbre = sinePlusStochasticTimbre::factory();
	
	owner.validateSubclassOwner();
	_state.isInitialized = true;
}

/*
 what we've got to work with
realVec_t _sineFreqs;
realVec_t _sineMags;
realVec_t _sinePhases;
realVec_t _stocEnv;

realVec_t outputFrame;
*/

void SynthesisContainer::sinePlusStochastic_svc::doInit(){
	_spsModelAlgo = owner.ess->factory.create("SpsModelSynth",
											  "fftSize", owner.sc_analysisData->frameSize,
											  "hopSize", owner.sc_analysisData->hopSize,
											  "sampleRate", owner.sc_analysisData->sampleRate,
											  "stocf", 0.1f);
	unsigned int nHarmonics = owner.sc_analysisData->numHarmonics;

	_sineFreqs.reserve(nHarmonics * 2);
	_sineMags.reserve(nHarmonics * 2);
	_sinePhases.reserve(nHarmonics * 2);
	
	_spsModelAlgo->input("magnitudes").set(_sineMags);
	_spsModelAlgo->input("frequencies").set(_sineFreqs);
	_spsModelAlgo->input("phases").set(_sinePhases);
	_spsModelAlgo->input("stocenv").set(_stocEnv);
	
	_spsModelAlgo->output("frame").set(outputFrame);
	_spsModelAlgo->output("sineframe").set(sineOutputFrame);
	_spsModelAlgo->output("stocframe").set(stocOutputFrame);

	owner.validateSubclassOwner();
}
void SynthesisContainer::sinePlusStochastic_svc::doSetPitch(float new_f0){
	_voice_f0 = new_f0;
}
void SynthesisContainer::sinePlusStochastic_svc::doSetLoudness(float new_loudness){
	_voice_loudness = new_loudness;
}
void SynthesisContainer::sinePlusStochastic_svc::doSetTimbre(baseTimbre &&new_timbre){
	sinePlusStochasticTimbre *spstPtr = static_cast<sinePlusStochasticTimbre*>(&new_timbre);
	const int minIdx = baseTimbre::pMapIdx::min_e;
	const int maxIdx = baseTimbre::pMapIdx::max_e;
	
	if(spstPtr){
		spstPtr->tranposeFreqs = std::min(spstPtr->tranposeFreqs,
			sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::transpose_e)[maxIdx]);
		spstPtr->tranposeFreqs = std::max(spstPtr->tranposeFreqs,
			 sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::transpose_e)[minIdx]);
		
		spstPtr->shiftFreqs = std::min(spstPtr->shiftFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::shift_e)[maxIdx]);
		spstPtr->shiftFreqs = std::max(spstPtr->shiftFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::shift_e)[minIdx]);
		
		spstPtr->tiltFreqs = std::min(spstPtr->tiltFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::tilt_e)[maxIdx]);
		spstPtr->tiltFreqs = std::max(spstPtr->tiltFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::tilt_e)[minIdx]);
		
		spstPtr->stretchFreqs = std::min(spstPtr->stretchFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::stretch_e)[maxIdx]);
		spstPtr->stretchFreqs = std::max(spstPtr->stretchFreqs,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::stretch_e)[minIdx]);
		
		spstPtr->tonalStochRatio = std::min(spstPtr->tonalStochRatio,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::tonalStochRatio_e)[maxIdx]);
		spstPtr->tonalStochRatio = std::max(spstPtr->tonalStochRatio,
				sinePlusStochasticTimbre::pMap.at(sinePlusStochasticTimbre::spsTimbreEnum::tonalStochRatio_e)[minIdx]);
		
		
		baseTimbre* tempBase = _voice_timbre.get();
		sinePlusStochasticTimbre* tempDerived =  dynamic_cast<sinePlusStochasticTimbre*>(tempBase);
		if (tempDerived){
			tempDerived->tranposeFreqs = spstPtr->tranposeFreqs;
			tempDerived->shiftFreqs = spstPtr->shiftFreqs;
			tempDerived->tiltFreqs = spstPtr->tiltFreqs;
			tempDerived->tonalStochRatio = spstPtr->tonalStochRatio;
		}
	}
}
void SynthesisContainer::sinePlusStochastic_svc::doResetPhases(bool shouldReset){
	shouldResetPhasesFromData = shouldReset;
}

void SynthesisContainer::sinePlusStochastic_svc::doProduceBlock(float* block){}

void SynthesisContainer::sinePlusStochastic_svc::prepareRequestedSynthesis(){
	if (!(_spsModelAlgo)){
		std::cerr << "algos nullptr\n";
		return;
	}
	if (!(owner.sc_state.isInitialized)){
		std::cerr << "unititialized \n";
		return;
	}
	
	_state.frameReady = false;
	
	const size_t nHarm = static_cast<size_t>(owner.sc_analysisData->numHarmonics);

	const float interp = owner.sc_frameInterp;
	//if (owner.getIsCaughtUpToCurrentFrame()){
		// then sc_frameInterp = 0.f
		// jassert ( interp == 0.f );
	
	// set next freqs, mags, phases, stoch to owner.sc_nextTargetFrameIdx
	
	const realVec_t &freqs = owner.sc_analysisData->spsModelFreqs[chanNum][owner.sc_nextTargetFrameIdx];
	const realVec_t &mags = owner.sc_analysisData->spsModelMags[chanNum][owner.sc_nextTargetFrameIdx];
	const realVec_t &phases = owner.sc_analysisData->spsModelPhases[chanNum][owner.sc_nextTargetFrameIdx];
	
	const realVec_t &stocenvs = owner.sc_analysisData->spsModelStocEnvs[chanNum][owner.sc_nextTargetFrameIdx];

	size_t last = std::min(nHarm, freqs.size());
	_sineFreqs.assign(freqs.begin(), freqs.begin() + last);
	_sineMags.assign(mags.begin(), mags.begin() + last);
	_sinePhases.assign(phases.begin(), phases.begin() + last);
	_stocEnv.assign(stocenvs.begin(), stocenvs.begin() + last);
	//} else {
		// not caught up, so interpolating
		
		// set next vectors to owner.sc_nextTargetFrameIdx, previous vecors to owner.sc_prevTargetFrameIdx
	//}
	
	if (!shouldResetPhasesFromData){
		_sinePhases.resize(0);
	}
	
	const sinePlusStochasticTimbre *spst = dynamic_cast<sinePlusStochasticTimbre*>(_voice_timbre.get());
	if(spst){
		// set timbral params
		const float ofst = spst->shiftFreqs;
		std::transform(_sineFreqs.begin(), _sineFreqs.end(), _sineFreqs.begin(), [&ofst](auto& c){
			return c+ofst;
		});
		
		const float semitones = spst->tranposeFreqs;
		const float ratio = tsaraCommon::semitones2ratio(semitones) * (_voice_f0 / concertPitch);
		

		const float tilt = spst->tiltFreqs;
		std::transform(_sineMags.begin(), _sineMags.end(), _sineFreqs.begin(), _sineMags.begin(), [tilt](float& mag, float& freq){
			freq = std::min(freq, 20000.f);
			float a0 = (freq / 20000.f);
			float a1 = 1.f - (freq / 20000.f);
			float amp = std::pow(10.f, (mag / 20.f));
			amp *= ((tilt * a0) + ((1.f - tilt) * a1));
			float val = 20.f * std::log10(amp);
			return val;
		});
		
		std::transform(_sineFreqs.begin(), _sineFreqs.end(), _sineFreqs.begin(), [&ratio](auto& c){
			auto val =  c*ratio;
			return val;
		});
		
		// should really remove_if too high magnitude
	}
	// compute algo
	_spsModelAlgo->compute();
	
	_state.frameReady = true;
	shouldResetPhasesFromData = false;
}
void SynthesisContainer::sinePlusStochastic_svc::doPrepareSynthesis(){
	CHECK_INITIALIZED;
	prepareRequestedSynthesis();
}
float SynthesisContainer::sinePlusStochastic_svc::doProduceSample(){
	if (! _state.frameReady){
		std::cerr << "FRAME NOT READY\n";
		return 0.f;
	}
	unsigned int idxWithinFrame =	owner.totalFramesSamplewiseIdx % owner.sc_analysisData->hopSize;
	float val = outputFrame[idxWithinFrame];
	return val;
}


#endif /* sinePlusStochastic_svc_impl_h */
