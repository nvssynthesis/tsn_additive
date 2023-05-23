//
//  sineModel_svc_impl.h
//  tsara_additive
//
//  Created by Nicholas Solem on 5/22/23.
//  Copyright © 2023 nvssynthesis. All rights reserved.
//

#ifndef sineModel_svc_impl_h
#define sineModel_svc_impl_h
#include "Synthesis.h"

SynthesisContainer::sineModel_svc::sineModel_svc(const SynthesisContainer &owningSynth, unsigned int i)
:   									abstractSynthVoiceChannel(owningSynth, i)
{
	_voice_timbre = sineModelTimbre::factory();
	
	owner.validateSubclassOwner();
	_state.isInitialized = true;

}
void SynthesisContainer::sineModel_svc::doInit() {
	_sineModelAlgo = owner.ess->factory.create("SineModelSynth",
											   "fftSize", owner.sc_analysisData->frameSize,
												   "hopSize", owner.sc_analysisData->hopSize,
												   "sampleRate", owner.sc_analysisData->sampleRate);
	_ifftAlgo = owner.ess->factory.create("IFFT", "normalize", true, "size", owner.sc_analysisData->frameSize);
	_overlapAddAlgo = owner.ess->factory.create("OverlapAdd", "frameSize", owner.sc_analysisData->frameSize, "hopSize", owner.sc_analysisData->hopSize);
	
	unsigned int nHarmonics = owner.sc_analysisData->numHarmonics;
	_sineMags.reserve(nHarmonics * 2);
	_sineFreqs.reserve(nHarmonics * 2);
	_sinePhases.reserve(nHarmonics * 2);
	
	size_t fftSize = (owner.sc_analysisData->frameSize / 2) + 1;
	_complexSpec.resize(fftSize);
	
	_sineModelAlgo->input("magnitudes").set(_sineMags);
	_sineModelAlgo->input("frequencies").set(_sineFreqs);
	_sineModelAlgo->input("phases").set(_sinePhases);
	_sineModelAlgo->output("fft").set(_complexSpec);
	
	_ifftAlgo->input("fft").set(_complexSpec);
	_ifftAlgo->output("frame").set(_ifftFrame);

	_overlapAddAlgo->input("signal").set(_ifftFrame);
	_overlapAddAlgo->output("signal").set(_overlapAddFrame);
	
	owner.validateSubclassOwner();
}

void SynthesisContainer::sineModel_svc::prepareRequestedSynthesis(){
	if (!(_sineModelAlgo && _ifftAlgo && _overlapAddAlgo)){
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
	if (owner.getIsCaughtUpToCurrentFrame()){
		// then sc_frameInterp = 0.f
		jassert ( interp == 0.f );
		const realVec_t &freqs = owner.sc_analysisData->sineModelFreqs[chanNum][owner.sc_nextTargetFrameIdx];
		const realVec_t &mags = owner.sc_analysisData->sineModelMags[chanNum][owner.sc_nextTargetFrameIdx];
		const realVec_t &phases = owner.sc_analysisData->sineModelPhases[chanNum][owner.sc_nextTargetFrameIdx];
		size_t last = std::min(nHarm, freqs.size());
		_sineFreqs.assign(freqs.begin(), freqs.begin() + last);
		_sineMags.assign(mags.begin(), mags.begin() + last);
		_sinePhases.assign(phases.begin(), phases.begin() + last);
	} else {
		// not caught up, so interpolating
		const realVec_t &freqsPrev = owner.sc_analysisData->sineModelFreqs[chanNum][owner.sc_prevTargetFrameIdx];
		const realVec_t &magsPrev = owner.sc_analysisData->sineModelMags[chanNum][owner.sc_prevTargetFrameIdx];
		const realVec_t &phasesPrev = owner.sc_analysisData->sineModelPhases[chanNum][owner.sc_prevTargetFrameIdx];
		const realVec_t &freqsNext = owner.sc_analysisData->sineModelFreqs[chanNum][owner.sc_nextTargetFrameIdx];
		const realVec_t &magsNext = owner.sc_analysisData->sineModelMags[chanNum][owner.sc_nextTargetFrameIdx];
		const realVec_t &phasesNext = owner.sc_analysisData->sineModelPhases[chanNum][owner.sc_nextTargetFrameIdx];
		
		size_t last = std::min(nHarm, freqsPrev.size());
		_sineFreqs.assign(freqsPrev.begin(), freqsPrev.begin() + last);
		_sineMags.assign(magsPrev.begin(), magsPrev.begin() + last);
		_sinePhases.assign(phasesPrev.begin(), phasesPrev.begin() + last);
		std::cout << interp << '\n';
		auto splitPoint = std::transform(_sineMags.begin(), _sineMags.begin() + last, _sineMags.begin(),
		[interp](float mag){
			// at interp = 0, we want no change.
			// at interp = 1, we want full attenuation.
			float amp = std::pow(10.f, (mag / 20.f));
			amp *= std::sqrt(1.f - interp);
			return 20.f * std::log10(amp);
		});

		last = std::min(nHarm, freqsNext.size());
		std::copy(freqsNext.begin(), freqsNext.begin() + last, std::back_inserter(_sineFreqs));
		std::copy(magsNext.begin(), magsNext.begin() + last, std::back_inserter(_sineMags));
		std::copy(phasesNext.begin(), phasesNext.begin() + last, std::back_inserter(_sinePhases));
		std::transform(splitPoint, splitPoint + last, splitPoint,
		[interp](float mag){
			// at interp = 0, we want full attenuation.
			// at interp = 1, we want no change.
			float amp = std::pow(10.f, (mag / 20.f));
			amp *= std::sqrt(interp);
			return 20.f * std::log10(amp);
		});
	}
	if (!shouldResetPhasesFromData){
		_sinePhases.resize(0);
	}
	
	const sineModelTimbre *smt = dynamic_cast<sineModelTimbre*>(_voice_timbre.get());
	if(smt){
		const float ofst = smt->shiftFreqs;
		std::transform(_sineFreqs.begin(), _sineFreqs.end(), _sineFreqs.begin(), [&ofst](auto& c){
			return c+ofst;
		});
		
		const float semitones = smt->tranposeFreqs;
		const float ratio = tsaraCommon::semitones2ratio(semitones) * (_voice_f0 / concertPitch);
		
		std::transform(_sineFreqs.begin(), _sineFreqs.end(), _sineFreqs.begin(), [&ratio](auto& c){
			auto val =  c*ratio;
			return val;
		});
		const float tilt = smt->tiltFreqs;
		std::transform(_sineMags.begin(), _sineMags.end(), _sineFreqs.begin(), _sineMags.begin(), [tilt](float& mag, float& freq){
			float a0 = (freq / 20000.f);
			float a1 = 1.f - (freq / 20000.f);
			float amp = std::pow(10.f, (mag / 20.f));
			amp *= ((tilt * a0) + ((1.f - tilt) * a1));
			float val = 20.f * std::log10(amp);
			if (val != val)
				val = -120.f;
			return val;
		});
	}

	_sineModelAlgo->compute();				// processes _sinePhases, outputs _complexSpec
	_ifftAlgo->compute();					// takes in _complexSpecTotal, outputs reals
	_overlapAddAlgo->compute();

	auto blocksPerHop = owner.getBlocksInAHop();
	unsigned int blockLimit = (unsigned int)ceil(blocksPerHop - 1);
	if (owner.sc_blockIdx == blockLimit){
		outputFrames[outputFrameToOverwrite] = _overlapAddFrame;
		++outputFrameToOverwrite;
		outputFrameToOverwrite %= getNumFramesToOverwritePerBlock();
	}
	_state.frameReady = true;
	shouldResetPhasesFromData = false;
}
// prep of 1st frame synthesis-should auto-set to 'prepareNext' after single call
void SynthesisContainer::sineModel_svc::doPrepareSynthesis() {
	CHECK_INITIALIZED;
	prepareRequestedSynthesis();
}
float SynthesisContainer::sineModel_svc::doProduceSample() {
	if (! _state.frameReady){
		std::cerr << "FRAME NOT READY\n";
		return 0.f;
	}
	unsigned int whichFrame = floor(owner.totalFramesSamplewiseIdx / owner.sc_analysisData->hopSize);
	unsigned int idxWithinFrame =	owner.totalFramesSamplewiseIdx % owner.sc_analysisData->hopSize;
	float val = outputFrames[whichFrame][idxWithinFrame];
	return val;
}

void SynthesisContainer::sineModel_svc::doProduceBlock(float* block){
	CHECK_INITIALIZED;
	_state.frameReady = false;

	float nOverlap = owner.getNumOverlap();
	int blockSize = owner.sc_blockSize;
	float hopsPerBlock = owner.getHopsPerBlock();

	std::fill_n(&(block[0]), blockSize, 0.f);

	for (int i = 0; i < blockSize; ++i){
		block[i] = sin(2.0 * M_PI * (double)i * 440.0 * owner.sc_playbackSampleRateInv);
	}
#if 1
	if (owner.sc_analysisData->frameSize < blockSize){
		// this var count to nOverlap - 1.
		int whichFrame = 0;
		for (int i = 0; i < hopsPerBlock; ++i){
			prepareSynthesis();
			
			if (i == 1)
				outputFrames[0] = _overlapAddFrame;
			else if (i == 3)
				outputFrames[1] = _overlapAddFrame;

			if (whichFrame == (nOverlap - 1)){
				int stIdx = static_cast<int> (((float)blockSize / hopsPerBlock) * (float)0);
				[[maybe_unused]]
				int endIdx = stIdx + owner.sc_analysisData->frameSize - 1;
				
				std::copy(&(outputFrames[0][0]),
						  &(outputFrames[0][owner.sc_analysisData->frameSize - 1]),
						  &(block[stIdx]));
				
				stIdx = static_cast<int> (((float)blockSize / hopsPerBlock) * (float)2);
				endIdx = stIdx + owner.sc_analysisData->frameSize - 1;
				std::copy(&(outputFrames[1][0]),
						  &(outputFrames[1][owner.sc_analysisData->frameSize - 1]),
						  &(block[stIdx]));
			}
			whichFrame += 1;
			whichFrame %= (int)ceil(nOverlap);
		}
	} else if (owner.sc_analysisData->frameSize >= blockSize){
#if 0
		if (owner.sc_blockIdx == 0){
			for (int i = 0; i < (int)nOverlap; ++i){
				prepareSynthesis();
				outputFrames[i] = _overlapAddFrame;
			}
		}
		int frameStIdx = blockSize * owner.sc_blockIdx;
		int frameEndIdx = frameStIdx + blockSize - 1;
		std::copy(&(outputFrames[nOverlap - 1][frameStIdx]),
				  &(outputFrames[nOverlap - 1][frameEndIdx]),
				  &(block[0]));
#endif
	}
#endif
}
void SynthesisContainer::sineModel_svc::doSetPitch(float new_f0) {
	_voice_f0 = new_f0;
}
void SynthesisContainer::sineModel_svc::doSetLoudness(float new_loudness) {
	_voice_loudness = new_loudness;
}
void SynthesisContainer::sineModel_svc::doSetTimbre(baseTimbre &&new_timbre) {
	sineModelTimbre *smtPtr = static_cast<sineModelTimbre*>(&new_timbre);
	const int minIdx = baseTimbre::pMapIdx::min_e;
	const int maxIdx = baseTimbre::pMapIdx::max_e;
	
	if(smtPtr){
		smtPtr->tranposeFreqs = std::min(smtPtr->tranposeFreqs,
			 sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::transpose_e)[maxIdx]);
		smtPtr->tranposeFreqs = std::max(smtPtr->tranposeFreqs,
			 sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::transpose_e)[minIdx]);
		
		smtPtr->shiftFreqs = std::min(smtPtr->shiftFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::shift_e)[maxIdx]);
		smtPtr->shiftFreqs = std::max(smtPtr->shiftFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::shift_e)[minIdx]);
		
		smtPtr->tiltFreqs = std::min(smtPtr->tiltFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::tilt_e)[maxIdx]);
		smtPtr->tiltFreqs = std::max(smtPtr->tiltFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::tilt_e)[minIdx]);
		
		smtPtr->stretchFreqs = std::min(smtPtr->stretchFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::stretch_e)[maxIdx]);
		smtPtr->stretchFreqs = std::max(smtPtr->stretchFreqs,
				sineModelTimbre::pMap.at(sineModelTimbre::smTimbreEnum::stretch_e)[minIdx]);
		
		
		baseTimbre* tempBase = _voice_timbre.get();
		sineModelTimbre* tempDerived =  dynamic_cast<sineModelTimbre*>(tempBase);
		if (tempDerived){
			tempDerived->tranposeFreqs = smtPtr->tranposeFreqs;
			tempDerived->shiftFreqs = smtPtr->shiftFreqs;
			tempDerived->tiltFreqs = smtPtr->tiltFreqs;
		}
	}
}
void SynthesisContainer::sineModel_svc::doResetPhases(bool shouldReset){
	shouldResetPhasesFromData = shouldReset;
}
#if 0
void SynthesisContainer::sineModel_svc::getPhasesFirst(){
	const realVec_t &phases = owner.sc_analysisData->sineModelPhases[chanNum][owner.sc_nextTargetFrame];
	_sinePhases.assign(phases.begin(), phases.end());
	getPhases = &sineModel_svc::getPhasesNext;
}
void SynthesisContainer::sineModel_svc::getPhasesNext(){
	_sinePhases.resize(0);
}
#endif



#endif /* sineModel_svc_impl_h */
