//
//  distance.cpp
//  probability_graph_testing
//
//  Created by Nicholas Solem on 5/15/23.
//  Copyright © 2023 Nicholas Solem. All rights reserved.
//

#include <stdio.h>
#include "Distance.h"
#include "graph.h"

//namespace nvs {
//namespace dst {

std::optional<float> nvs::dst::euc_distance(std::vector<float> const &x0, std::vector<float> const &x1){
	if (x0.size() > x1.size()){
		std::cerr << "euc_distance exit: erroneous size mismatch\n";
		return {};
	}
	std::vector<float> diff(x0.size());
	float d = 0.f;
	auto op = [](const float& a, const float& b){
		float dif = (a - b);
		return dif * dif;
	};
	std::transform(x0.begin(), x0.end(), x1.begin(), diff.begin(), op);
	d = std::accumulate(diff.begin(), diff.end(), 0.f);
	d = std::sqrtf(d);
	return d;
}

std::optional<float> nvs::dst::kullback_leibler_divergence(std::vector<float> const &x0, std::vector<float> const &x1){
	if (x0.size() != x1.size()){
		std::cerr << "kl divergence exit: size mismatch\n";
		return {};
	}
	std::vector<float> logRatio(x0.size());
	float d = 0.f;
	auto op = [](const float& p, const float& q){
		float q_lim = std::max(q, 0.0000001f);
		float ratio = (p / q_lim);
		return p * std::log(ratio);
	};
	std::transform(x0.begin(), x0.end(), x1.begin(), logRatio.begin(), op);
	d = std::accumulate(logRatio.begin(), logRatio.end(), 0.f);
	return d;
}

float nvs::dst::loudnessDistance(float src_loud, float dst_loud){
	float dif = src_loud - dst_loud;
	dif *= dif;
	return std::sqrtf(dif);
}

float nvs::dst::pitchDistanceRatio(float src_freq, float dst_freq){
	src_freq = std::abs(src_freq);
	dst_freq = std::abs(dst_freq);
	src_freq = std::max(src_freq, 0.0000001f);
	dst_freq = std::max(dst_freq, 0.0000001f);
	float val = src_freq / dst_freq;
	val = std::log(val);
	val = std::abs(val);
	return val;
}

bool nvs::dst::pitchDistanceInRange(float src_freq, float dst_freq, nvs::sgt::connectionHeuristics const& settings){
	constexpr float cutoff = 500.f;
	src_freq = std::abs(src_freq);
	dst_freq = std::abs(dst_freq);
	float val;
	if ((src_freq >= cutoff) & (dst_freq >= cutoff)){ // ratio based
		if (dst_freq > src_freq){	// swap
			float tmp;
			tmp = dst_freq;
			dst_freq = src_freq;
			src_freq = tmp;
		}
		dst_freq = std::max(dst_freq, 0.0000001f);
		val = src_freq / dst_freq;
		return (val < settings.maxPitchDeviationRatio);
	} else {	// hertz difference based
		val = std::abs(src_freq - dst_freq);
		return (val < settings.maxPitchDeviationAbsolute);
	}
}

//}	// namespace dst
//}	// namespace nvs
