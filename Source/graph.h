//
//  graph.h
//  pool2graph
//
//  Created by Nicholas Solem on 4/3/23.
//  Copyright © 2023 Nicholas Solem. All rights reserved.
//

#ifndef graph_h
#define graph_h
#include <cstring>
#include <iostream>
#include <fstream>
#include <random>
#include <array>
#include <limits>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/core/span.hpp>
#include "Analysis.h"

namespace nvs{
typedef std::vector<std::vector<float>> pca_matrix_t;
 class sound_representation {
private:
	pca_matrix_t PCAmat;
	std::vector<float> fundamental_freq_vec;
	std::vector<float> voicedness_vec;
	std::vector<float> loudness_vec;
public:
	bool verifyEqualLengths(){
		size_t size = PCAmat.size();
		bool equalLengths = (size == fundamental_freq_vec.size()) && (size = loudness_vec.size());
		return equalLengths;
	}
	pca_matrix_t getPCAmat() const {return PCAmat;}
	std::vector<float> getFundamentalFreqVec() const {return fundamental_freq_vec;}
	std::vector<float> getVoicednessVec() const {return voicedness_vec;}
	std::vector<float> getLoudnessVec() const {return loudness_vec;}
	void setPCAmat(const pca_matrix_t &other){
		PCAmat = other;
	}
	void setFundamentalFreqVec(std::vector<float> const &other){
		fundamental_freq_vec = other;
	}
	void setVoicednessVec(std::vector<float> const &other){
		voicedness_vec = other;
	}
	void setLoudnessVec(std::vector<float> const &other){
		loudness_vec = other;
	}
};
struct vertex_property{
	double timestamp;
	unsigned int cluster;
	float fundamentalFreq;
	float voicedness;
	float loudness;
};
struct edge_property{
    float timbral_distance;
    float temporal_distance;
	float pitch_distance;
	float loudness_distance;
};
//typedef boost::property < boost::edge_weight_t, double > weight_prop_t;
typedef boost::adjacency_list < boost::vecS, 		// out-edges, a Sequence or an AssociativeContainer
						boost::vecS, 		// vertices, a Sequence or a RandomAccessContainer
						boost::directedS,
						vertex_property,	// VertexProperty
						edge_property	 	// EdgeProperty
						> DirectedGraph_t;
typedef boost::graph_traits<DirectedGraph_t>::edge_iterator edge_iterator_t;
typedef DirectedGraph_t::adjacency_iterator adjacency_iter_t;
typedef DirectedGraph_t::vertex_descriptor vertex_descriptor_t;
typedef DirectedGraph_t::vertex_iterator vertex_iterator_t;
typedef DirectedGraph_t::edge_descriptor edge_descriptor_t;

typedef boost::random::mt19937 random_gen_t;
typedef boost::random::discrete_distribution<size_t,double> discrete_distr_t;
typedef boost::normal_distribution<double> normal_distr_t;

static random_gen_t gen;

static std::vector<float> processBarks(std::vector<float> const &x, size_t N);

static sound_representation getSoundRepresentationFromPool(essentia::Pool const &p){
	nvs::tsaraCommon common;
	
	const std::string globals_stereo_mode_str = tsaraCommon::getGlobalsString(tsaraCommon::globals_e::stereo_mode);
	const std::string mono_str = p.getSingleStringPool().at(std::move(globals_stereo_mode_str));
	tsaraCommon::stereoModes_e stMode_e = tsaraCommon::chanDescr2stereoModeEnumMap.at(std::move(mono_str));

	std::string const barkStr = tsaraCommon::getLowlevelString(tsaraCommon::algo_e::bark, 0, stMode_e);

	auto vReal = p.getVectorRealPool();
	
	std::vector<std::vector<float>> barkVec = vReal.at(barkStr);
	essentia::Pool tmp_p;
	for (auto &barkSpec : barkVec){
		barkSpec = processBarks(barkSpec, 23UL);	// takes 4th root for loudness of bark bands
	
		tmp_p.add(barkStr, barkSpec);
	}
	
	std::string const nsOut {"pca.bark"};
	essentia::standard::Algorithm* pca = common.getFactory().create("PCA",
							"dimensions", 11,
							"namespaceIn", std::move(barkStr),
							"namespaceOut", std::move(nsOut));
	essentia::Pool pca_pool;
	pca->input("poolIn").set(tmp_p);
	pca->output("poolOut").set(pca_pool);
	pca->compute();
	// load essentia pool into data struct
	// std::unique_ptr<nvs::AnalysisData> a_data = std::make_unique<nvs::AnalysisData>(p);	// do we even need this once we have PCA below?
	const pca_matrix_t pca_mat = pca_pool.getVectorRealPool().at(nsOut);
	
	const auto freqs = p.getRealPool().at(tsaraCommon::getLowlevelString(tsaraCommon::algo_e::pitchYinProbabilistic_freqs, 0,stMode_e));
	const auto voicedness = p.getRealPool().at(tsaraCommon::getLowlevelString(tsaraCommon::algo_e::pitchYinProbabilistic_voicedness ,0,stMode_e));
	const auto loudness = p.getRealPool().at(tsaraCommon::getLowlevelString(tsaraCommon::algo_e::loudness,0 ,stMode_e));

	sound_representation rep;
	rep.setPCAmat(pca_mat);
	rep.setFundamentalFreqVec(freqs);
	rep.setVoicednessVec(voicedness);
	rep.setLoudnessVec(loudness);
	return rep;
}
// start essentia. this can be in small scope, we just need to fill a_data.
static sound_representation getSoundRepresentationFromFile(std::string filename  = "/Users/nicholassolem/development/audio for analysis/mr sandman tiny.yaml"){
	nvs::tsaraCommon common;
	essentia::Pool p = common.file2pool(filename, "yaml");
	
	return getSoundRepresentationFromPool(p);
}

static std::optional<float> euc_distance(std::vector<float> const &x0, std::vector<float> const &x1){
	if (x0.size() != x1.size()){
		std::cerr << "euc_distance exit: size mismatch\n";
		return {};
	}
	std::vector<float> diff(x0.size());
	float d = 0.f;
	auto op = [](const float& a, const float& b){
		float dif = (a - b);
		return dif * dif;
	};
	std::transform(std::begin(x0), std::end(x0), std::begin(x1), std::begin(diff), op);
	d = std::accumulate(std::begin(diff), std::end(diff), 0.f);
	d = std::sqrtf(d);
	return d;
}
static float kullback_leibler_divergence(std::vector<float> const &x0, std::vector<float> const &x1){
	if (x0.size() != x1.size()){
		std::cerr << "kl divergence exit: size mismatch\n";
		return -1.f;
	}
	std::vector<float> logRatio(x0.size());
	float d = 0.f;
	auto op = [](const float& p, const float& q){
		float q_lim = std::max(q, 0.0000001f);
		float ratio = (p / q_lim);
		return p * std::log(ratio);
	};
	std::transform(std::begin(x0), std::end(x0), std::begin(x1), std::begin(logRatio), op);
	d = std::accumulate(std::begin(logRatio), std::end(logRatio), 0.f);
	return d;
}
static std::vector<float> splitBarks(std::vector<float> const &x, size_t N){
	N = std::min(N, x.size());
	return std::vector<float>(x.begin(), x.begin() + N);
}
static std::vector<float> barkLoudness(std::vector<float> const &x){
	std::vector<float> y(x.size());
	std::transform(std::begin(x), std::end(x), std::begin(y), [](const float& a){
// "In each band we estimate a loudness contribution as the fourth root of the power on the band; this is close to a loudness measure sug- gested in (Rossing, Moore, and Wheeler 2002)". (puckette)
		return std::sqrtf(std::sqrtf(a));
	});
	return y;
}
static std::vector<float> processBarks(std::vector<float> const &x, size_t N){
	std::vector<float> v = splitBarks(x, N);
	v = barkLoudness(v);
	return v;
}

static std::vector<size_t> getRandomIndices(size_t maxValInclusive, size_t numVals = 100){
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<size_t> distr(0, maxValInclusive); // define the range
	
	std::vector<size_t> randVec (numVals, 0);
	for(size_t n = 0; n < numVals; ++n){
		size_t val = distr(gen);
		randVec[n] = val;
	}
	return randVec;
}

static float getMaxTimbralDistanceForConnection(pca_matrix_t const &pca_mat, float percentile = 0.15f, size_t nSearch = 100UL)	{
	const size_t numFrames = pca_mat.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	for (int src = 0; src < nSearch; ++src){
		auto x0 = pca_mat[randIdx[src]];
		for (int dst = 0; dst < nSearch; ++dst){
			if (src == dst)
				continue;
			float d = euc_distance(x0, pca_mat[randIdx[dst]]).value_or(-1.f);
			distances.push_back(d);
		}
	}
	std::sort(distances.begin(), distances.end());
	size_t pctlIdx = size_t((distances.size()-1) * percentile);
	float maxDistanceForConnection = distances[pctlIdx];
	std::cout << "max pca DistanceForConnection: " << maxDistanceForConnection << '\n';
	
	nth_element( distances.begin(), distances.begin() + pctlIdx, distances.end());
	std::cout << "					nth_element: " << *(distances.begin() + pctlIdx) << '\n';

	return maxDistanceForConnection;
}
static float getMaxTimbralDistanceForConnectionFromVertex(const vertex_descriptor_t v, pca_matrix_t const &pca_mat, float percentile = 0.15f, size_t nSearch = 100UL){
	auto src_timbre = pca_mat[v];
	const size_t numFrames = pca_mat.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	distances.reserve(nSearch);
	for (auto dst = 0; dst < nSearch; ++dst){
		if (dst == v)
			continue;
		float d = euc_distance(src_timbre, pca_mat[randIdx[dst]]).value_or(-1.f);
		distances.push_back(d);
	}
	size_t pctlIdx = size_t((distances.size()-1) * percentile);
	nth_element( distances.begin(), distances.begin() + pctlIdx, distances.end());
	
	return *(distances.begin() + pctlIdx);
}

static float loudnessDistance(float src_loud, float dst_loud);
static float getMaxLoudnessDistanceForConnection(std::vector<float> const &loudnessVec, float percentile = 0.15f, size_t nSearch = 100UL){
	const size_t numFrames = loudnessVec.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	for (int src = 0; src < nSearch; ++src){
		auto x0 = loudnessVec[randIdx[src]];
		for (int dst = 0; dst < nSearch; ++dst){
			if (src == dst)
				continue;
			float d = loudnessDistance(x0, loudnessVec[randIdx[dst]]);
			distances.push_back(d);
		}
	}
	std::sort(distances.begin(), distances.end());
	size_t pctlIdx = size_t((distances.size()-1) * percentile);
	float maxDistanceForConnection = distances[pctlIdx];
	std::cout << "max loudness DistanceForConnection: " << maxDistanceForConnection << '\n';
	return maxDistanceForConnection;
}


static float loudnessDistance(float src_loud, float dst_loud){
	float dif = src_loud - dst_loud;
	dif *= dif;
	return std::sqrtf(dif);
}

struct connectionHeuristics{
	float timbralPercentile {0.12f};
	size_t timbralNumSearch {100UL};
	float loudnessPercentile {0.15f};
	size_t loudnessNumSearch {100UL};
	float maxPitchDeviationRatio {5.f / 4.f};	// over 500 Hz
	float maxPitchDeviationAbsolute {100.f};	// under 500 Hz
	
	size_t minimumConnections {6};
	size_t minimumConnectionCandidates {18};
};
static float pitchDistanceRatio(float src_freq, float dst_freq){
	src_freq = std::abs(src_freq);
	dst_freq = std::abs(dst_freq);
	src_freq = std::max(src_freq, 0.0000001f);
	dst_freq = std::max(dst_freq, 0.0000001f);
	float val = src_freq / dst_freq;
	val = std::log(val);
	val = std::abs(val);
	return val;
}
static bool pitchDistanceInRange(float src_freq, float dst_freq, connectionHeuristics const& settings){
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
static void reduceEdgesBasedOnDegree(DirectedGraph_t &dg){
	// vertex iterator
	
	// get out degree of vertex
}
// create directed graph based on numFrames
static DirectedGraph_t createGraphFromSoundRepresentation(sound_representation const &sound_rep, connectionHeuristics const settings = connectionHeuristics()){
	const pca_matrix_t &pca_mat = sound_rep.getPCAmat();
	const auto &freqs = sound_rep.getFundamentalFreqVec();
	const auto &voicedness = sound_rep.getVoicednessVec();
	const auto &loudness = sound_rep.getLoudnessVec();
	const size_t numFrames = pca_mat.size();
	DirectedGraph_t dg(numFrames);

	const size_t minimumConnections = settings.minimumConnections;
	const size_t minimumConnectionCandidates = settings.minimumConnectionCandidates;

	const float maxLoudnessDev = getMaxLoudnessDistanceForConnection(loudness, settings.loudnessPercentile, settings.loudnessNumSearch);
//========================================================================================
	// this stuff is used in case we do not attain minimumConnections
	typedef std::pair<size_t, edge_property> idxAndEdgeProp;
	auto idxAndEdgeProp_CmpByTimbre = [](idxAndEdgeProp a, idxAndEdgeProp b){
		return a.second.timbral_distance < b.second.timbral_distance;
	};
	auto idxAndEdgeProp_CmpByPitch = [](idxAndEdgeProp a, idxAndEdgeProp b){
		return a.second.pitch_distance < b.second.pitch_distance;
	};
	edge_property maximallyDistantEdgeProp = {
		.timbral_distance = std::numeric_limits<float>::max(),
		.temporal_distance = std::numeric_limits<float>::max(),
		.pitch_distance = std::numeric_limits<float>::max(),
		.loudness_distance = std::numeric_limits<float>::max()
	};
//========================================================================================

	for (int src = 0; src < numFrames; ++src){
		dg[src].timestamp = src;
		dg[src].cluster = 0;
		dg[src].fundamentalFreq = freqs[src];
		dg[src].voicedness = voicedness[src];
		dg[src].loudness = loudness[src];

		float maxTimbralDistanceForConnection = getMaxTimbralDistanceForConnectionFromVertex(src, pca_mat, settings.timbralPercentile, settings.timbralNumSearch);
		
		auto x0 = pca_mat[src];
		size_t num_connections_from_node = 0;

		// prepare array of minima in case nothing no connections are made
		// possible to use std::nth_element??
		// yes it is possible, but my way may be better because the list of distances is
		// 	only ever built up if the distance is a minimum. if i used nth_element, i
		//  would need the whole list (or more than the minumum few) just to throw away what's above the nth.
		std::vector<idxAndEdgeProp> minimaIdxAndVals(minimumConnectionCandidates,
											{ (size_t) 0, maximallyDistantEdgeProp } );
		for (int dst = 0; dst < numFrames; ++dst){
			if (src == dst)
				continue;
			auto x1 = pca_mat[dst];
			edge_property ep;
			ep.temporal_distance = (float)std::abs(dst - src);
			ep.timbral_distance = euc_distance(x0, x1).value_or(-1.f);
			ep.pitch_distance = pitchDistanceRatio(dg[src].fundamentalFreq, freqs[dst]);
			ep.loudness_distance = loudnessDistance(dg[src].loudness, loudness[dst]);
			
			auto maxOfMinVec = std::max_element(minimaIdxAndVals.begin(), minimaIdxAndVals.end(), idxAndEdgeProp_CmpByTimbre);
			if ((ep.timbral_distance < maxTimbralDistanceForConnection)
				& (pitchDistanceInRange(dg[src].fundamentalFreq, freqs[dst], settings))
				& (ep.loudness_distance < maxLoudnessDev)
				)
			{
				boost::add_edge(src, dst, ep, dg);
				++num_connections_from_node;
			}
			else if (ep.timbral_distance < maxOfMinVec->second.timbral_distance){
				// only add edge property to minimal list if it's not connected
				maxOfMinVec->first = dst;
				maxOfMinVec->second = ep;
			}
		}
#if 0
		if (num_connections_from_node < minimumConnections){
			// fill 'minima to be used' with 'minimaIdxAndVals' whose indices correspond with smallest pitch differences.
			// these will be the minima that are used.
			std::sort(minimaIdxAndVals.begin(), minimaIdxAndVals.end(), idxAndEdgeProp_CmpByPitch);
			// connect the smallest pitch differences from 'minima to be used'
			for (size_t i = 0; i < (minimumConnections - num_connections_from_node); ++i){
				boost::add_edge(src, minimaIdxAndVals[i].first, minimaIdxAndVals[i].second, dg);
			}
		}
#endif
	}
	return dg;
}
static void printGraph(DirectedGraph_t const &dg){
	vertex_iterator_t vit, vend;
	std::tie(vit, vend) = boost::vertices(dg);
	for (; vit != vend; ++vit){
		DirectedGraph_t::out_edge_iterator eit, eend;
		std::tie(eit, eend) = boost::out_edges(*vit, dg);
		std::for_each(eit, eend,
		  [&dg](edge_descriptor_t it)
			{ std::cout << "           " << boost::target(it, dg) <<  '\n'; });
	}
}
static void printEdges(DirectedGraph_t &dg){
	std::cout << "printing edges: \n";
	for (std::pair<edge_iterator_t, edge_iterator_t>edgePair = edges(dg); edgePair.first != edgePair.second; ++edgePair.first)
	{
		dg[*(edgePair.first)].timbral_distance = 42;
	   std::cout << *edgePair.first << " " << dg[*(edgePair.first)].timbral_distance << std::endl;
	}
	std::cout << "done printing edges\n";
}
static double getOutEdgeWeightPercentile(DirectedGraph_t &dg, vertex_iterator_t vit, double percentile = 0.5){
	DirectedGraph_t::out_edge_iterator eit, eend;
	std::tie(eit, eend) = boost::out_edges(*vit, dg);
	
	std::vector<double> distances;
	for (; eit != eend; ++eit){
		edge_descriptor_t edge_dscr = *eit;

		distances.push_back( dg[edge_dscr].timbral_distance );
	}
	std::sort(distances.begin(), distances.end());
	size_t idxOfPrctl = size_t((double)distances.size() * percentile);
	double edgeWeightPercentile = 0.0;
	if (!distances.empty()) {
		edgeWeightPercentile = distances[idxOfPrctl];
	}
	
	return edgeWeightPercentile;
}
static vertex_iterator_t vertexDescriptorToIterator(DirectedGraph_t const &dg, vertex_descriptor_t const v_dscr){
	vertex_iterator_t vit, vend;
	std::tie(vit, vend) = boost::vertices(dg);
	
	vertex_iterator_t v_target_it = std::find (vit, vend, v_dscr);
	return v_target_it;
}
static vertex_descriptor_t traverseToNearestVertex(DirectedGraph_t &dg, vertex_iterator_t vit){
	DirectedGraph_t::out_edge_iterator eit, eend;
	std::tie(eit, eend) = boost::out_edges(*vit, dg);
	double leastDistance = 1e15;
	vertex_descriptor_t closestVertex = *vit;
	for (; eit != eend; ++eit){
		edge_descriptor_t edge_dscr = *eit;
		float distance = dg[edge_dscr].timbral_distance ;
		if (distance < leastDistance){
			leastDistance = distance;
			closestVertex = boost::target(*eit, dg);
		}
	}
	return closestVertex;
}
static void printProbs(std::vector<double> const &probs){
	const char on = '|';
	const char off = '_';
	for (auto d : probs){
		const unsigned int val = static_cast<unsigned int>((d * 10.0) + 0.5);
		for (unsigned i = 0; i < val; ++i)
			std::cout << on;
		for (unsigned i = 0; i < (10 - val); ++i)
			std::cout << off;
		std::cout << ' ';
	}
	std::cout << '\n';
}
static inline size_t rollWeightedDie(std::vector<double> const &probabilities, double power = 1.0) {
	std::vector<double> tmp(probabilities.size());
	double min = *(std::min_element(probabilities.begin(), probabilities.end()));
	double max = *(std::max_element(probabilities.begin(), probabilities.end()));
	std::transform(probabilities.begin(), probabilities.end(), tmp.begin(), [min, max](double x){
		const double rng = max - min;
		if (rng > 0.0){
			x -= min;
			x /= rng;
		} else {
			std::cerr << "0 range?!\n";
		}
		return x;
	});
	
	const static auto func1 = [power](const double x){
		double xpt = std::pow(x, power);
		double mxpt = std::pow((1.0 - x), power);
		return xpt / (xpt + mxpt);
	};
	const static auto func2 = [power](const double x){
		return std::exp(x * power);
	};
	std::transform(tmp.begin(), tmp.end(), tmp.begin(), func2);
	
	double sum = std::accumulate(tmp.begin(), tmp.end(), 0.0);
	if (sum > 0.0){
		std::transform(tmp.begin(), tmp.end(), tmp.begin(), [sum](const double d){
			return d / sum;
		});
	}
	printProbs(tmp);
	
	discrete_distr_t dist(std::begin(tmp), std::end(tmp));
	// HERE is where it matters that i'm changing the gaussian
	boost::variate_generator< random_gen_t&, discrete_distr_t > weightsampler(gen, dist);
	return weightsampler();
}
static inline std::vector<double> getProbabilitiesFromCurrentNode(DirectedGraph_t const &dg, vertex_descriptor_t current_vertex, const float C_kernelScaling){
	int i = 0;
	auto es = boost::out_edges(current_vertex, dg);
	const size_t num_options = es.second - es.first;
	std::vector<double> probabilities(num_options);
	
	for (auto eit = es.first; eit != es.second; ++eit) {
//            std::cout << boost::source(*eit, g) << ' ' << boost::target(*eit, g) << std::endl;
		edge_property *e_prop = (edge_property *)eit->get_property();
//            std::cout << hgf->edge_weight << std::endl;
//            probabilities.push_back(hgf->edge_weight);
		double distance = e_prop->timbral_distance;
//		distance = std::max(distance, 0.000000001);
		probabilities[i] = std::exp(-(distance*distance) / 2.0 * (C_kernelScaling*C_kernelScaling));
		i++;
	}
	
	double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
	if (sum > 0.0){
		std::transform(probabilities.begin(), probabilities.end(), probabilities.begin(),
		[sum](const double d){
			return d / sum;
		});
	}
	
	return probabilities;
}
static vertex_descriptor_t traverseToRandomVertex(DirectedGraph_t &dg, vertex_descriptor_t current_vertex, const float C_kernelScaling = 12.f, double probPower = 1.0){
//        std::cout << "current node: " << g[current_node].name << "\n";
		//this->previous_node = this->current_node;
	std::pair<adjacency_iter_t, adjacency_iter_t> adjacent_verts = boost::adjacent_vertices(current_vertex, dg);
	adjacency_iter_t current_adjacent = adjacent_verts.first;
	adjacency_iter_t last_adjacent = adjacent_verts.second;

	std::vector<double> probabilities = getProbabilitiesFromCurrentNode(dg, current_vertex, C_kernelScaling);

	// IF NUM OPTIONS IS 0, 2 PROBLEMS:
	// WE GET A SIZE 0 VECTOR
	// THE TRAVERSAL WILL ALWAYS SELF-TRANSITION
	size_t num_options = last_adjacent - current_adjacent;
	if (num_options == 0){
		return current_vertex;
	}
	std::vector<vertex_descriptor_t> next_vertices_options(num_options);

	for (; current_adjacent != last_adjacent; ++current_adjacent){
		// this method counts down
		size_t idx = last_adjacent - current_adjacent - 1;
//            next_vertices_options.push_back(*current_adjacent);
//            std::cout << "adj diff: " << last_adjacent - current_adjacent << std::endl;
		next_vertices_options[idx] = *current_adjacent;
//            std::cout << "adjacent to node" << current_node << ": " << *current_adjacent << "\n";
	}

	size_t new_idx = rollWeightedDie(probabilities, probPower);
	//    this->current_node = next_vertices_options[new_idx];
	vertex_descriptor_t next_node = next_vertices_options[new_idx];
//        if (next_node > num_nodes){ std::cerr << "NODE ABOVE LIMITS" << "\n"; }
//        std::cout << "next node: " << next_node << "new idx: " << new_idx << "\n";
	return next_node;
}
static void reduceEdges(DirectedGraph_t &dg, double percentile = 0.5){
	percentile = 1.0 - percentile;	// in terms of amount to remove now
	std::cout << "reduceEdges\n";
	vertex_iterator_t vit, vend;
	std::tie(vit, vend) = boost::vertices(dg);

	for (; vit != vend; ++vit){
		double avgDist = getOutEdgeWeightPercentile(dg, vit, percentile);
		auto pred = [&dg, avgDist](edge_descriptor_t ed){
			double dist = dg[ed].timbral_distance;
			return (dist > avgDist) || (dist < 0.1);// if true it will be deleted
		};
		vertex_descriptor_t u = *vit;
		remove_out_edge_if(u, pred, dg);
	}
	std::cout << "reduceEdges complete\n";
}
static void removeProportionOfVertices(DirectedGraph_t &dg, double proportion){
	unsigned long  N = boost::num_vertices(dg);
	unsigned long numToRemove = (unsigned long)((double)N * proportion);
	unsigned long c = N / numToRemove;
	for (unsigned long i = 0; i < numToRemove; i ++){
		auto idx = N - 1 - (i * c);
		clear_vertex(idx, dg);	// clears vertex's edges first, otherwise we may get loops or other unexpected results
		remove_vertex(idx, dg);
	}
}
	  
static void exportGraphAsDot(DirectedGraph_t &dg, std::string const &filename){
	std::ostream& strm = (std::cout);
	std::ofstream file(filename);
	
	// save output buffer of the stream
	std::streambuf* strm_buffer = strm.rdbuf();

	// redirect ouput into the file
	strm.rdbuf (file.rdbuf());

	const std::vector<std::string> colors = {
		"\"#40e0d0\"",
		"\"#000000\"",
		"\"#ff0000\"",
		"\"#40e0d0\"",
		"\"#a0522d\"",
		"\"#6557d2\"",
		"\"#9b6da4\""
	};
	unsigned short colorIdx = 0;
	
	auto vw = [](std::ostream& out, const vertex_descriptor_t& v){};
	auto ew = [&dg, &colors, &colorIdx](std::ostream& out, const edge_descriptor_t& ed){
//		out << "[label=\"" << "bo" << "\"]";
		auto val = dg[ed].timbral_distance;
		auto label = val;
		label = ((int)(label * 100.0)) / 100.0;
		val *= val;
		val = std::max(val, 0.1f);
		val = 1.0 / val;
		out << std::fixed << "[weight=" << val << ",label=" << label << ",color=" << colors[colorIdx] << "]";
		++colorIdx;
		colorIdx %= colors.size();
	};
	write_graphviz(std::cout, dg, vw, ew);
	  
	// restore old output buffer
	strm.rdbuf (strm_buffer);
}


} // end namespace nvs

#endif /* graph_h */
