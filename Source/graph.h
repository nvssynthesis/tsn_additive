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
namespace sgt{	// stochastic graph traversal
typedef std::vector<std::vector<essentia::Real>> multiDimRealVec_t;

struct vertex_property{
	double timestamp;
	unsigned int cluster;
	float const *fundamentalFreq;
	float const *voicedness;
	float const *loudness;
	std::vector<float> const *PCAvec;
};
struct edge_property{
    float timbral_distance;
    float temporal_distance;
	float pitch_distance;
	float loudness_distance;
};
typedef boost::adjacency_list < boost::vecS, // out-edges, a Sequence or an AssociativeContainer
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
	std::transform(x0.begin(), x0.end(), x1.begin(), diff.begin(), op);
	d = std::accumulate(diff.begin(), diff.end(), 0.f);
	d = std::sqrtf(d);
	return d;
}
static std::optional<float> kullback_leibler_divergence(std::vector<float> const &x0, std::vector<float> const &x1){
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

static float getMaxTimbralDistanceForConnection(multiDimRealVec_t const &pca_mat, float percentile = 0.15f, size_t nSearch = 100UL)	{
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
static float getMaxTimbralDistanceForConnectionFromVertex(const vertex_descriptor_t v, multiDimRealVec_t const &pca_mat, float percentile = 0.15f, size_t nSearch = 100UL){
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

// create directed graph based on numFrames
static DirectedGraph_t createGraphFromAnalysisData(AnalysisData const &data, connectionHeuristics const settings = connectionHeuristics()){
	const multiDimRealVec_t &pca_mat = data.getPCAmat();
	
	
	const auto &freqs = data.f0;
	const auto &voicedness = data.voicedProbability;
	const auto &loudness = data.loudness[0];
	const auto &PCAs = data.getPCAmat();
	const size_t numFrames = pca_mat.size();
	
	assert(freqs.size() == numFrames);
	assert(voicedness.size() == numFrames);
	assert(loudness.size() == numFrames);
	assert(PCAs.size() == numFrames);

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
		dg[src].fundamentalFreq = &freqs[src];
		dg[src].voicedness = &voicedness[src];
		dg[src].loudness = &loudness[src];
		dg[src].PCAvec = &PCAs[src];

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
		for (int dst = 0; dst < numFrames; dst += int(1 + (numFrames / 1000))){
			if (src == dst)
				continue;
			auto x1 = pca_mat[dst];
			edge_property ep;
			ep.temporal_distance = (float)std::abs(dst - src);
			ep.timbral_distance = euc_distance(x0, x1).value_or(-1.f);
			ep.pitch_distance = pitchDistanceRatio(*(dg[src].fundamentalFreq), freqs[dst]);
			ep.loudness_distance = loudnessDistance(*(dg[src].loudness), loudness[dst]);
			
			if ((ep.timbral_distance < maxTimbralDistanceForConnection)
				& (pitchDistanceInRange(*(dg[src].fundamentalFreq), freqs[dst], settings))
				& (ep.loudness_distance < maxLoudnessDev)
				)
			{
				boost::add_edge(src, dst, ep, dg);
				++num_connections_from_node;
			}
			else if (auto maxOfMins = std::max_element(minimaIdxAndVals.begin(), minimaIdxAndVals.end(), idxAndEdgeProp_CmpByTimbre); ep.timbral_distance < maxOfMins->second.timbral_distance){
				// only add edge property to minimal list if it's not connected
				maxOfMins->first = dst;
				maxOfMins->second = ep;
			}
		}
#if 1
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
	for (auto [edge_st, edge_end] = edges(dg); edge_st != edge_end; ++edge_st)
	{
		dg[*edge_st].timbral_distance = 42;
	   std::cout << *edge_st << " " << dg[*edge_st].timbral_distance << '\n';
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
	for (auto &p : probs)
		std::cout << p << ' ';
	std::cout << ";\n";
}
static inline void exaggerateProbabilities(std::vector<double> &probs, double power){
	std::transform(probs.begin(), probs.end(), probs.begin(), [power](const double x){
		double xpt = std::pow(x, power);
		double mxpt = std::pow((1.0 - x), power);
		return xpt / (xpt + mxpt);
	});
}
static inline void normalizeProbabilities(std::vector<double> &probs){
	double sum = std::accumulate(probs.begin(), probs.end(), 0.0);
	if (sum > 0.0){
		std::transform(probs.begin(), probs.end(), probs.begin(), [sum](const double d){
			return d / sum;
		});
	}
}
static inline size_t rollWeightedDie(std::vector<double> const &probs) {
	discrete_distr_t dist(probs.begin(), probs.end());
	// HERE is where it matters that i'm changing the gaussian
	boost::variate_generator< random_gen_t&, discrete_distr_t > weightsampler(gen, dist);
	return weightsampler();
}
static inline std::vector<double> getProbabilitiesFromCurrentNode(DirectedGraph_t const &dg, vertex_descriptor_t current_vertex, const float C_kernelScaling){
	int i = 0;
	auto [out_edge_st, out_edge_end] = boost::out_edges(current_vertex, dg);
	const size_t num_options = out_edge_end - out_edge_st;
	std::vector<double> probabilities(num_options);
	
	for (auto eit = out_edge_st; eit != out_edge_end; ++eit) {
		edge_property *e_prop = (edge_property *)eit->get_property();
		double distance = e_prop->timbral_distance;
		probabilities[i] = std::exp(-(distance*distance) / (2.0 * (C_kernelScaling*C_kernelScaling)));
		i++;
	}
	return probabilities;
}
static vertex_descriptor_t traverseToRandomVertex(DirectedGraph_t &dg, vertex_descriptor_t current_vertex, const float C_kernelScaling = .01f, double probPower = 1.0){
	
	auto [current_adjacent, last_adjacent] = boost::adjacent_vertices(current_vertex, dg);
	// IF NUM OPTIONS IS 0, 2 PROBLEMS:
	// WE GET A SIZE 0 VECTOR
	// THE TRAVERSAL WILL ALWAYS SELF-TRANSITION
	size_t num_options = last_adjacent - current_adjacent;
	if (num_options == 0){
		return current_vertex;
	}
	std::vector<double> probs = getProbabilitiesFromCurrentNode(dg, current_vertex, C_kernelScaling);
	exaggerateProbabilities(probs, probPower);
	normalizeProbabilities(probs);
	
	std::vector<vertex_descriptor_t> next_vertices_options(num_options);
	for (; current_adjacent != last_adjacent; ++current_adjacent){
		// this method counts down
		size_t idx = last_adjacent - current_adjacent - 1;
		next_vertices_options[idx] = *current_adjacent;
	}
	size_t new_idx = rollWeightedDie(probs);
	vertex_descriptor_t next_node = next_vertices_options[new_idx];
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

} // end namespace sgt
} // end namespace nvs

#endif /* graph_h */
