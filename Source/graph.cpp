//
//  graph.cpp
//
//  Created by Nicholas Solem on 5/15/23.
//  Copyright © 2023 Nicholas Solem. All rights reserved.
//

#include <stdio.h>

#include "graph.h"

 std::vector<size_t> nvs::sgt::getRandomIndices(size_t maxValInclusive, size_t numVals) {
//	std::random_device rd; // obtain a random number from hardware
//	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<size_t> distr(0, maxValInclusive); // define the range
	
	std::vector<size_t> randVec (numVals, 0);
	for(size_t n = 0; n < numVals; ++n){
		size_t val = distr(gen);
		randVec[n] = val;
	}
	return randVec;
}

 float nvs::sgt::getMaxTimbralDistanceForConnection(multiDimRealVec_t const &pca_mat, float percentile, size_t nSearch){
	const size_t numFrames = pca_mat.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	for (int src = 0; src < nSearch; ++src){
		auto x0 = pca_mat[randIdx[src]];
		for (int dst = 0; dst < nSearch; ++dst){
			if (src == dst)
				continue;
			float d = nvs::dst::euc_distance(x0, pca_mat[randIdx[dst]]).value_or(-1.f);
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

 float nvs::sgt::getMaxTimbralDistanceForConnectionFromVertex(const vertex_descriptor_t v, multiDimRealVec_t const &pca_mat, float percentile, size_t nSearch){
	auto src_timbre = pca_mat[v];
	const size_t numFrames = pca_mat.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	distances.reserve(nSearch);
	for (auto dst = 0; dst < nSearch; ++dst){
		if (dst == v)
			continue;
		float d = nvs::dst::euc_distance(src_timbre, pca_mat[randIdx[dst]]).value_or(-1.f);
		distances.push_back(d);
	}
	size_t pctlIdx = size_t((distances.size()-1) * percentile);
	nth_element( distances.begin(), distances.begin() + pctlIdx, distances.end());
	
	return *(distances.begin() + pctlIdx);
}

 float nvs::sgt::getMaxLoudnessDistanceForConnection(std::vector<float> const &loudnessVec, float percentile, size_t nSearch){
	const size_t numFrames = loudnessVec.size();
	nSearch = std::min(nSearch, numFrames);
	std::vector<size_t> randIdx = nvs::sgt::getRandomIndices(numFrames-1, nSearch);
	std::vector<float> distances;
	for (int src = 0; src < nSearch; ++src){
		auto x0 = loudnessVec[randIdx[src]];
		for (int dst = 0; dst < nSearch; ++dst){
			if (src == dst)
				continue;
			float d = nvs::dst::loudnessDistance(x0, loudnessVec[randIdx[dst]]);
			distances.push_back(d);
		}
	}
	std::sort(distances.begin(), distances.end());
	size_t pctlIdx = size_t((distances.size()-1) * percentile);
	float maxDistanceForConnection = distances[pctlIdx];
	std::cout << "max loudness DistanceForConnection: " << maxDistanceForConnection << '\n';
	return maxDistanceForConnection;
}


 nvs::sgt::DirectedGraph_t nvs::sgt::createGraphFromAnalysisData(AnalysisData const &data, connectionHeuristics const settings){
	const multiDimRealVec_t &pca_mat = data.PCAmat;
	
	const auto &freqs = data.f0;
	const auto &voicedness = data.voicedProbability;
	const auto &loudness = data.loudness[0];
	const size_t numFrames = data.numFrames;
	
	assert(freqs.size() == numFrames);
	assert(voicedness.size() == numFrames);
	assert(loudness.size() == numFrames);
	assert(pca_mat.size() == numFrames);

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
		dg[src].PCAvec = &pca_mat[src];

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
			ep.timbral_distance = nvs::dst::euc_distance(x0, x1).value_or(-1.f);
			ep.pitch_distance = nvs::dst::pitchDistanceRatio(*(dg[src].fundamentalFreq), freqs[dst]);
			ep.loudness_distance = nvs::dst::loudnessDistance(*(dg[src].loudness), loudness[dst]);
			
			if ((ep.timbral_distance < maxTimbralDistanceForConnection)
				& (nvs::dst::pitchDistanceInRange(*(dg[src].fundamentalFreq), freqs[dst], settings))
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



 void nvs::sgt::printEdges(nvs::sgt::DirectedGraph_t &dg){
	std::cout << "printing edges: \n";
	for (auto [edge_st, edge_end] = edges(dg); edge_st != edge_end; ++edge_st)
	{
		std::cout << *edge_st << " " << dg[*edge_st].timbral_distance << '\n';
	}
	std::cout << "done printing edges\n";
}
 void nvs::sgt::printGraph(nvs::sgt::DirectedGraph_t const &dg){
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

 double nvs::sgt::getOutEdgeWeightPercentile(nvs::sgt::DirectedGraph_t &dg, vertex_iterator_t vit, double percentile){
	nvs::sgt::DirectedGraph_t::out_edge_iterator eit, eend;
	std::tie(eit, eend) = boost::out_edges(*vit, dg);
	
	std::vector<double> distances;
	for (; eit != eend; ++eit){
		nvs::sgt::edge_descriptor_t edge_dscr = *eit;

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
 nvs::sgt::vertex_iterator_t nvs::sgt::vertexDescriptorToIterator(nvs::sgt::DirectedGraph_t const &dg, nvs::sgt::vertex_descriptor_t const v_dscr){
	nvs::sgt::vertex_iterator_t vit, vend;
	std::tie(vit, vend) = boost::vertices(dg);
	
	nvs::sgt::vertex_iterator_t v_target_it = std::find (vit, vend, v_dscr);
	return v_target_it;
}

 nvs::sgt::vertex_descriptor_t nvs::sgt::traverseToNearestVertex(nvs::sgt::DirectedGraph_t &dg, nvs::sgt::vertex_iterator_t vit){
	nvs::sgt::DirectedGraph_t::out_edge_iterator eit, eend;
	std::tie(eit, eend) = boost::out_edges(*vit, dg);
	double leastDistance = 1e15;
	nvs::sgt::vertex_descriptor_t closestVertex = *vit;
	for (; eit != eend; ++eit){
		nvs::sgt::edge_descriptor_t edge_dscr = *eit;
		float distance = dg[edge_dscr].timbral_distance ;
		if (distance < leastDistance){
			leastDistance = distance;
			closestVertex = boost::target(*eit, dg);
		}
	}
	return closestVertex;
}

 void nvs::sgt::printProbs(std::vector<double> const &probs){
	for (auto &p : probs)
		std::cout << p << ' ';
	std::cout << ";\n";
}

 inline void nvs::sgt::exaggerateProbabilities(std::vector<double> &probs, double power){
	std::transform(probs.begin(), probs.end(), probs.begin(), [power](const double x){
		double xpt = std::pow(x, power);
		double mxpt = std::pow((1.0 - x), power);
		return xpt / (xpt + mxpt);
	});
}
 inline void nvs::sgt::normalizeProbabilities(std::vector<double> &probs){
	double sum = std::accumulate(probs.begin(), probs.end(), 0.0);
	if (sum > 0.0){
		std::transform(probs.begin(), probs.end(), probs.begin(), [sum](const double d){
			return d / sum;
		});
	}
}
 inline size_t nvs::sgt::rollWeightedDie(std::vector<double> const &probs) {
	nvs::sgt::discrete_distr_t dist(probs.begin(), probs.end());
	// HERE is where it matters that i'm changing the gaussian
	boost::variate_generator< nvs::sgt::random_gen_t&, nvs::sgt::discrete_distr_t > weightsampler(gen, dist);
	return weightsampler();
}

 inline std::vector<double> nvs::sgt::getProbabilitiesFromCurrentNode(nvs::sgt::DirectedGraph_t const &dg, nvs::sgt::vertex_descriptor_t current_vertex, const float C_kernelScaling){
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

 nvs::sgt::vertex_descriptor_t nvs::sgt::traverseToRandomVertex(nvs::sgt::DirectedGraph_t &dg, nvs::sgt::vertex_descriptor_t current_vertex, const float C_kernelScaling, double probPower){
	
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

 void nvs::sgt::reduceEdges(nvs::sgt::DirectedGraph_t &dg, double percentile){
	percentile = 1.0 - percentile;	// in terms of amount to remove now
	std::cout << "reduceEdges\n";
	nvs::sgt::vertex_iterator_t vit, vend;
	std::tie(vit, vend) = boost::vertices(dg);

	for (; vit != vend; ++vit){
		double avgDist = nvs::sgt::getOutEdgeWeightPercentile(dg, vit, percentile);
		auto pred = [&dg, avgDist](nvs::sgt::edge_descriptor_t ed){
			double dist = dg[ed].timbral_distance;
			return (dist > avgDist) || (dist < 0.1);// if true it will be deleted
		};
		nvs::sgt::vertex_descriptor_t u = *vit;
		remove_out_edge_if(u, pred, dg);
	}
	std::cout << "reduceEdges complete\n";
}

 void nvs::sgt::removeProportionOfVertices(nvs::sgt::DirectedGraph_t &dg, double proportion){
	unsigned long  N = boost::num_vertices(dg);
	unsigned long numToRemove = (unsigned long)((double)N * proportion);
	unsigned long c = N / numToRemove;
	for (unsigned long i = 0; i < numToRemove; i ++){
		auto idx = N - 1 - (i * c);
		clear_vertex(idx, dg);	// clears vertex's edges first, otherwise we may get loops or other unexpected results
		remove_vertex(idx, dg);
	}
}

 void nvs::sgt::exportGraphAsDot(nvs::sgt::DirectedGraph_t &dg, std::string const &filename){
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

class nvs::sgt::StatefulDirectedGraph{
	nvs::sgt::DirectedGraph_t m_dg;
};

