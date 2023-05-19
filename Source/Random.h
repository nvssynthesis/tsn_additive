//
//  Random.h
//  tsara_additive
//
//  Created by Nicholas Solem on 5/17/23.
//  Copyright © 2023 nvssynthesis. All rights reserved.
//

#ifndef Random_h
#define Random_h
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>

namespace nvs	{
namespace rand	{
typedef boost::random::mt19937 random_gen_t;
typedef boost::random::discrete_distribution<size_t,double> discrete_distr_t;
typedef boost::normal_distribution<double> normal_distr_t;

static random_gen_t gen;

static std::vector<size_t> getRandomIndices(size_t maxValInclusive, size_t numVals = 100){
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

inline size_t rollWeightedDie(std::vector<double> const &probs){
	discrete_distr_t dist(probs.begin(), probs.end());
	// HERE is where it matters that i'm changing the gaussian
	boost::variate_generator< random_gen_t&, discrete_distr_t > weightsampler(gen, dist);
	return weightsampler();
}


}	// namespace rand
}	// namespace nvs

#endif /* Random_h */
