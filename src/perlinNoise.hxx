/*! \file perlinNoise.hxx
 *  \brief breastPhantom Perlin Noise header file
 *  \author Christian G. Graff
 *  \version 1.0
 *  \date 2018
 *
 *  \copyright To the extent possible under law, the author(s) have
 *  dedicated all copyright and related and neighboring rights to this
 *  software to the public domain worldwide. This software is
 *  distributed without any warranty.  You should have received a copy
 *  of the CC0 Public Domain Dedication along with this software.
 *  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 */

#ifndef __PERLINNOISE_HXX__
#define __PERLINNOISE_HXX__

#include <cstdint>
#include <boost/program_options.hpp>


class perlinNoise final {
    using variables_map = boost::program_options::variables_map;

private:
    double frequency;
    double lacunarity;
    double persistence;
    unsigned numOctaves;
    int32_t seed;
    int32_t xNoiseGen, yNoiseGen, zNoiseGen, seedNoiseGen, shiftNoiseGen;

public:
    perlinNoise(const variables_map &vm, int32_t inSeed, double freq, double lac, double pers, unsigned oct);
    perlinNoise(const variables_map &vm, int32_t inSeed, const char type[]);
    perlinNoise(const variables_map &vm, const char type[]);

    inline void setSeed(int32_t inSeed) {
        seed = inSeed;
    }

    double getNoise(const double r[3]) const;
};

#endif /* __PERLINNOISE_HXX__ */
