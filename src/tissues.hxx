/*! \file tissueStruct.hxx
 *  \brief breastPhantom tissue header file
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

#ifndef __TISSUES_HXX__
#define __TISSUES_HXX__

// tissue types

namespace tissue {
    const static unsigned char bg = 0;
    const static unsigned char skin = 2;
    const static unsigned char nipple = 33;
    const static unsigned char fat = 1;
    const static unsigned char cooper = 88;
    const static unsigned char gland = 29;
    const static unsigned char TDLU = 95;
    const static unsigned char duct = 125;
    const static unsigned char artery = 150;
    const static unsigned char vein = 225;
    const static unsigned char muscle = 40;

    // unligamented tissue classes
    const static unsigned char ufat = 60;
    const static unsigned char ugland = 61;
    const static unsigned char uTDLU = 62;
    const static unsigned char uduct = 63;
}

#endif /* __TISSUES_HXX__ */
