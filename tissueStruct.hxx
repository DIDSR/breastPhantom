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

// tissue types

typedef struct {
	unsigned char bg;
	unsigned char skin;
	unsigned char nipple;
	unsigned char fat;
	unsigned char cooper;
	unsigned char gland;
	unsigned char TDLU;
	unsigned char duct;
	unsigned char artery;
	unsigned char vein;
	unsigned char muscle;
} tissueStruct;


