/*
 * breastPhantom.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

// create volumetric breast

#include "breastPhantom.hxx"

namespace po = boost::program_options;

unsigned int ductTree::num = 0;
unsigned int arteryTree::num = 0;
unsigned int veinTree::num = 0;

int main(int argc, char* argv[]){

	double pi = vtkMath::Pi();
	
	// timing
	time_t currentTime, previousTime;
	time(&previousTime);
	time_t startTime, endTime;
	time(&startTime);

	// configuration file variables
	po::options_description baseOpt("base options");
	baseOpt.add_options()
		("base.outputDir",po::value<std::string>()->default_value("."),"output directory")
		("base.imgRes",po::value<double>()->default_value(0.25),"phantom resolution (mm)")
		("base.skinThick",po::value<double>()->default_value(0.75),"skin thickness (mm)")
		("base.nippleLen",po::value<double>()->default_value(4.0),"nipple length (mm)")
		("base.nippleRad",po::value<double>()->default_value(4.0),"nipple radius (mm)")
		("base.areolaRad",po::value<double>()->default_value(8.0),"areola radius (mm)")
		("base.leftBreast",po::value<bool>()->default_value(true),"left side breast (boolean)")
		("base.targetFatFrac",po::value<double>()->default_value(0.75),"desired fraction of breast to be fat")
		("base.seed",po::value<unsigned int>(),"random number generator seed")
	;

	po::options_description shapeOpt("breast shape options");
	shapeOpt.add_options()
		("shape.ures",po::value<double>()->default_value(0.02),"u resolution of base shape")
		("shape.vres",po::value<double>()->default_value(0.02),"v resolution of base shape")
		("shape.pointSep",po::value<double>()->default_value(0.01),"minimum point separation (mm)")
		("shape.ringWidth",po::value<double>()->default_value(10.0),"thickness of back ring (mm)")
		("shape.ringSep",po::value<double>()->default_value(0.5),"ring node step size (mm)")
		("shape.featureAngle",po::value<double>()->default_value(20.0),"angle to preserve while smoothing (degrees)")
		("shape.targetReduction",po::value<double>()->default_value(0.05),"fraction of triangles to decimate")
		("shape.a1b",po::value<double>()->default_value(1.0),"bottom size")
		("shape.a1t",po::value<double>()->default_value(1.0),"top size")
		("shape.a2l",po::value<double>()->default_value(1.0),"left size")
		("shape.a2r",po::value<double>()->default_value(1.0),"right size")
		("shape.a3",po::value<double>()->default_value(1.0),"outward size")
		("shape.eps1",po::value<double>()->default_value(1.0),"u quadric exponent")
		("shape.eps2",po::value<double>()->default_value(1.0),"v quadric exponent")
		("shape.doPtosis",po::value<bool>()->default_value(false),"ptosis deformation (boolean)")
		("shape.ptosisB0",po::value<double>()->default_value(0.4),"ptosis b0")
		("shape.ptosisB1",po::value<double>()->default_value(0.3),"ptosis b1")
		("shape.doTurn",po::value<bool>()->default_value(false),"turn deformation (boolean)")
		("shape.turnC0",po::value<double>()->default_value(-0.4),"turn c0")
		("shape.turnC1",po::value<double>()->default_value(0.2),"turn c1")
		("shape.doTopShape",po::value<bool>()->default_value(false),"top shape deformation (boolean)")
		("shape.topShapeS0",po::value<double>()->default_value(1.54),"top shape s0")
		("shape.topShapeS1",po::value<double>()->default_value(1.3),"top shape s1")
		("shape.topShapeT0",po::value<double>()->default_value(-3.0),"top shape t0")
		("shape.topShapeT1",po::value<double>()->default_value(-1.0),"top shape t1")
		("shape.doFlattenSide",po::value<bool>()->default_value(false),"flatten side deformation (boolean)")
		("shape.flattenSideG0",po::value<double>()->default_value(-0.4),"flatten side g0")
		("shape.flattenSideG1",po::value<double>()->default_value(0.2),"flatten side g1")
		("shape.doTurnTop",po::value<bool>()->default_value(false),"turn top deformation (boolean)")
		("shape.turnTopH0",po::value<double>()->default_value(-2.7),"turn top h0")
		("shape.turnTopH1",po::value<double>()->default_value(-3.7),"turn top h1")
	;

	po::options_description compartOpt("breast compartment options");
	compartOpt.add_options()
		("compartments.num",po::value<int>()->default_value(8),"number of breast compartments")
		("compartments.seedBaseDist",po::value<double>()->default_value(12.5),"distance along nipple line of compartment seed base (mm)")
		("compartments.backFatBufferFrac",po::value<double>()->default_value(0.2),"fraction of phantom in nipple direction forced to be fat")
		("compartments.numBackSeeds",po::value<int>()->default_value(100),"number of backplane seed points")
		("compartments.angularJitter",po::value<double>()->default_value(0.125),"maximum seed jitter (fraction of subtended angle)")
		("compartments.zJitter",po::value<double>()->default_value(5.0),"maximum seed jitter in nipple direction (mm)")
		("compartments.maxFracRadialDist",po::value<double>()->default_value(0.5),"maximum radial distance from base seed as a fraction of distance to breast surface")
		("compartments.minFracRadialDist",po::value<double>()->default_value(0.25),"minimum radial distance from base seed as a fraction of distance to breast surface")
		("compartments.minScaleNippleDir",po::value<double>()->default_value(0.01),"minimum scale in nipple direction")
		("compartments.maxScaleNippleDir",po::value<double>()->default_value(0.01),"maximum scale in nipple direction")
		("compartments.minScale",po::value<double>()->default_value(50.0),"minimum scale in non-nipple direction")
		("compartments.maxScale",po::value<double>()->default_value(60.0),"maximum scale in non-nipple direction")
		("compartments.minGlandStrength",po::value<double>()->default_value(20.0),"minimum gland strength")
		("compartments.maxGlandStrength",po::value<double>()->default_value(20.0),"maximum gland strength")
		("compartments.maxDeflect",po::value<double>()->default_value(0.01),"maximum compartment deflection angle from pointing towards nipple (fraction of pi)")
		("compartments.minSkinScaleNippleDir",po::value<double>()->default_value(0.01),"minimum scale skin seeds in nipple direction")
		("compartments.maxSkinScaleNippleDir",po::value<double>()->default_value(0.01),"maximum scale skin seeds in nipple direction")
		("compartments.minSkinScale",po::value<double>()->default_value(70.0),"minimum scale skin in non-nipple direction")
		("compartments.maxSkinScale",po::value<double>()->default_value(70.0),"maximum scale skin in non-nipple direction")
		("compartments.skinStrength",po::value<double>()->default_value(20.0),"skin strength")
		("compartments.backScale",po::value<double>()->default_value(50.0),"back scale")
		("compartments.backStrength",po::value<double>()->default_value(50.0),"back strength")
		("compartments.nippleScale",po::value<double>()->default_value(50.0),"nipple scale")
		("compartments.nippleStrength",po::value<double>()->default_value(50.0),"nipple strength")
		("compartments.voronSeedRadius",po::value<double>()->default_value(40.0),"check all seeds in radius (mm)")
	;

	po::options_description TDLUOpt("TDLU options");
	TDLUOpt.add_options()
		("TDLU.maxLength",po::value<double>()->default_value(4.0),"maximum TDLU length")
		("TDLU.minLength",po::value<double>()->default_value(2.0),"minimum TDLU length")
		("TDLU.maxWidth",po::value<double>()->default_value(2.0),"maximum TDLU width")
		("TDLU.minWidth",po::value<double>()->default_value(1.0),"minimum TDLU width")
	;

	po::options_description perlinOpt("Perlin noise options");
	perlinOpt.add_options()
		("perlin.maxDeviation",po::value<double>()->default_value(0.1),"maximum perlin perturbation fraction of radius")
		("perlin.frequency",po::value<double>()->default_value(0.5),"starting frequency")
		("perlin.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
		("perlin.persistence",po::value<double>()->default_value(0.2),"octave signal decay")
		("perlin.numOctaves",po::value<int>()->default_value(6),"number of frequency octaves")
		("perlin.xNoiseGen",po::value<int>()->default_value(683),"x direction noise generation seed")
		("perlin.yNoiseGen",po::value<int>()->default_value(4933),"y direction noise generation seed")
		("perlin.zNoiseGen",po::value<int>()->default_value(23),"z direction noise generation seed")
		("perlin.seedNoiseGen",po::value<int>()->default_value(3095),"seed noise generation")
		("perlin.shiftNoiseGen",po::value<int>()->default_value(11),"shift noise generation seed")
	;

	po::options_description boundaryOpt("Boundary noise options");
	boundaryOpt.add_options()
		("boundary.maxDeviation",po::value<double>()->default_value(0.25),"maximum perlin perturbation fraction of radius")
		("boundary.frequency",po::value<double>()->default_value(0.2),"starting frequency")
		("boundary.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
		("boundary.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
	;

	po::options_description perturbOpt("Lobule perturbation noise options");
	perturbOpt.add_options()
		("perturb.maxDeviation",po::value<double>()->default_value(0.25),"maximum perlin perturbation fraction of radius")
		("perturb.frequency",po::value<double>()->default_value(0.2),"starting frequency")
		("perturb.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
		("perturb.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
	;
	
	po::options_description bufferOpt("Lobule glandular buffer options");
	bufferOpt.add_options()
		("buffer.maxDeviation",po::value<double>()->default_value(0.1),"maximum perlin buffer fraction of radius")
		("buffer.frequency",po::value<double>()->default_value(0.2),"starting frequency")
		("buffer.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
		("buffer.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
	;
	
	po::options_description voronOpt("Voronoi options");
	voronOpt.add_options()
		("voronoi.fatInFatSeedDensity",po::value<double>()->default_value(0.0),"fat voronoi seed density (mm^-3)")
		("voronoi.fatInGlandSeedDensity",po::value<double>()->default_value(0.0),"fat voronoi seed in glandular tissue density (mm^-3)")
		("voronoi.glandInGlandSeedDensity",po::value<double>()->default_value(0.0),"glandular voronoi seed density (mm^-3)")
		("voronoi.TDLUDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
		("voronoi.minScaleLenTDLU",po::value<double>()->default_value(0.0),"minimum length scale")
		("voronoi.maxScaleLenTDLU",po::value<double>()->default_value(0.0),"maximum length scale")
		("voronoi.minScaleWidTDLU",po::value<double>()->default_value(0.0),"minimum width scale")
		("voronoi.maxScaleWidTDLU",po::value<double>()->default_value(0.0),"maximum width scale")
		("voronoi.minStrTDLU",po::value<double>()->default_value(0.0),"minimum strength")
		("voronoi.maxStrTDLU",po::value<double>()->default_value(0.0),"maximum strength")
		("voronoi.fatInFatDeflectMax",po::value<double>()->default_value(0.15),"maximum deflection (fraction of pi)")
		("voronoi.minScaleLenFatInFat",po::value<double>()->default_value(0.0),"minimum length scale")
		("voronoi.maxScaleLenFatInFat",po::value<double>()->default_value(0.0),"maximum length scale")
		("voronoi.minScaleWidFatInFat",po::value<double>()->default_value(0.0),"minimum width scale")
		("voronoi.maxScaleWidFatInFat",po::value<double>()->default_value(0.0),"maximum width scale")
		("voronoi.minStrFatInFat",po::value<double>()->default_value(0.0),"minimum strength")
		("voronoi.maxStrFatInFat",po::value<double>()->default_value(0.0),"maximum strength")
		("voronoi.fatInGlandDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
		("voronoi.minScaleLenFatInGland",po::value<double>()->default_value(0.0),"minimum length scale")
		("voronoi.maxScaleLenFatInGland",po::value<double>()->default_value(0.0),"maximum length scale")
		("voronoi.minScaleWidFatInGland",po::value<double>()->default_value(0.0),"minimum width scale")
		("voronoi.maxScaleWidFatInGland",po::value<double>()->default_value(0.0),"maximum width scale")
		("voronoi.minStrFatInGland",po::value<double>()->default_value(0.0),"minimum strength")
		("voronoi.maxStrFatInGland",po::value<double>()->default_value(0.0),"maximum strength")
		("voronoi.glandInGlandDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
		("voronoi.minScaleLenGlandInGland",po::value<double>()->default_value(0.0),"minimum length scale")
		("voronoi.maxScaleLenGlandInGland",po::value<double>()->default_value(0.0),"maximum length scale")
		("voronoi.minScaleWidGlandInGland",po::value<double>()->default_value(0.0),"minimum width scale")
		("voronoi.maxScaleWidGlandInGland",po::value<double>()->default_value(0.0),"maximum width scale")
		("voronoi.minStrGlandInGland",po::value<double>()->default_value(0.0),"minimum strength")
		("voronoi.maxStrGlandInGland",po::value<double>()->default_value(0.0),"maximum strength")
		("voronoi.seedRadius",po::value<double>()->default_value(10),"check seeds in radius (mm)")
	;

	po::options_description fatOpt("Fat lobule options");
	fatOpt.add_options()
		("fat.minLobuleAxis",po::value<double>()->default_value(7.0),"min lobule axis length (mm)")
		("fat.maxLobuleAxis",po::value<double>()->default_value(12.0),"max lobule axis length (mm)")
		("fat.minAxialRatio",po::value<double>()->default_value(0.33),"min axial ratio")
		("fat.maxAxialRatio",po::value<double>()->default_value(0.53),"max axial ratio")
		("fat.minLobuleGap",po::value<double>()->default_value(0.2),"minimum ligament separation between lobules")
		("fat.maxCoeffStr",po::value<double>()->default_value(0.08),"maximum of absolute value of Fourier coefficient as fraction of main radius")
		("fat.minCoeffStr",po::value<double>()->default_value(0.02),"minimum of absolute value of Fourier coefficient as fraction of main radius")
		("fat.maxLobuleTry",po::value<int>()->default_value(600),"maximum number of trial lobules")
	;

	po::options_description ductTreeOpt("Duct tree options");
	ductTreeOpt.add_options()
		("ductTree.maxBranch",po::value<uint>()->default_value(100),"Maximum number of branches")
		("ductTree.maxGen",po::value<uint>()->default_value(15),"Maximum generation")
		("ductTree.initRad",po::value<double>()->default_value(2.0),"tree start radius")
		("ductTree.nFillX",po::value<uint>()->default_value(100),"number x voxels for density map")
		("ductTree.nFillY",po::value<uint>()->default_value(100),"number y voxels for density map")
		("ductTree.nFillZ",po::value<uint>()->default_value(100),"number z voxels for density map")
	;

	po::options_description ductBrOpt("Duct Branch options");
	ductBrOpt.add_options()
		("ductBr.minLen0",po::value<double>()->default_value(2.0),"min length level 0")
		("ductBr.maxLen0",po::value<double>()->default_value(10.0),"max length level 0")
		("ductBr.minLen1",po::value<double>()->default_value(1.0),"min length level 1")
		("ductBr.maxLen1",po::value<double>()->default_value(5.0),"max length level 1")
		("ductBr.minLen2",po::value<double>()->default_value(0.5),"min length level 2")
		("ductBr.maxLen2",po::value<double>()->default_value(4.0),"max length level 2")
		("ductBr.minLenDefault",po::value<double>()->default_value(0.5),"min length default")
		("ductBr.maxLenDefault",po::value<double>()->default_value(1.0),"max length default")
		("ductBr.maxChild",po::value<uint>()->default_value(4),"max number of children")
		("ductBr.childMinRad",po::value<double>()->default_value(0.25),"min radius to have children")
		("ductBr.childLevBound",po::value<uint>()->default_value(4),"max level for child probability")
		("ductBr.child00",po::value<double>()->default_value(0.0),"cumulative probability 0->0")
		("ductBr.child01",po::value<double>()->default_value(0.0),"cumulative probability 0->1")
		("ductBr.child02",po::value<double>()->default_value(0.5),"cumulative probability 0->2")
		("ductBr.child10",po::value<double>()->default_value(0.0),"cumulative probability 1->0")
		("ductBr.child11",po::value<double>()->default_value(0.1),"cumulative probability 1->1")
		("ductBr.child12",po::value<double>()->default_value(0.75),"cumulative probability 1->2")
		("ductBr.child20",po::value<double>()->default_value(0.1),"cumulative probability 2->0")
		("ductBr.child21",po::value<double>()->default_value(0.5),"cumulative probability 2->1")
		("ductBr.child22",po::value<double>()->default_value(0.9),"cumulative probability 2->2")
		("ductBr.child30",po::value<double>()->default_value(0.2),"cumulative probability 3->0")
		("ductBr.child31",po::value<double>()->default_value(0.7),"cumulative probability 3->1")
		("ductBr.child32",po::value<double>()->default_value(0.95),"cumulative probability 3->2")
		("ductBr.child40",po::value<double>()->default_value(0.5),"cumulative probability 3->0")
		("ductBr.child41",po::value<double>()->default_value(0.9),"cumulative probability 3->1")
		("ductBr.child42",po::value<double>()->default_value(0.99),"cumulative probability 3->2")
		("ductBr.minRadFrac",po::value<double>()->default_value(0.6),"min start radius as fraction of parent radius")
		("ductBr.maxRadFrac",po::value<double>()->default_value(0.85),"max start radius as fraction of parent radius")
		("ductBr.radFrac0",po::value<double>()->default_value(0.9),"start radius as fraction of parent radius for level 0")
		("ductBr.minAngle",po::value<double>()->default_value(0.1),"min angle between parent and child")
		("ductBr.maxAngle",po::value<double>()->default_value(0.5),"max angle between parent and child")
	;

	po::options_description ductSegOpt("Duct Segment options");
	ductSegOpt.add_options()
		("ductSeg.lengthBetaA",po::value<double>()->default_value(2.0),"length distribution shape parameter A")
		("ductSeg.lengthBetaB",po::value<double>()->default_value(2.0),"length distribution shape parameter B")
		("ductSeg.radiusBetaA",po::value<double>()->default_value(2.0),"radius distribution shape parameter A")
		("ductSeg.radiusBetaB",po::value<double>()->default_value(2.0),"radius distribution shape parameter B")
		("ductSeg.minLen",po::value<double>()->default_value(1.0),"min segment length")
		("ductSeg.maxLen",po::value<double>()->default_value(5.0),"max segment length")
		("ductSeg.maxCurvRad",po::value<double>()->default_value(20.0),"max radius of curvature")
		("ductSeg.maxCurvFrac",po::value<double>()->default_value(0.33),"max fraction of circle")
		("ductSeg.minEndRad",po::value<double>()->default_value(0.85),"min end radius as fraction of start radius")
		("ductSeg.maxEndRad",po::value<double>()->default_value(1.05),"max end radius as fraction of start radius")
		("ductSeg.angleWt",po::value<double>()->default_value(1.0),"cost function preferential angle weighting")
		("ductSeg.densityWt",po::value<double>()->default_value(5e-5),"cost function density weighting")
		("ductSeg.numTry",po::value<uint>()->default_value(10),"number of trial segments")
		("ductSeg.maxTry",po::value<uint>()->default_value(100),"max number of trial segments before reducing length")
		("ductSeg.absMaxTry",po::value<uint>()->default_value(10000),"max number of trial segments before completely giving up")
		// check if this needs to be changed
		("ductSeg.roiStep",po::value<double>()->default_value(0.1),"step size for checking segment validity")
	;

	po::options_description vesselTreeOpt("Vessel tree options");
	vesselTreeOpt.add_options()
		("vesselTree.maxBranch",po::value<uint>()->default_value(100),"Maximum number of branches")
		("vesselTree.maxGen",po::value<uint>()->default_value(15),"Maximum generation")
		("vesselTree.initRad",po::value<double>()->default_value(2.0),"tree start radius")
		("vesselTree.nFillX",po::value<uint>()->default_value(100),"number x voxels for density map")
		("vesselTree.nFillY",po::value<uint>()->default_value(100),"number y voxels for density map")
		("vesselTree.nFillZ",po::value<uint>()->default_value(100),"number z voxels for density map")
	;

	po::options_description vesselBrOpt("Vessel branch options");
	vesselBrOpt.add_options()
		("vesselBr.minLen0",po::value<double>()->default_value(2.0),"min length level 0")
		("vesselBr.maxLen0",po::value<double>()->default_value(10.0),"max length level 0")
		("vesselBr.minLen1",po::value<double>()->default_value(1.0),"min length level 1")
		("vesselBr.maxLen1",po::value<double>()->default_value(5.0),"max length level 1")
		("vesselBr.minLen2",po::value<double>()->default_value(0.5),"min length level 2")
		("vesselBr.maxLen2",po::value<double>()->default_value(4.0),"max length level 2")
		("vesselBr.minLenDefault",po::value<double>()->default_value(0.5),"min length default")
		("vesselBr.maxLenDefault",po::value<double>()->default_value(1.0),"max length default")
		("vesselBr.maxChild",po::value<uint>()->default_value(4),"max number of children")
		("vesselBr.childMinRad",po::value<double>()->default_value(0.25),"min radius to have children")
		("vesselBr.childLevBound",po::value<uint>()->default_value(4),"max level for child probability")
		("vesselBr.child00",po::value<double>()->default_value(0.0),"cumulative probability 0->0")
		("vesselBr.child01",po::value<double>()->default_value(0.0),"cumulative probability 0->1")
		("vesselBr.child02",po::value<double>()->default_value(0.5),"cumulative probability 0->2")
		("vesselBr.child10",po::value<double>()->default_value(0.0),"cumulative probability 1->0")
		("vesselBr.child11",po::value<double>()->default_value(0.1),"cumulative probability 1->1")
		("vesselBr.child12",po::value<double>()->default_value(0.75),"cumulative probability 1->2")
		("vesselBr.child20",po::value<double>()->default_value(0.1),"cumulative probability 2->0")
		("vesselBr.child21",po::value<double>()->default_value(0.5),"cumulative probability 2->1")
		("vesselBr.child22",po::value<double>()->default_value(0.9),"cumulative probability 2->2")
		("vesselBr.child30",po::value<double>()->default_value(0.2),"cumulative probability 3->0")
		("vesselBr.child31",po::value<double>()->default_value(0.7),"cumulative probability 3->1")
		("vesselBr.child32",po::value<double>()->default_value(0.95),"cumulative probability 3->2")
		("vesselBr.child40",po::value<double>()->default_value(0.5),"cumulative probability 3->0")
		("vesselBr.child41",po::value<double>()->default_value(0.9),"cumulative probability 3->1")
		("vesselBr.child42",po::value<double>()->default_value(0.99),"cumulative probability 3->2")
		("vesselBr.minRadFrac",po::value<double>()->default_value(0.6),"min start radius as fraction of parent radius")
		("vesselBr.maxRadFrac",po::value<double>()->default_value(0.85),"max start radius as fraction of parent radius")
		("vesselBr.radFrac0",po::value<double>()->default_value(0.9),"start radius as fraction of parent radius for level 0")
		("vesselBr.minAngle",po::value<double>()->default_value(0.1),"min angle between parent and child")
		("vesselBr.maxAngle",po::value<double>()->default_value(0.5),"max angle between parent and child")
	;

	po::options_description vesselSegOpt("Vessel segment options");
	vesselSegOpt.add_options()
		("vesselSeg.lengthBetaA",po::value<double>()->default_value(2.0),"length distribution shape parameter A")
		("vesselSeg.lengthBetaB",po::value<double>()->default_value(2.0),"length distribution shape parameter B")
		("vesselSeg.radiusBetaA",po::value<double>()->default_value(2.0),"radius distribution shape parameter A")
		("vesselSeg.radiusBetaB",po::value<double>()->default_value(2.0),"radius distribution shape parameter B")
		("vesselSeg.minLen",po::value<double>()->default_value(1.0),"min segment length")
		("vesselSeg.maxLen",po::value<double>()->default_value(5.0),"max segment length")
		("vesselSeg.maxCurvRad",po::value<double>()->default_value(20.0),"max radius of curvature")
		("vesselSeg.maxCurvFrac",po::value<double>()->default_value(0.33),"max fraction of circle")
		("vesselSeg.minEndRad",po::value<double>()->default_value(0.85),"min end radius as fraction of start radius")
		("vesselSeg.maxEndRad",po::value<double>()->default_value(1.05),"max end radius as fraction of start radius")
		("vesselSeg.angleWt",po::value<double>()->default_value(1.0),"cost function preferential angle weighting")
		("vesselSeg.densityWt",po::value<double>()->default_value(5e-5),"cost function density weighting")
		("vesselSeg.numTry",po::value<uint>()->default_value(10),"number of trial segments")
		("vesselSeg.maxTry",po::value<uint>()->default_value(100),"max number of trial segments before reducing length")
		("vesselSeg.absMaxTry",po::value<uint>()->default_value(10000),"max number of trial segments before completely giving up")
		// check if this needs to be changed
		("vesselSeg.roiStep",po::value<double>()->default_value(0.1),"step size for checking segment validity")
	;

	// config file options
	po::options_description configFileOpt("Configuration file options");
	configFileOpt.add(baseOpt).add(shapeOpt);
	configFileOpt.add(ductTreeOpt).add(ductBrOpt).add(ductSegOpt);
	configFileOpt.add(vesselTreeOpt).add(vesselBrOpt).add(vesselSegOpt);
	configFileOpt.add(compartOpt).add(TDLUOpt).add(fatOpt);
	configFileOpt.add(voronOpt).add(perlinOpt).add(boundaryOpt);
	configFileOpt.add(bufferOpt).add(perturbOpt);

	// all of the options
	po::options_description all("All options");
	all.add_options()
		("config,c", po::value<std::string>()->required(), "name of configuration file")
	;
	all.add(configFileOpt);

	po::variables_map vm;

	// get configuration filename from command line
	po::store(parse_command_line(argc,argv,all), vm);
	std::string configFile = vm["config"].as<std::string>();

	// read configuration file
	ifstream inConfig(configFile.c_str());
	if(!inConfig){
		cout << "Can not open configuration file: " << configFile << "\n";
		return(1);
	} else {
		po::store(parse_config_file(inConfig, configFileOpt), vm);
		inConfig.close();
	};

	// load resolution variables

	double ures = vm["shape.ures"].as<double>();	// spacing in u-v space
	double vres = vm["shape.vres"].as<double>();
	double pointSep = vm["shape.pointSep"].as<double>(); // remove close points
	double featureAngle = vm["shape.featureAngle"].as<double>(); // angle to preserve when smoothing
	double targetReduction = vm["shape.targetReduction"].as<double>(); // fraction of triangles to decimate
	double imgRes = vm["base.imgRes"].as<double>();	// size of image voxels (mm)
	double skinThick = vm["base.skinThick"].as<double>(); // thickness of skin (mm)
	std::string baseShapeFilename = "baseShape.vtp";
	std::string breastVoxelFilename = "breastVoxel.vti";
	std::string breastPointsFilename = "breastPoints.vtp";
	std::string segmentedVoxelFilename = "segmentedVoxel.vti";
	std::string ductVoxelFilename = "final.vti";
	std::string skinMeshFilename = "skinMesh.vtp";

	double scaleFactor = 50.0;	// scale voxel size to millimeters

	double nippleLen = vm["base.nippleLen"].as<double>();	// length of nipple (mm)
	double nippleRad = vm["base.nippleRad"].as<double>(); // radius of nipple (mm)
	double areolaRad = vm["base.areolaRad"].as<double>(); // radius of areola (mm)

	double voxelVol = imgRes*imgRes*imgRes;	// volume of a voxel (cubic mm)
	
	tissueStruct tissue;
	
	tissue.bg = 0;
	tissue.skin = 2;
	tissue.nipple = 33;
	tissue.fat = 1;
	tissue.cooper = 88;
	tissue.gland = 29;
	tissue.TDLU = 95;
	tissue.duct = 125;
	tissue.artery = 150;
	tissue.vein = 225 ;
	tissue.muscle = 40;

	const unsigned char glandBuffer = 222;
	const unsigned char innerVal = 175; // not needed after compartment creation
	const unsigned char boundVal = 50;  // edge of inner breast volume, not needed after skin generation

	const unsigned char compartmentVal[20] = {12, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};

	const unsigned char compMax = 28;	// max compartment value for cooper's ligament part
	const unsigned char compMin = 12;


	// random number generator seed
	// use seed specified in config file, otherwise random seed
	int randSeed;
	if(vm.count("base.seed")){
		// seed specified in configuration file
		randSeed = vm["base.seed"].as<unsigned int>();
	} else {
		// generate random seed
		FILE *randFile;
		randFile = fopen("/dev/urandom","r");
		fread((char*)(&randSeed), sizeof(int),1,randFile);
		fclose(randFile);
	}

	// output base directory
	std::string outputDir = vm["base.outputDir"].as<std::string>();

	// exists?
	if(access(outputDir.c_str(),F_OK)){
		// directory doesn't exist
		// try to create it
		if(mkdir(outputDir.c_str(),S_IRWXU) == -1){
			// couldn't create directory
			cerr << "Could not create directory " << outputDir << "\n";
			cerr << "Exiting...\n";
			return(1);
		}
	}

	// is it a directory?
	struct stat fileInfo;

	stat(outputDir.c_str(), &fileInfo);

	if(S_ISDIR(fileInfo.st_mode)){
		// it's a directory
		// make it the working directory
		if(chdir(outputDir.c_str()) == -1){
			cerr << "Could not access directory " << outputDir << "\n";
			cerr << "Exiting...\n";
			return(1);
		}
	} else {
		// error, not a directory
		cerr << "Specified path is not a valid directory\n";
		cerr << "Exiting...\n";
		return(1);
	}


	// copy over config file
	FILE *cfgCopy = fopen("my.cfg","wt");
	FILE *cfgRead = fopen(configFile.c_str(),"rt");
	
	int readChar = getc(cfgRead);
	while(readChar != EOF){
		putc(readChar, cfgCopy);
		readChar = getc(cfgRead);
	}
	fclose(cfgCopy);
	fclose(cfgRead);


	// shape parameters

	// base shape coefficients
	double a1b =  vm["shape.a1b"].as<double>();
	double a1t = vm["shape.a1t"].as<double>();
	double a2l = vm["shape.a2l"].as<double>();
	double a2r = vm["shape.a2r"].as<double>();
	double a3 = vm["shape.a3"].as<double>();
	double eps1 = vm["shape.eps1"].as<double>();
	double eps2 = vm["shape.eps2"].as<double>();

	// breast side
	bool leftSide = vm["base.leftBreast"].as<bool>();

	// ptosis parameters
	bool doPtosis = vm["shape.doPtosis"].as<bool>();
	double b0 = vm["shape.ptosisB0"].as<double>();
	double b1 = vm["shape.ptosisB1"].as<double>();

	// turn parameters
	bool doTurn = vm["shape.doTurn"].as<bool>();
	double c0 = vm["shape.turnC0"].as<double>();
	double c1 = vm["shape.turnC1"].as<double>();

	// top shape parameters
	bool doTopShape = vm["shape.doTopShape"].as<bool>();
	double s0 = vm["shape.topShapeS0"].as<double>();
	double t0 = vm["shape.topShapeT0"].as<double>();
	double s1 = vm["shape.topShapeS1"].as<double>();
	double t1 = vm["shape.topShapeT1"].as<double>();

	// derived top shape parameters
	double At = -0.5*t0-3.0*s0-3.0*s1+0.5*t1;
	double Bt = 1.5*t0+8.0*s0+7.0*s1-t1;
	double Ct = -1.5*t0-6.0*s0-4.0*s1+0.5*t1;
	double Dt = 0.5*t0;
	double Et = s0;
	double Ft = 1.0;

	// flatten side parameters
	bool doFlattenSide = vm["shape.doFlattenSide"].as<bool>();
	double g0 = vm["shape.flattenSideG0"].as<double>();
	double g1 = vm["shape.flattenSideG1"].as<double>();

	// derived parameters for flatten side
	double Af = g1+2.0-2.0*g0;
	double Bf = -g1-3.0+3.0*g0;
	double Cf = 0.0;
	double Df = 1.0;

	// turn top parameters
	bool doTurnTop = vm["shape.doTurnTop"].as<bool>();
	double h0 = vm["shape.turnTopH0"].as<double>();
	double h1 = vm["shape.turnTopH1"].as<double>();

	// start a random number generator
	vtkSmartPointer<vtkMinimalStandardRandomSequence> rgen =
		vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();

	rgen->SetSeed(randSeed);

	/***********************
	Shape
	***********************/

	// timing
	time(&currentTime);
	//cout << "\nTime interval 1 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 

	// create base shape
	cout << "Initializing point cloud...";
	// point positions
	double uval,vval,xval,yval,zval;
	vtkIdType myId;

	vtkSmartPointer<vtkPoints> brightFront =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> brightBack =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> brightRing =
		vtkSmartPointer<vtkPoints>::New();
		
	vtkSmartPointer<vtkPoints> trightFront =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> trightBack =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> trightRing =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkPoints> tleftFront =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> tleftBack =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> tleftRing =
		vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkPoints> bleftFront =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> bleftBack =
		vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> bleftRing =
		vtkSmartPointer<vtkPoints>::New();
	
	vtkIdType centerPt;
	
	#pragma omp sections private(uval,vval,xval,yval,zval,myId)
	{
		#pragma omp section
		{
		// bottom right breast
		vval = -1.0*pi + vres;
		while(vval <= -0.5*pi){
			uval = 0.0;
			// ring is uval == 0
			xval = 0.0;
			yval = pow(a2r,eps1)*pow(sin(vval),eps2);
			zval = pow(a1b,eps1)*pow(cos(vval),eps2);
			myId = brightRing->InsertNextPoint(xval,yval,zval);
			uval += ures;
			
			while(uval <= 0.5*pi){
				xval = pow(a3*sin(uval),eps1);
				yval = pow(a2r*cos(uval),eps1)*pow(sin(vval),eps2);
				zval = pow(a1b*cos(uval),eps1)*pow(cos(vval),eps2);
				// add point to front
				myId = brightFront->InsertNextPoint(xval,yval,zval);
				// add backplane point
				myId = brightBack->InsertNextPoint(0.0,yval,zval);
				uval += ures;
			}
			vval += vres;
		}
	
		// add center point to bottom-right breast front and back and skip
		centerPt = brightFront->InsertNextPoint(pow(a3,eps1),0.0,0.0);
		centerPt = brightBack->InsertNextPoint(0.0,0.0,0.0);
		}

		#pragma omp section
		{
		// top right breast
	
		vval = -0.5*pi + vres;
		while(vval <= 0.0){
			uval = 0.0;
			xval = 0.0;
			yval = pow(a2r,eps1)*pow(sin(vval),eps2);
			zval = pow(a1t,eps1)*pow(cos(vval),eps2);
			myId = trightRing->InsertNextPoint(xval,yval,zval);
			uval += ures;
			
			while(uval <= 0.5*pi){
				xval = pow(a3*sin(uval),eps1);
				yval = pow(a2r*cos(uval),eps1)*pow(sin(vval),eps2);
				zval = pow(a1t*cos(uval),eps1)*pow(cos(vval),eps2);
				// add the point
				myId = trightFront->InsertNextPoint(xval,yval,zval);
				// add backplane point
				myId = trightBack->InsertNextPoint(0.0,yval,zval);
				uval += ures;
			}
			vval += vres;
		}
		}
		
		#pragma omp section
		{
		// top left breast
		
		vval = 0.0 + vres;
		while(vval <= 0.5*pi){
			uval = 0.0;
			xval = 0.0;
			yval = pow(a2l,eps1)*pow(sin(vval),eps2);
			zval = pow(a1t,eps1)*pow(cos(vval),eps2);
			myId = tleftRing->InsertNextPoint(xval,yval,zval);
			uval += ures;
			
			while(uval <= 0.5*pi){
				xval = pow(a3*sin(uval),eps1);
				yval = pow(a2l*cos(uval),eps1)*pow(sin(vval),eps2);
				zval = pow(a1t*cos(uval),eps1)*pow(cos(vval),eps2);
				// add the point
				myId = tleftFront->InsertNextPoint(xval,yval,zval);
				// add backplane point
				myId = tleftBack->InsertNextPoint(0.0,yval,zval);
				uval += ures;
			}
			vval += vres;
		}
		}

		#pragma omp section
		{
		// bottom left breast
	
		vval = 0.5*pi + vres;
		while(vval <= pi){
			uval = 0.0;
			xval = 0.0;
			yval = pow(a2l,eps1)*pow(sin(vval),eps2);
			zval = pow(a1b,eps1)*pow(cos(vval),eps2);
			myId = bleftRing->InsertNextPoint(xval,yval,zval);
			uval += ures;
			
			while(uval <= 0.5*pi){
				xval = pow(a3*sin(uval),eps1);
				yval = pow(a2l*cos(uval),eps1)*pow(sin(vval),eps2);
				zval = pow(a1b*cos(uval),eps1)*pow(cos(vval),eps2);
				// add the point
				myId = bleftFront->InsertNextPoint(xval,yval,zval);
				// add backplane point
				myId = bleftBack->InsertNextPoint(0.0,yval,zval);
				uval += ures;
			}
			vval += vres;
		}
		}

	}  // end omp sections
	
	cout << "done.\n";

	// have base shape
	// do deformations

	cout << "Doing deformation...";

	// timing
	time(&currentTime);
	//cout << "\nTime interval 2 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 

	// top shape
	if(doTopShape){
		vtkIdType npts = trightFront->GetNumberOfPoints();
		uval = ures; // need original u value
		for(int i=0; i<npts; i++){
			double pval[3];
			double u2 = uval*2.0/pi;
			// front
			trightFront->GetPoint(i, pval);
			pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
				Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
			trightFront->SetPoint(i,pval);
			// back
			trightBack->GetPoint(i, pval);
			pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
				Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
			trightBack->SetPoint(i,pval);
			
			uval += ures;
			if(uval > 0.5*pi){
				uval = ures;
			}
		}
		
		// ring
		npts = trightRing->GetNumberOfPoints();
		for(int i=0; i<npts; i++){
			double pval[3];
			trightRing->GetPoint(i, pval);
			pval[2] = pval[2]*Ft;
			trightRing->SetPoint(i,pval);
		}
		
		npts = tleftFront->GetNumberOfPoints();
		uval = ures; // need original u value
		for(int i=0; i<npts; i++){
			double pval[3];
			double u2 = uval*2.0/pi;
			// front
			tleftFront->GetPoint(i, pval);
			pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
				Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
			tleftFront->SetPoint(i,pval);
			// back
			tleftBack->GetPoint(i, pval);
			pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
				Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
			tleftBack->SetPoint(i,pval);
			
			uval += ures;
			if(uval > 0.5*pi){
				uval = ures;
			}
		}
		
		// ring
		npts = tleftRing->GetNumberOfPoints();
		for(int i=0; i<npts; i++){
			double pval[3];
			tleftRing->GetPoint(i, pval);
			pval[2] = pval[2]*Ft;
			tleftRing->SetPoint(i,pval);
		}
		
		
	}
	
	// flatten side
	if(doFlattenSide){
		if(leftSide){
			trightFront->ComputeBounds();
			double bound[6];
			trightFront->GetBounds(bound);
			double scale;
			if(fabs(bound[2]) > fabs(bound[3])){
				scale = fabs(bound[2]);
			} else {
				scale = fabs(bound[3]);
			}
			vtkIdType npts = trightFront->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				// front
				trightFront->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				trightFront->SetPoint(i,pval);
				// back
				trightBack->GetPoint(i,pval);
				yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				trightBack->SetPoint(i,pval);
			}
			// ring
			npts = trightRing->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				trightRing->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				trightRing->SetPoint(i,pval);
			}
			
			
			brightFront->ComputeBounds();
			brightFront->GetBounds(bound);
			if(fabs(bound[2]) > fabs(bound[3])){
				scale = fabs(bound[2]);
			} else {
				scale = fabs(bound[3]);
			}
			npts = brightFront->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				//front
				brightFront->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				brightFront->SetPoint(i,pval);
				// back
				brightBack->GetPoint(i,pval);
				yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				brightBack->SetPoint(i,pval);
			}
			// ring
			npts = brightRing->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				brightRing->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				brightRing->SetPoint(i,pval);
			}
		} else {
			// right side
			tleftFront->ComputeBounds();
			double bound[6];
			tleftFront->GetBounds(bound);
			double scale;
			if(fabs(bound[2]) > fabs(bound[3])){
				scale = fabs(bound[2]);
			} else {
				scale = fabs(bound[3]);
			}
			vtkIdType npts = tleftFront->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				// front
				tleftFront->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				tleftFront->SetPoint(i,pval);
				// back
				tleftBack->GetPoint(i,pval);
				yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				tleftBack->SetPoint(i,pval);
			}
			// ring
			npts = tleftRing->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				tleftRing->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				tleftRing->SetPoint(i,pval);
			}
			
			bleftFront->ComputeBounds();
			bleftFront->GetBounds(bound);
			if(fabs(bound[2]) > fabs(bound[3])){
				scale = fabs(bound[2]);
			} else {
				scale = fabs(bound[3]);
			}
			npts = bleftFront->GetNumberOfPoints();
			#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				// front
				bleftFront->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				bleftFront->SetPoint(i,pval);
				// back
				bleftBack->GetPoint(i,pval);
				yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				bleftBack->SetPoint(i,pval);
			}
			// ring
			npts = bleftRing->GetNumberOfPoints();
			//#pragma omp parallel for
			for(int i=0; i<npts; i++){
				double pval[3];
				bleftRing->GetPoint(i,pval);
				double yv = fabs(pval[1]/scale);
				pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
					Cf*yv + Df);
				bleftRing->SetPoint(i,pval);
			}
		}
	}
	
	// turn top
	if(doTurnTop){
		trightFront->ComputeBounds();
		double bound[6];
		trightFront->GetBounds(bound);
		double scale;
		if(fabs(bound[4]) > fabs(bound[5])){
			scale = fabs(bound[4]);
		} else {
			scale = fabs(bound[5]);
		}
		vtkIdType npts = trightFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			trightFront->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			trightFront->SetPoint(i,pval);
			// back
			trightBack->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			trightBack->SetPoint(i,pval);
		}
		// ring
		npts = trightRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			trightRing->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			trightRing->SetPoint(i,pval);
		}
		
		tleftFront->ComputeBounds();
		tleftFront->GetBounds(bound);
		if(fabs(bound[4]) > fabs(bound[5])){
			scale = fabs(bound[4]);
		} else {
			scale = fabs(bound[5]);
		}
		npts = tleftFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			tleftFront->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			tleftFront->SetPoint(i,pval);
			// back
			tleftBack->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			tleftBack->SetPoint(i,pval);
		}
		// ring
		npts = tleftRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			tleftRing->GetPoint(i,pval);
			pval[1] = pval[1] - h0*pval[2]/scale -
				h1*pval[2]*pval[2]/scale/scale;
			tleftRing->SetPoint(i,pval);
		}
	}

	// ptosis
	if(doPtosis){
		vtkIdType npts = trightFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			trightFront->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			trightFront->SetPoint(i,pval);
			// back
			trightBack->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			trightBack->SetPoint(i,pval);
		}
		// ring
		npts = trightRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			trightRing->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			trightRing->SetPoint(i,pval);
		}
		
		npts = tleftFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			tleftFront->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			tleftFront->SetPoint(i,pval);
			// back
			tleftBack->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			tleftBack->SetPoint(i,pval);
		}
		// ring
		npts = tleftRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			tleftRing->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			tleftRing->SetPoint(i,pval);
		}
		
		npts = brightFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			brightFront->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			brightFront->SetPoint(i,pval);
			// back
			brightBack->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			brightBack->SetPoint(i,pval);
		}
		// ring
		npts = brightRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			brightRing->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			brightRing->SetPoint(i,pval);
		}
		
		npts = bleftFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			bleftFront->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			bleftFront->SetPoint(i,pval);
			// back
			bleftBack->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			bleftBack->SetPoint(i,pval);
		}
		// ring
		npts = bleftRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			bleftRing->GetPoint(i,pval);
			pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
			bleftRing->SetPoint(i,pval);
		}
	}

	// turn
	if(doTurn){
		vtkIdType npts = trightFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			trightFront->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			trightFront->SetPoint(i,pval);
			// back
			trightBack->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			trightBack->SetPoint(i,pval);
		}
		// ring
		npts = trightRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			trightRing->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			trightRing->SetPoint(i,pval);
		}
		
		npts = tleftFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			tleftFront->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			tleftFront->SetPoint(i,pval);
			// back
			tleftBack->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			tleftBack->SetPoint(i,pval);
		}
		// ring
		npts = tleftRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			tleftRing->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			tleftRing->SetPoint(i,pval);
		}
		
		npts = brightFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			brightFront->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			brightFront->SetPoint(i,pval);
			// back
			brightBack->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			brightBack->SetPoint(i,pval);
		}
		// ring
		npts = brightRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			brightRing->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			brightRing->SetPoint(i,pval);
		}
		
		npts = bleftFront->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			// front
			bleftFront->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			bleftFront->SetPoint(i,pval);
			// back
			bleftBack->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			bleftBack->SetPoint(i,pval);
		}
		// ring
		npts = bleftRing->GetNumberOfPoints();
		#pragma omp parallel for
		for(int i=0; i<npts; i++){
			double pval[3];
			bleftRing->GetPoint(i,pval);
			pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
			bleftRing->SetPoint(i,pval);
		}
	}

	// timing
	time(&currentTime);
	//cout << "\nTime interval 3 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 


	// get nipple position
	double nipplePos[3];
	brightFront->GetPoint(centerPt,nipplePos);

	// scale nipple position
	for(int i=0; i<3; i++){
		nipplePos[i] = nipplePos[i]*scaleFactor;
	}

	// aggregate front points
	vtkSmartPointer<vtkPoints> frontPts =
		vtkSmartPointer<vtkPoints>::New();
	
	for(vtkIdType i=0; i<brightFront->GetNumberOfPoints(); i++){
		double t[3];
		brightFront->GetPoint(i, t);
		frontPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<trightFront->GetNumberOfPoints(); i++){
		double t[3];
		trightFront->GetPoint(i, t);
		frontPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<bleftFront->GetNumberOfPoints(); i++){
		double t[3];
		bleftFront->GetPoint(i, t);
		frontPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<tleftFront->GetNumberOfPoints(); i++){
		double t[3];
		tleftFront->GetPoint(i, t);
		frontPts->InsertNextPoint(t);
	}
	vtkIdType numFrontPts = frontPts->GetNumberOfPoints();
	
	// decimate front points
	vtkSmartPointer<vtkPolyData> frontPoly = 
		vtkSmartPointer<vtkPolyData>::New();
	frontPoly->SetPoints(frontPts);
	vtkSmartPointer<vtkVertexGlyphFilter> frontVertAdd =
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	frontVertAdd->SetInput(frontPoly);
#else
	frontVertAdd->SetInputData(frontPoly);
#endif


	frontVertAdd->Update();
	vtkSmartPointer<vtkCleanPolyData> cleanFront =
		vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	cleanFront->SetInput(frontVertAdd->GetOutput());
#else
	cleanFront->SetInputConnection(frontVertAdd->GetOutputPort());
#endif
	cleanFront->SetTolerance(pointSep);
	cleanFront->Update();	
	
	// aggregate and shift back points
	vtkSmartPointer<vtkPoints> backPts =
		vtkSmartPointer<vtkPoints>::New();
	
	// ringWidth and ringSep converted to non-physical coordinates
	double ringWidth = vm["shape.ringWidth"].as<double>();
	double ringSep = vm["shape.ringSep"].as<double>();
	double ringWidthOrig = ringWidth/scaleFactor;
	double ringSepOrig = ringSep/scaleFactor;
	double backPos = (floor(ringWidthOrig/ringSepOrig)+1)*ringSepOrig;
	
	for(vtkIdType i=0; i<brightBack->GetNumberOfPoints(); i++){
		double t[3];
		brightBack->GetPoint(i, t);
		t[0] = -backPos;
		backPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<trightBack->GetNumberOfPoints(); i++){
		double t[3];
		trightBack->GetPoint(i, t);
		t[0] = -backPos;
		backPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<bleftBack->GetNumberOfPoints(); i++){
		double t[3];
		bleftBack->GetPoint(i, t);
		t[0] = -backPos;
		backPts->InsertNextPoint(t);
	}
	for(vtkIdType i=0; i<tleftBack->GetNumberOfPoints(); i++){
		double t[3];
		tleftBack->GetPoint(i, t);
		t[0] = -backPos;
		backPts->InsertNextPoint(t);
	}
	vtkIdType numBackPts = backPts->GetNumberOfPoints();
	
	// decimate back points
	vtkSmartPointer<vtkPolyData> backPoly = 
		vtkSmartPointer<vtkPolyData>::New();
	backPoly->SetPoints(backPts);
	vtkSmartPointer<vtkVertexGlyphFilter> backVertAdd =
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	backVertAdd->SetInput(backPoly);
#else
	backVertAdd->SetInputData(backPoly);
#endif
	backVertAdd->Update();
	vtkSmartPointer<vtkCleanPolyData> cleanBack =
		vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	cleanBack->SetInput(backVertAdd->GetOutput());
#else
	cleanBack->SetInputConnection(backVertAdd->GetOutputPort());
#endif
	cleanBack->SetTolerance(pointSep);
	cleanBack->Update();	
	
	// aggregate ring points and create thickness
	vtkSmartPointer<vtkPoints> ringPts =
		vtkSmartPointer<vtkPoints>::New();
	
	for(vtkIdType i=0; i<brightRing->GetNumberOfPoints(); i++){
		double t[3];
		brightRing->GetPoint(i, t);
		ringPts->InsertNextPoint(t);
		double xstep = ringSepOrig;
		while(xstep < ringWidthOrig){
			t[0] = -xstep;
			ringPts->InsertNextPoint(t);
			xstep += ringSepOrig;
		}
	}
	for(vtkIdType i=0; i<trightRing->GetNumberOfPoints(); i++){
		double t[3];
		trightRing->GetPoint(i, t);
		ringPts->InsertNextPoint(t);
		double xstep = ringSepOrig;
		while(xstep < ringWidthOrig){
			t[0] = -xstep;
			ringPts->InsertNextPoint(t);
			xstep += ringSepOrig;
		}
	}
	for(vtkIdType i=0; i<bleftRing->GetNumberOfPoints(); i++){
		double t[3];
		bleftRing->GetPoint(i, t);
		ringPts->InsertNextPoint(t);
		double xstep = ringSepOrig;
		while(xstep < ringWidthOrig){
			t[0] = -xstep;
			ringPts->InsertNextPoint(t);
			xstep += ringSepOrig;
		}
	}
	for(vtkIdType i=0; i<tleftRing->GetNumberOfPoints(); i++){
		double t[3];
		tleftRing->GetPoint(i, t);
		ringPts->InsertNextPoint(t);
		double xstep = ringSepOrig;
		while(xstep < ringWidthOrig){
			t[0] = -xstep;
			ringPts->InsertNextPoint(t);
			xstep += ringSepOrig;
		}
	}
	vtkIdType numRingPts = ringPts->GetNumberOfPoints();
	
	// timing
	time(&currentTime);
	//cout << "\nTime interval 4 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 
	
	// decimate ring points
	vtkSmartPointer<vtkPolyData> ringPoly = 
		vtkSmartPointer<vtkPolyData>::New();
	ringPoly->SetPoints(ringPts);
	vtkSmartPointer<vtkVertexGlyphFilter> ringVertAdd =
		vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	ringVertAdd->SetInput(ringPoly);
#else
	ringVertAdd->SetInputData(ringPoly);
#endif
	ringVertAdd->Update();
	vtkSmartPointer<vtkCleanPolyData> cleanRing =
		vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	cleanRing->SetInput(ringVertAdd->GetOutput());
#else
	cleanRing->SetInputConnection(ringVertAdd->GetOutputPort());
#endif
	cleanRing->SetTolerance(pointSep);
	cleanRing->Update();	
	
	// timing
	time(&currentTime);
	//cout << "\nTime interval 5 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 
	
	// allocate space for breast points
	vtkSmartPointer<vtkPoints> breastPts =
		vtkSmartPointer<vtkPoints>::New();

	vtkIdType numBreastPts = cleanFront->GetOutput()->GetNumberOfPoints();
	numBreastPts += cleanBack->GetOutput()->GetNumberOfPoints();
	numBreastPts += cleanRing->GetOutput()->GetNumberOfPoints();
	breastPts->SetNumberOfPoints(numBreastPts);
	
	// cell array for verticies
	vtkSmartPointer<vtkCellArray> breastVerts =
		vtkSmartPointer<vtkCellArray>::New();
	
	// scale to physical units
	vtkIdType totalCount = 0;
	vtkIdType currentPoints = cleanFront->GetOutput()->GetNumberOfPoints();

	for(vtkIdType i=0; i<currentPoints; i++){
		double t[3];
		
		cleanFront->GetOutput()->GetPoint(i,t);
		for(int j=0; j<3; j++){
			t[j] = t[j]*scaleFactor;
		}
		breastPts->SetPoint(totalCount,t);
		breastVerts->InsertNextCell(1,&totalCount); 
		totalCount++;
	}

	// calculate scaled back points
	currentPoints = cleanBack->GetOutput()->GetNumberOfPoints();
	for(vtkIdType i=0; i<currentPoints; i++){
		double t[3];
		cleanBack->GetOutput()->GetPoint(i,t);
		for(int j=0; j<3; j++){
			t[j] = t[j]*scaleFactor;
		}
		breastPts->SetPoint(totalCount,t);
		breastVerts->InsertNextCell(1,&totalCount); 
		totalCount++;
	}
	
	// calculate scaled ring points and normals
	currentPoints = cleanRing->GetOutput()->GetNumberOfPoints();

	for(vtkIdType i=0; i<currentPoints; i++){
		double t[3];
		cleanRing->GetOutput()->GetPoint(i,t);
		for(int j=0; j<3; j++){
			t[j] = t[j]*scaleFactor;
		}
		breastPts->SetPoint(totalCount,t);
		breastVerts->InsertNextCell(1,&totalCount); 
		totalCount++;
	}

	cout << "done.\n";

	cout << "Meshing base shape...";

	// timing
	time(&currentTime);
	//cout << "\nTime interval 6 " << difftime(currentTime,previousTime) << "\n";
	previousTime = currentTime; 

	// create polydata
	vtkSmartPointer<vtkPolyData> breastPoly =
		vtkSmartPointer<vtkPolyData>::New();
	breastPoly->SetPoints(breastPts);
	
	// add verticies to polydata
	breastPoly->SetVerts(breastVerts);
	
	// debug
	vtkSmartPointer<vtkXMLPolyDataWriter> writerd =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writerd->SetFileName(breastPointsFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
	writerd->SetInput(breastPoly);
#else
	writerd->SetInputData(breastPoly);
#endif
	writerd->Write();

	// surface mesh reconstruction
	vtkSmartPointer<vtkSurfaceReconstructionFilter> breastFilter =
		vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
		
#if VTK_MAJOR_VERSION <= 5
	breastFilter->SetInput(breastPoly);
#else
	breastFilter->SetInputData(breastPoly);
#endif
	breastFilter->Update();
	
	// fix bug in vtkSurfaceReconstructionFilter
	double pointSetBounds[6];
	breastPts->GetBounds(pointSetBounds);
	
	vtkSmartPointer<vtkContourFilter> breastContourFilter =
	vtkSmartPointer<vtkContourFilter>::New();
		
#if VTK_MAJOR_VERSION <= 5
	breastContourFilter->SetInput(breastFilter->GetOutput());
#else
	breastContourFilter->SetInputConnection(breastFilter->GetOutputPort());
#endif
	breastContourFilter->ComputeNormalsOn();
	breastContourFilter->SetValue(0, 0.0);
	breastContourFilter->Update();
	
	vtkSmartPointer<vtkReverseSense> breastSenseFilter =
		vtkSmartPointer<vtkReverseSense>::New();
	
#if VTK_MAJOR_VERSION <= 5
	breastSenseFilter->SetInput(breastContourFilter->GetOutput());
#else
	breastSenseFilter->SetInputConnection(breastContourFilter->GetOutputPort());
#endif
	breastSenseFilter->ReverseCellsOn();
	breastSenseFilter->ReverseNormalsOn();
	breastSenseFilter->Update();

	double surfaceBounds[6];
	breastSenseFilter->GetOutput()->GetBounds(surfaceBounds);
	
	double scaling[3];
	
	scaling[0] = (pointSetBounds[1]-pointSetBounds[0])/(surfaceBounds[1]-surfaceBounds[0]);
	scaling[1] = (pointSetBounds[3]-pointSetBounds[2])/(surfaceBounds[3]-surfaceBounds[2]);
	scaling[2] = (pointSetBounds[5]-pointSetBounds[4])/(surfaceBounds[5]-surfaceBounds[4]);

	vtkSmartPointer<vtkTransform> transp = 
		vtkSmartPointer<vtkTransform>::New();
	transp->Translate(pointSetBounds[0], pointSetBounds[2], pointSetBounds[4]);
	transp->Scale(scaling[0], scaling[1], scaling[2]);
	transp->Translate(-surfaceBounds[0], -surfaceBounds[2], -surfaceBounds[4]);

	vtkSmartPointer<vtkTransformPolyDataFilter> fixed = 
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	fixed->SetInputData(breastSenseFilter->GetOutput());
#if VTK_MAJOR_VERSION <= 5
	fixed->SetInput(breastSenseFilter->GetOutput());
#else
	fixed->SetInputData(breastSenseFilter->GetOutput());
#endif
	fixed->SetTransform(transp);
	fixed->Update();

	// clean mesh
	vtkSmartPointer<vtkCleanPolyData> cleanPoly =
		vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	cleanPoly->SetInput(fixed->GetOutput());
#else
	cleanPoly->SetInputConnection(fixed->GetOutputPort());
#endif
	
	// decimate mesh
	vtkSmartPointer<vtkDecimatePro> decimate = 
		vtkSmartPointer<vtkDecimatePro>::New();
	// debug
	decimate->SetTargetReduction(0.25);
#if VTK_MAJOR_VERSION <= 5
	decimate->SetInput(cleanPoly->GetOutput());
#else
	decimate->SetInputConnection(cleanPoly->GetOutputPort());
#endif
	decimate->Update();

	vtkSmartPointer<vtkPolyData> innerPoly =
		vtkSmartPointer<vtkPolyData>::New();
	innerPoly->ShallowCopy(decimate->GetOutput());

	// at this point have base surface
	// save base surface to disk
	vtkSmartPointer<vtkXMLPolyDataWriter> writer1 =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer1->SetFileName(baseShapeFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
	writer1->SetInput(decimate->GetOutput());
#else
	writer1->SetInputConnection(decimate->GetOutputPort());
#endif
	writer1->Write();

	cout << "done.\n";

	// timing
	time(&currentTime);
	cout << "Time interval 7: " << difftime(currentTime,previousTime) << " seconds.\n";
	previousTime = currentTime; 

	cout << "Creating voxelization...";
	// convert inner breast to voxel representation

	// create 3d imagedata
	vtkSmartPointer<vtkImageData> breast =
		vtkSmartPointer<vtkImageData>::New();

	double baseBound[6];
	innerPoly->GetBounds(baseBound);

	// extend bounds slightly to make room for nipple
	double baseFOV[3];
	for(int i=0; i<3; i++){
		baseFOV[i] = fabs(baseBound[2*i+1] - baseBound[2*i]);
	}

	double finalBound[6];
	for(int i=0; i<6; i++){
		finalBound[i] = baseBound[i];
	}
	finalBound[0] = finalBound[0] + 1.0;	// move backplane in by 1 mm to avoid surface edge
	finalBound[1] = finalBound[1] + nippleLen + 2.0;	// make room for nipple
	finalBound[2] = finalBound[2] - 0.01*baseFOV[1];
	finalBound[3] = finalBound[3] + 0.01*baseFOV[1];
	finalBound[4] = finalBound[4] - 0.01*baseFOV[2];
	finalBound[5] = finalBound[5] + 0.01*baseFOV[2];

	double spacing[3];
	for(int i=0; i<3; i++){
		spacing[i] = imgRes;
	}
	breast->SetSpacing(spacing);

	int dim[3];
	for(int i=0; i<3; i++){
		dim[i] = static_cast<int>(ceil((finalBound[2*i+1]-finalBound[2*i])/spacing[i]));
	}

	breast->SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);

	double origin[3];
	for(int i=0; i<3; i++){
		origin[i] = finalBound[2*i]+spacing[i]/2.0;
	}
	breast->SetOrigin(origin);

	// allocate unsigned char
#if VTK_MAJOR_VERSION <= 5
	breast->SetNumberOfScalarComponents(1);
	breast->SetScalarTypeToUnsignedChar();
	breast->AllocateScalars();
#else
	breast->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

	int originIndex[3] = {0, 0, 0};
	double originCoords[3];
	breast->GetPoint(breast->ComputePointId(originIndex),originCoords);

	// initialize to zero
	unsigned char* voxVal = static_cast<unsigned char *>(breast->GetScalarPointer());
	const int numElements = dim[0]*dim[1]*dim[2];
	for(int i=0; i<numElements; i++){
		*voxVal = tissue.bg;
		voxVal++;
	}

	// voxelize

	// Create the tree
	vtkSmartPointer<vtkCellLocator> innerLocator =
		vtkSmartPointer<vtkCellLocator>::New();
	innerLocator->SetDataSet(innerPoly);
	innerLocator->BuildLocator();


	// list of boundary voxels
	vtkSmartPointer<vtkIdList> boundaryList =
		vtkSmartPointer<vtkIdList>::New();

	// find intersect with top and bottom surface to voxelize breast
	// iterate over x and y values
	//#pragma omp parallel for
	// note: cell locator is not thread safe
	for(int i=0; i<dim[0]; i++){
		double xpos = origin[0]+i*spacing[0];
		int ijk[3];
		ijk[0] = i;
		for(int j=0; j<dim[1]; j++){

			double ypos = origin[1]+j*spacing[1];
			ijk[1] = j;
			// calculate z position of top surface
			double lineStart[3]; // end points of line
			double lineEnd[3];
			
			double tol = 0.005;
			double tval;
			vtkIdType intersectCell;
			int subId;

			double intersect[3]; // output position
			double pcoords[3];

			lineStart[0] = xpos;
			lineStart[1] = ypos;
			lineStart[2] = baseBound[4];

			lineEnd[0] = xpos;
			lineEnd[1] = ypos;
			lineEnd[2] = baseBound[5];

			if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
				tval, intersect, pcoords, subId, intersectCell)){

				// found intersection
				double topZ = intersect[2];

				// do other direction
				lineStart[2] = baseBound[5];
				lineEnd[2] = baseBound[4];

				if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					tval, intersect, pcoords, subId, intersectCell)){

					double bottomZ = intersect[2];

					// find nearest voxels to intersections
					int indexTop = static_cast<int>(floor((topZ-origin[2])/spacing[2]));
					int indexBottom = static_cast<int>(ceil((bottomZ-origin[2])/spacing[2]));

					// set edge voxels
					unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,indexTop));
					p[0] = boundVal;
					ijk[2] = indexTop;
					boundaryList->InsertUniqueId(breast->ComputePointId(ijk));
					
					p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,indexBottom));
					p[0] = boundVal;
					ijk[2] = indexBottom;
					boundaryList->InsertUniqueId(breast->ComputePointId(ijk));
					
					// set voxels between these 2 points to inner value;
					for(int k=indexTop+1; k<indexBottom; k++){
						unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
						p[0] = innerVal;
					}
				}
			}
		}
	}

	// find intersect with second set of directions
	// iterate over x and z values
	//#pragma omp parallel for
	for(int i=0; i<dim[0]; i++){
		double xpos = origin[0]+i*spacing[0];
		int ijk[3];
		ijk[0] = i;
		for(int j=0; j<dim[2]; j++){

			double zpos = origin[2]+j*spacing[2];
			ijk[2] = j;
			// calculate y position of top surface
			double lineStart[3]; // end points of line
			double lineEnd[3];

			double tol = 0.005;
			double tval; 
			vtkIdType intersectCell;
			int subId;  

			double intersect[3]; // output position
			double pcoords[3];

			lineStart[0] = xpos;
			lineStart[1] = baseBound[2];
			lineStart[2] = zpos;

			lineEnd[0] = xpos;
			lineEnd[1] = baseBound[3];
			lineEnd[2] = zpos;

			if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
				tval, intersect, pcoords, subId, intersectCell)){

				// found intersection
				double topY = intersect[1];

				// do other direction
				lineStart[1] = baseBound[3];
				lineEnd[1] = baseBound[2];

				if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					tval, intersect, pcoords, subId, intersectCell)){

					double bottomY = intersect[1];

					// find nearest voxels to intersections
					int indexTop = static_cast<int>(floor((topY-origin[1])/spacing[1]));
					int indexBottom = static_cast<int>(ceil((bottomY-origin[1])/spacing[1]));

					// set edge voxels
					unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,indexTop,j));
					p[0] = boundVal;
					ijk[1] = indexTop;
					boundaryList->InsertUniqueId(breast->ComputePointId(ijk));
					
					p = static_cast<unsigned char*>(breast->GetScalarPointer(i,indexBottom,j));
					p[0] = boundVal;
					ijk[1] = indexBottom;
					boundaryList->InsertUniqueId(breast->ComputePointId(ijk));
					
					// set voxels between these 2 points to inner value;
					for(int k=indexTop+1; k<indexBottom; k++){
						unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,k,j));
						p[0] = innerVal;
					}
				}
			}
		}
	}

	// find intersect with final set of directions
	// iterate over y and z values
	//#pragma omp parallel for
	for(int i=0; i<dim[1]; i++){
		double ypos = origin[1]+i*spacing[1];
		int ijk[3];
		ijk[1] = i;
		for(int j=0; j<dim[2]; j++){

			double zpos = origin[2]+j*spacing[2];
			ijk[2] = j;

			// calculate x position of top surface
			double lineStart[3]; // end points of line
			double lineEnd[3];

			double tol = 0.005;
			double tval; 
			vtkIdType intersectCell;
			int subId; 

			double intersect[3]; // output position
			double pcoords[3];

			lineStart[0] = baseBound[0];
			lineStart[1] = ypos;
			lineStart[2] = zpos;

			lineEnd[0] = baseBound[1];;
			lineEnd[1] = ypos;
			lineEnd[2] = zpos;

			if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
				tval, intersect, pcoords, subId, intersectCell)){

				// found intersection
				double topX = intersect[0];

				// do other direction
				lineStart[0] = baseBound[1];
				lineEnd[0] = baseBound[0];

				if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					tval, intersect, pcoords, subId, intersectCell)){

					double bottomX = intersect[0];

					// find nearest voxels to intersections
					int indexTop = static_cast<int>(floor((topX-origin[0])/spacing[0]));
					indexTop = (indexTop < 0) ? 0 : indexTop; 
					int indexBottom = static_cast<int>(ceil((bottomX-origin[0])/spacing[0]));
					indexBottom = (indexBottom > dim[0]-1) ? dim[0]-1 : indexBottom;
					indexBottom = (indexBottom < 0) ? 0 : indexBottom;

					// set edge voxels on front side only
					unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(indexBottom,i,j));
					p[0] = boundVal;
					ijk[0] = indexBottom;
					boundaryList->InsertUniqueId(breast->ComputePointId(ijk));

					// set voxels between these 2 points to inner value;
					for(int k=indexTop+1; k<indexBottom; k++){
						unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(k,i,j));
						p[0] = innerVal;
					}
				}
			}
		}
	}
	cout << "done.\n";

	/***********************
	Skin
	* **********************/
	// timing
	time(&currentTime);
	cout << "Time interval 8: " << difftime(currentTime,previousTime) << " seconds.\n";
	previousTime = currentTime; 

	// grow a skin layer
	cout << "Growing " << skinThick << " mm skin layer...";

	// create list of surrounding voxels to check
	vtkSmartPointer<vtkIntArray> checkVoxels =
		vtkSmartPointer<vtkIntArray>::New();

	checkVoxels->SetNumberOfComponents(3);

	

	int voxelThick = (int)ceil(skinThick/imgRes);

	for(int i=-1*voxelThick; i<=voxelThick; i++){
		for(int j=-1*voxelThick; j<=voxelThick; j++){
			for(int k=-1*voxelThick; k<=voxelThick; k++){
				if(sqrt((double)(i*i+j*j+k*k))*imgRes <= skinThick){
					// found a voxel to check
					checkVoxels->InsertNextTuple3(i,j,k);
				}
			}
		}
	}



	vtkIdType numCheck = checkVoxels->GetNumberOfTuples();

	int breastExtent[6];
	breast->GetExtent(breastExtent);

	// breast volume in voxels
	int breastVol = 0;

	// only add skin for x>0
	int minSkinXVox = static_cast<int>(ceil(-origin[0]/imgRes));
	int areolaVoxel = static_cast<int>(ceil(areolaRad/imgRes));

	// iterate over boundary voxels
	#pragma omp parallel for collapse(3) reduction(+:breastVol)
	for(int i=0; i<dim[0]; i++){
		for(int j=0; j<dim[1]; j++){
			for(int k=0; k<dim[2]; k++){
				// process now if not near nipple
				int ijk[3] = {i,j,k};
				double loc[3];
				breast->GetPoint(breast->ComputePointId(ijk),loc);
				double nipDist2 = vtkMath::Distance2BetweenPoints(loc,nipplePos);
				unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(p[0] == boundVal && nipDist2 > 4*areolaRad*areolaRad){
					if(i >= minSkinXVox){
						// we have a boundary voxel for skinning
						for(vtkIdType m=0; m<numCheck; m++){
							double offset[3];
							checkVoxels->GetTuple(m,offset);
							int a,b,c;
							a = i+(int)offset[0];
							b = j+(int)offset[1];
							c = k+(int)offset[2];
							if(a>= breastExtent[0] && a<=breastExtent[1] && b>=breastExtent[2] && b<=breastExtent[3] &&
								c>=breastExtent[4] && c<=breastExtent[5]){
								unsigned char* q =
									static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
								if(q[0] == tissue.bg){
									q[0] = tissue.skin;
								}
							}
						}
					}
					// can now change boundary to interior
					p[0] = innerVal;
				}
				if(p[0] == innerVal){
					breastVol += 1;
				}
			}
		}
	}

	// add skin thickness near nipple
	double skinThick2 = skinThick*2.0;
	
	// find point on mesh closest to center
	vtkSmartPointer<vtkPointLocator> locator =
		vtkSmartPointer<vtkPointLocator>::New();
	locator->SetDataSet(innerPoly);
	locator->BuildLocator();

	vtkIdType nipplePt;
	nipplePt = locator->FindClosestPoint(nipplePos);
	
	int nippleVoxel[3];	// coordinates of nipple base
	double nipplePCoords[3]; // parametric coordinates
	breast->ComputeStructuredCoordinates(nipplePos,nippleVoxel,nipplePCoords);
	
	for(int i=nippleVoxel[0]-2*areolaVoxel; i<=dim[0]-1; i++){
		for(int j=nippleVoxel[1]-2*areolaVoxel; j<=nippleVoxel[1]+2*areolaVoxel; j++){
			for(int k=nippleVoxel[2]-2*areolaVoxel; k<=nippleVoxel[2]+2*areolaVoxel; k++){
				int ijk[3] = {i,j,k};
				double loc[3];
				breast->GetPoint(breast->ComputePointId(ijk),loc);
				double nipDist2 = vtkMath::Distance2BetweenPoints(loc,nipplePos);
				unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(p[0] == boundVal && nipDist2 <= 4*areolaRad*areolaRad){
					double mySkinThick = skinThick + (skinThick2-skinThick)/(1+exp(12/areolaRad*(sqrt(nipDist2)-areolaRad)));
					int mySearchRad = static_cast<int>(ceil(mySkinThick/imgRes));
					for(int a=i-mySearchRad; a<=i+mySearchRad; a++){
						for(int b=j-mySearchRad; b<=j+mySearchRad; b++){
							for(int c=k-mySearchRad; c<=k+mySearchRad; c++){
								unsigned char* q = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
								if(q[0] == tissue.bg){
									// check distance
									double skinDist = imgRes*sqrt(static_cast<double>((a-i)*(a-i)+(b-j)*(b-j)+(c-k)*(c-k)));
									if(skinDist <= mySkinThick){
										q[0] = tissue.skin;
									}
								}
							}
						}
					}
					// can now change boundary to interior
					p[0] = innerVal;
					breastVol += 1;
				}
			}
		}
	}					
				

	cout << "done.\n";

	cout << "Breast volume: " << breastVol*pow(imgRes,3.0)/1000 << " cc ("<< breastVol << " voxels).\n";

	/***********************
	Nipple
	* **********************/

	// debug
	time(&currentTime);
	cout << "Time interval 9: " << difftime(currentTime,previousTime) << " seconds.\n";
	previousTime = currentTime; 

	// create nipple structure
	cout << "Creating nipple structure...";

	// squared radius
	double nippleRad2 = nippleRad*nippleRad;

	
	double nippleNorm[3];
	// extract normal
	vtkSmartPointer<vtkFloatArray> normals2 =
		vtkFloatArray::SafeDownCast(innerPoly->GetPointData()->GetNormals());


	normals2->GetTuple(nipplePt, nippleNorm);
	double nippleNormLen = sqrt(nippleNorm[0]*nippleNorm[0] +
		nippleNorm[1]*nippleNorm[1] + nippleNorm[2]*nippleNorm[2]);
	for(int i=0; i<3; i++){
		nippleNorm[i] = nippleNorm[i]/nippleNormLen;
	}
	// nippleNorm may be inward pointing
	if(nippleNorm[0] < 0.0){
		for(int i=0; i<3; i++){
			nippleNorm[i] = -1.0*nippleNorm[i];
		}
	}

	// create nipple

	// nipple function
	// superquadric (rad/nippleRad)^t+abs(len/nippleLen)^t <= 1 t = 2.5 - 8
	double nippleShape = 3.0;

	double searchRad = sqrt(2*nippleRad*nippleRad+nippleLen*nippleLen);

	

	int searchSpace[6];	// search region for nipple
	searchSpace[0] = nippleVoxel[0] - (int)ceil(searchRad/spacing[0]);
	if(searchSpace[0] < breastExtent[0]){
		searchSpace[0] = breastExtent[0];
	}
	searchSpace[1] = nippleVoxel[0] + (int)ceil(searchRad/spacing[0]);
	if(searchSpace[1] > breastExtent[1]){
		searchSpace[1] = breastExtent[1];
	}
	searchSpace[2] = nippleVoxel[1] - (int)ceil(searchRad/spacing[1]);
	if(searchSpace[2] < breastExtent[2]){
		searchSpace[2] = breastExtent[2];
	}
	searchSpace[3] = nippleVoxel[1] + (int)ceil(searchRad/spacing[1]);
	if(searchSpace[3] > breastExtent[3]){
		searchSpace[3] = breastExtent[3];
	}
	searchSpace[4] = nippleVoxel[2] - (int)ceil(searchRad/spacing[2]);
	if(searchSpace[4] < breastExtent[4]){
		searchSpace[4] = breastExtent[4];
	}
	searchSpace[5] = nippleVoxel[2] + (int)ceil(searchRad/spacing[2]);
	if(searchSpace[5] > breastExtent[5]){
		searchSpace[5] = breastExtent[5];
	}

	#pragma omp parallel for
	for(int i=searchSpace[0]; i<= searchSpace[1]; i++){
		for(int j=searchSpace[2]; j<= searchSpace[3]; j++){
			for(int k=searchSpace[4]; k<= searchSpace[5]; k++){
				// search cube, project onto normal and check if in nipple bound
				// only if outside breast interior
				unsigned char* q =
					static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(q[0] != innerVal){ // not in interior
					// get id for point
					vtkIdType id;
					int coord[3];
					coord[0] = i;
					coord[1] = j;
					coord[2] = k;
					id = breast->ComputePointId(coord);
					// get spatial coordinates of point
					double pos[3];
					breast->GetPoint(id,pos);
					// compute distance to nipple line and length along line
					double dist=0.0;
					double len=0.0;
					// project onto nipple line
					for(int m=0; m<3; m++){
						len += nippleNorm[m]*(pos[m]-nipplePos[m]);
					}
					// distance from nipple line
					for(int m=0; m<3; m++){
						dist += (pos[m]-nipplePos[m]-len*nippleNorm[m])*(pos[m]-nipplePos[m]-len*nippleNorm[m]);
					}
					dist = sqrt(dist);
					
					if(pow(dist/nippleRad,nippleShape)+pow(fabs(len)/nippleLen,nippleShape) <= 1.0){
						q[0] = tissue.nipple;
					}
				}
			}
		}
	}

	cout << "done.\n";


	// add chest muscle
	
	for(int j=0; j<dim[1]; j++){
	
		int muscleThick;
		
		if(leftSide){
			muscleThick = static_cast<int>(ceil((minSkinXVox-1)*(1-static_cast<double>(j*j)/(dim[1]*dim[1]))));
		}else{
			muscleThick = static_cast<int>(ceil((minSkinXVox-1)*(1-static_cast<double>((j-dim[1])*(j-dim[1]))/(dim[1]*dim[1]))));
		}
		
		for(int k=0; k<dim[2]; k++){
			for(int i=0; i<=muscleThick; i++){
				unsigned char* p =
					static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(p[0] == innerVal){
					p[0] = tissue.muscle;
				}
			}
		}
	}

	// save the volumetric data
	// save volume image data
	vtkSmartPointer<vtkXMLImageDataWriter> writer2 =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writer2->SetFileName(breastVoxelFilename.c_str());
	
#if VTK_MAJOR_VERSION <= 5
	writer2->SetInput(breast);
#else
	writer2->SetInputData(breast);
#endif
	writer2->Write();

	// for blank phantom make interior fat
	
	#pragma omp parallel for
	for(int i=0; i< dim[0]; i++){
		for(int j=0; j< dim[1]; j++){
			for(int k=0; k< dim[2]; k++){
				// search cube, project onto normal and check if in nipple bound
				// only if outside breast interior
				unsigned char* q =
					static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(q[0] == innerVal){ // in interior
					q[0] = tissue.fat;
				}
			}
		}
	}
	
	/************
	 * Save stuff
	 *
	 * */

	// debug
	time(&currentTime);
	cout << "Time interval 20: " << difftime(currentTime,previousTime) << " seconds.\n";
	previousTime = currentTime; 

	// save breast
	vtkSmartPointer<vtkXMLImageDataWriter> writerSeg5 =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writerSeg5->SetFileName(ductVoxelFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
	writerSeg5->SetInput(breast);
#else
	writerSeg5->SetInputData(breast);
#endif
	writerSeg5->Write();

	// save some projections in raw format

	unsigned int* skinProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* fatProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* nippleProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* muscleProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));

	unsigned int* skin2DProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* fat2DProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* nipple2DProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
	unsigned int* muscle2DProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));

	double sliceThick = 3.0;	// slice thickness (mm)

	// do near-center slice
	int sliceThickInd = (int)(ceil(sliceThick/imgRes));
	int nippleInd[3];
	double lcoords[3];
	breast->ComputeStructuredCoordinates(nipplePos, nippleInd, lcoords);

	int startSlice = nippleInd[1]-(int)(ceil(sliceThick/imgRes/2));

	for(int i=0; i<dim[0]; i++){
		for(int k=0; k<dim[2]; k++){
			for(int j=0; j<dim[1]; j++){
				unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				
				if(p[0] == tissue.fat){
					fatProj[i*dim[2]+k] = fatProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.skin){
					skinProj[i*dim[2]+k] = skinProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.nipple){
					nippleProj[i*dim[2]+k] = nippleProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.muscle){
					muscleProj[i*dim[2]+k] = muscleProj[i*dim[2]+k] + 1;
				}
			}
			for(int j=startSlice; j<startSlice+sliceThickInd; j++){
				unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
				if(p[0] == tissue.fat){
					fat2DProj[i*dim[2]+k] = fat2DProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.skin){
					skin2DProj[i*dim[2]+k] = skin2DProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.nipple){
					nipple2DProj[i*dim[2]+k] = nipple2DProj[i*dim[2]+k] + 1;
				} else if(p[0] == tissue.muscle){
					muscle2DProj[i*dim[2]+k] = muscle2DProj[i*dim[2]+k] + 1;
				}
			}
		}
	}

	cout << "Output projection dimensions = " << dim[0] << " by " << dim[2] << "\n";

	// save in raw format

	FILE *projFile;

	projFile = fopen("skinProj.dat","wb");
	fwrite(skinProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("fatProj.dat","wb");
	fwrite(fatProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("muscleProj.dat","wb");
	fwrite(muscleProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("nippleProj.dat","wb");
	fwrite(nippleProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);

	// output slice
	projFile = fopen("skin2DProj.dat","wb");
	fwrite(skin2DProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("fat2DProj.dat","wb");
	fwrite(fat2DProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("muscle2DProj.dat","wb");
	fwrite(muscle2DProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);
	projFile = fopen("nipple2DProj.dat","wb");
	fwrite(nipple2DProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
	fclose(projFile);

	free(skinProj);
	free(fatProj);
	free(muscleProj);
	free(nippleProj);

	free(skin2DProj);
	free(fat2DProj);
	free(muscle2DProj);
	free(nipple2DProj);

	// debug
	vtkSmartPointer<vtkTIFFWriter> tiffWriter =
		vtkSmartPointer<vtkTIFFWriter>::New();
	tiffWriter->SetCompressionToPackBits();
	tiffWriter->SetFileName("final.tiff");
#if VTK_MAJOR_VERSION <= 5
	tiffWriter->SetInput(breast);
#else
	tiffWriter->SetInputData(breast);
#endif
	tiffWriter->Write();
	

	time(&endTime);
	cout << "Phantom generation time: " << difftime(endTime,startTime) << " seconds\n";

	return EXIT_SUCCESS;
}

inline void statusBar(unsigned int current, unsigned int total, unsigned int width, unsigned int numUpdate){
	
	double prevNum = floor((double)(current-1)/(double)(total)*(double)(numUpdate));
	double curNum = floor((double)(current)/(double)(total)*(double)(numUpdate));

	if(curNum != prevNum){
		// draw new status bar
		cout << setw(3) << (int)((double)(current)/(double)(total)*100) << "% [";
		for(int x=0; x<=(int)curNum; x++){
			cout << "=";
		}
		for(int x=(int)curNum+1; x<=width; x++){
			cout << " ";
		}
		cout << "]\r" << std::flush;
	}
}


