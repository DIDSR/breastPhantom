Configuration Files Parameters
==============================

The configuration file contains all the parameters necessary to generate a breast phantom.  They are divided into several categories, primarily by tissue type.
In the following tables each parameter is given a brief description.

base parameters
----------------
================== ========== =========================================================
Name               Type       Notes
================== ========== =========================================================
base.outputDir     string     directory output files are saved to
base.imgRes        float (mm) voxel size
base.skinThick     float (mm) thickness of skin
base.nippleLen     float (mm) length of nipple
base.nippleRad     float (mm) radius of nipple
base.areolaRad     float (mm) radius of thicker skin around nipple
base.leftBreast    Boolean    true for left breast, false for right breast
base.targetFatFrac float (mm) desired fraction of interior breast volume containing fat
base.seed          integer    random number seed (chosen randomly if not specified)
================== ========== =========================================================

shape parameters
----------------

===================== =============== ================================================
Name                  Type            Notes
===================== =============== ================================================
shape.ures            float           mesh resolution for breast shape
shape.vres            float           mesh resolution for breast shape
shape.pointSep        float (mm)      minimum point separation for point cloud
shape.ringWidth       float (mm)      thickness of muscle backing layer
shape.ringSep         float (mm)      mesh resolution for muscle layer
shape.featureAngle    float (degrees) minimum angle to preserve during mesh smoothing
shape.targetReduction float           fraction of triangles to remain after decimation
shape.a1b             float           scale of breast bottom
shape.a1t             float           scale of breast top
shape.a2l             float           scale of breast left side
shape.a2r             float           scale of breast right side
shape.a3              float           breast outward scale
shape.eps1            float           u direction quadric shape exponent
shape.eps2            float           v direction quadric shape exponent
shape.doPtosis        Boolean         if true include ptosis in shape
shape.ptosisB0        float           ptosis parameter B0
shape.ptosisB1        float           ptosis parameter B1
shape.doTurn          Boolean         if true include turn deformation
shape.turnC0          float           turn deformation parameter C0
shape.turnC1          float           turn deformation parameter C1
shape.doTopShape      Boolean         if true include top shape deformation
shape.topShapeS0      float           top shape deformation parameter S0
shape.topShapeS1      float           top shape deformation parameter S1
shape.topShapeT0      float           top shape deformation parameter T0
shape.topShapeT1      float           top shape deformation parameter T1
shape.doFlattenSide   Boolean         if true include flatten side deformation
shape.flattenSideG0   float           flatten side deformation parameter G0
shape.flattenSideG1   float           flatten side deformation parameter G1
shape.doTurnTop       Boolean         if true include turn top deformation
shape.turnTopH0       float           turn top deformation parameter H0
shape.turnTopH1       float           turn top deformation parameter H1
===================== =============== ================================================

glandular compartment parameters
--------------------------------

================================== ========== ==================================================================================
Name                               Type       Notes
================================== ========== ==================================================================================
compartments.num                   integer    number of glandular compartments
compartments.seedBaseDist          float (mm) distance along nipple line of compartment seed base
compartments.backFatBufferFrac     float      fraction of breast at muscle layer forced to be fat
compartments.numBackSeeds          integer    number of back plane Voronoi seeds
compartments.angularJitter         float      compartment Voronoi seed jitter (fraction of subtended angle)
compartments.zJitter               float (mm) maximum compartment Voronoi seed jitter in direction of nipple
compartments.maxFracRadialDist     float      maximum radial distance from base seed as a fraction of distance to breast surface
compartments.minFracRadialDist     float      minimum radial distance from base seed as a fraction of distance to breast surface
compartments.minScaleNippleDir     float      minimum scale in nipple direction
compartments.maxScaleNippleDir     float      maximum scale in nipple direction
compartments.minScale              float      minimum scale in non-nipple direction
compartments.maxScale              float      maximum scale in non-nipple direction
compartments.minGlandStrength      float      minimum gland strength
compartments.maxGlandStrength      float      maximum gland strength
compartments.maxDeflect            float      maximum compartment deflection angle from pointing towards nipple (fraction of pi)
compartments.minSkinScaleNippleDir float      minimum scale skin seeds in nipple direction
compartments.maxSkinScaleNippleDir float      maximum scale skin seeds in nipple direction
compartments.minSkinScale          float      minimum scale skin in non-nipple direction
compartments.maxSkinScale          float      maximum scale skin in non-nipple direction
compartments.skinStrength          float      skin strength
compartments.backScale             float      back scale
compartments.backStrength          float      back strength
compartments.nippleScale           float      nipple scale
compartments.nippleStrength        float      nipple strength
compartments.voronSeedRadius       float (mm) radius from point to check Voronoi seed distance
================================== ========== ==================================================================================

TDLU parameters
---------------

============== ========== ===================
Name           Type       Notes
============== ========== ===================
TDLU.maxLength float (mm) maximum TDLU length
TDLU.minLength float (mm) minimum TDLU length
TDLU.maxWidth  float (mm) maximum TDLU width
TDLU.minWidth  float (mm) minimum TDLU width
============== ========== ===================


Perlin noise parameters
-----------------------

==================== ======= ====================================
Name                 Type    Notes
==================== ======= ====================================
perlin.maxDeviation  float   maximum fraction of radius deviation
perlin.frequency     float   starting frequency
perlin.lacunarity    float   octave frequency multiplier
perlin.persistence   float   octave signal decay
perlin.numOctaves    integer number of frequency octaves
perlin.xNoiseGen     integer x direction noise generation seed
perlin.yNoiseGen     integer y direction noise generation seed
perlin.zNoiseGen     integer z direction noise generation seed
perlin.seedNoiseGen  integer seed noise generation
perlin.shiftNoiseGen integer shift noise generation seed
==================== ======= ====================================


Compartment boundary noise parameters
-------------------------------------

===================== ===== ======================================
Name                  Type  Notes
===================== ===== ======================================
boundary.maxDeviation float maximum fraction of distance deviation
boundary.frequency    float starting frequency
boundary.lacunarity   float octave frequency multiplier
boundary.persistence  float octave signal decay
===================== ===== ======================================


fat lobule boundary perturbation noise parameters
-------------------------------------------------

==================== ===== ======================================
Name                 Type  Notes
==================== ===== ======================================
perturb.maxDeviation float maximum fraction of distance deviation
perturb.frequency    float starting frequency
perturb.lacunarity   float octave frequency multiplier
perturb.persistence  float octave signal decay
==================== ===== ======================================


fat glandular buffer noise parameters
-------------------------------------

=================== ===== ======================================
Name                Type  Notes
=================== ===== ======================================
buffer.maxDeviation float maximum fraction of distance deviation
buffer.frequency    float starting frequency
buffer.lacunarity   float octave frequency multiplier
buffer.persistence  float octave signal decay
=================== ===== ======================================


Voronoi segmentation variables
------------------------------

=============================== ============= ============================================
Name                            Type          Notes
=============================== ============= ============================================
voronoi.fatInFatSeedDensity     float (mm^-3) fat voronoi seed density
voronoi.fatInGlandSeedDensity   float (mm^-3) fat voronoi seed in glandular tissue density
voronoi.glandInGlandSeedDensity float (mm^-3) glandular voronoi seed density
voronoi.TDLUDeflectMax          float         maximum deflection (fraction of pi)
voronoi.minScaleLenTDLU         float         minimum length scale
voronoi.maxScaleLenTDLU         float         maximum length scale
voronoi.minScaleWidTDLU         float         minimum width scale
voronoi.maxScaleWidTDLU         float         maximum width scale
voronoi.minStrTDLU              float         minimum strength
voronoi.maxStrTDLU              float         maximum strength
voronoi.fatInFatDeflectMax      float         maximum deflection (fraction of pi)
voronoi.minScaleLenFatInFat     float         minimum length scale
voronoi.maxScaleLenFatInFat     float         maximum length scale
voronoi.minScaleWidFatInFat     float         minimum width scale
voronoi.maxScaleWidFatInFat     float         maximum width scale
voronoi.minStrFatInFat          float         minimum strength
voronoi.maxStrFatInFat          float         maximum strength
voronoi.fatInGlandDeflectMax    float         maximum deflection (fraction of pi)
voronoi.minScaleLenFatInGland   float         minimum length scale
voronoi.maxScaleLenFatInGland   float         maximum length scale
voronoi.minScaleWidFatInGland   float         minimum width scale
voronoi.maxScaleWidFatInGland   float         maximum width scale
voronoi.minStrFatInGland        float         minimum strength
voronoi.maxStrFatInGland        float         maximum strength
voronoi.glandInGlandDeflectMax  float         maximum deflection (fraction of pi)
voronoi.minScaleLenGlandInGland float         minimum length scale
voronoi.maxScaleLenGlandInGland float         maximum length scale
voronoi.minScaleWidGlandInGland float         minimum width scale
voronoi.maxScaleWidGlandInGland float         maximum width scale
voronoi.minStrGlandInGland      float         minimum strength
voronoi.maxStrGlandInGland      float         maximum strength
voronoi.seedRadius              float (mm)    check seeds in radius
=============================== ============= ============================================

fat lobule parameters
---------------------

================= ========== ===========================================================================
Name              Type       Notes
================= ========== ===========================================================================
fat.minLobuleAxis float (mm) min lobule axis length
fat.maxLobuleAxis float (mm) max lobule axis length
fat.minAxialRatio float      axial ratio min
fat.maxAxialRatio float      axial ratio max
fat.minLobuleGap  float      minimum ligament separation between lobules
fat.maxCoeffStr   float      maximum of absolute value of Fourier coefficient as fraction of main radius
fat.minCoeffStr   float      minimum of absolute value of Fourier coefficient as fraction of main radius
fat.maxLobuleTry  integer    maximum number of trial lobules
================= ========== ===========================================================================

ductal tree parameters
----------------------

================== ========== ==========================================
Name               Type       Notes
================== ========== ==========================================
ductTree.maxBranch integer    target number of branches
ductTree.maxGen    integer    maximum generation
ductTree.initRad   float (mm) initial radius of tree
ductTree.nFillX    integer    number of voxels for tree density tracking
ductTree.nFillY    integer    number of voxels for tree density tracking
ductTree.nFillZ    integer    number of voxels for tree density tracking
================== ========== ==========================================

ductal branch parameters
------------------------

.. tabularcolumns:: |l|c|p{8cm}|

==================== ========== ==================================================================================================================
Name                 Type       Notes
==================== ========== ==================================================================================================================
ductBr.minLen0       float (mm) minimum and maximum branch lengths per level
ductBr.maxLen0       float (mm) minimum and maximum branch lengths per level
ductBr.minLen1       float (mm) minimum and maximum branch lengths per level
ductBr.maxLen1       float (mm) minimum and maximum branch lengths per level
ductBr.minLen2       float (mm) minimum and maximum branch lengths per level
ductBr.maxLen2       float (mm) minimum and maximum branch lengths per level
ductBr.minLenDefault float (mm) minimum and maximum branch lengths per level
ductBr.maxLenDefault float (mm) minimum and maximum branch lengths per level
ductBr.maxChild      integer    maximum number of children
ductBr.childMinRad   float (mm) minimum branch radius to have children
ductBr.childLevBound integer
ductBr.child00       float
ductBr.child01       float
ductBr.child02       float
ductBr.child10       float
ductBr.child11       float
ductBr.child12       float
ductBr.child20       float
ductBr.child21       float
ductBr.child22       float
ductBr.child30       float
ductBr.child31       float
ductBr.child32       float
ductBr.child40       float
ductBr.child41       float
ductBr.child42       float
ductBr.minRadFrac    float      minimum starting radius as a fraction of parent end radius
ductBr.maxRadFrac    float      maximum starting radius as a fraction of parent end radius
ductBr.radFrac0      float      starting radius as fraction of parent end radius for first child
ductBr.minAngle      float      min angle between parent end direction and child start direction for children after first (fraction of pi radians)
ductBr.maxAngle      float      max angle between parent end direction and child start direction for children after first (fraction of pi radians)
==================== ========== ==================================================================================================================

ductal segment parameters
-------------------------

========================= ========== ===========================================================================
Name                      Type       Notes
========================= ========== ===========================================================================
ductSeg.lengthBetaA       float      length distribution shape parameters
ductSeg.lengthBetaB       float      length distribution shape parameters
ductSeg.radiusBetaA       float      radius distribution shape parameters
ductSeg.radiusBetaB       float      radius distribution shape parameters
ductSeg.minLen            float (mm) min and max segment length
ductSeg.maxLen            float (mm) min and max segment length
ductSeg.maxCurvRad        float (mm) maximum radius of curvature
ductSeg.maxCurvFrac       float      maximum length of segment based on  curvature (fraction of pi radians)
ductSeg.ductSeg.minEndRad float      min and max end radius as fraction of start radius
ductSeg.maxEndRad         float      min and max end radius as fraction of start radius
ductSeg.angleWt           float      cost function preferential angle weighting
ductSeg.densityWt         float      cost function density weighting
ductSeg.numTry            integer    number of trial segments to generate
ductSeg.maxTry            integer    maximum number of segments to generate before giving up and reducing length
ductSeg.absMaxTry         integer    total number of segment tries before completely giving up
ductSeg.roiStep           float (mm) step size for checking segment is valid
========================= ========== ===========================================================================

vessel tree parameters
----------------------
==================== ========== ==========================================
Name                 Type       Notes
==================== ========== ==========================================
vesselTree.maxBranch integer    target number of branches
vesselTree.maxGen    integer    maximum generation
vesselTree.initRad   float (mm) initial radius of tree
vesselTree.nFillX    integer    number of voxels for tree density tracking
vesselTree.nFillY    integer    number of voxels for tree density tracking
vesselTree.nFillZ    integer    number of voxels for tree density tracking
==================== ========== ==========================================

vessel branch parameters
------------------------

.. tabularcolumns:: |l|c|p{8cm}|

====================== ========== ==================================================================================================================
Name                   Type       Notes
====================== ========== ==================================================================================================================
vesselBr.minLen0       float (mm) minimum and maximum branch lengths per level
vesselBr.maxLen0       float (mm) minimum and maximum branch lengths per level
vesselBr.minLen1       float (mm) minimum and maximum branch lengths per level
vesselBr.maxLen1       float (mm) minimum and maximum branch lengths per level
vesselBr.minLen2       float (mm) minimum and maximum branch lengths per level
vesselBr.maxLen2       float (mm) minimum and maximum branch lengths per level
vesselBr.minLenDefault float (mm) minimum and maximum branch lengths per level
vesselBr.maxLenDefault float (mm) minimum and maximum branch lengths per level
vesselBr.maxChild      int        maximum number of children
vesselBr.childMinRad   float      minimum branch radius to have children (mm)
vesselBr.childLevBound integer
vesselBr.child00       float
vesselBr.child01       float
vesselBr.child02       float
vesselBr.child10       float
vesselBr.child11       float
vesselBr.child12       float
vesselBr.child20       float
vesselBr.child21       float
vesselBr.child22       float
vesselBr.child30       float
vesselBr.child31       float
vesselBr.child32       float
vesselBr.child40       float
vesselBr.child41       float
vesselBr.child42       float
vesselBr.minRadFrac    float      minimum starting radius as a fraction of parent end radius
vesselBr.maxRadFrac    float      maximum starting radius as a fraction of parent end radius
vesselBr.radFrac0      float      starting radius as fraction of parent end radius for first child
vesselBr.minAngle      float      min angle between parent end direction and child start direction
                                  for children after first (fraction of pi radians)
vesselBr.maxAngle      float      max angle between parent end direction and child start direction
                                  for children after first (fraction of pi radians)
====================== ========== ==================================================================================================================


vessel segment parameters
-------------------------

===================== ========== ===========================================================================
Name                  Type       Notes
===================== ========== ===========================================================================
vesselSeg.lengthBetaA float      length distribution shape parameters
vesselSeg.lengthBetaB float      length distribution shape parameters
vesselSeg.radiusBetaA float      radius distribution shape parameters
vesselSeg.radiusBetaB float      radius distribution shape parameters
vesselSeg.minLen      float (mm) min and max segment length
vesselSeg.maxLen      float (mm) min and max segment length
vesselSeg.maxCurvRad  float (mm) maximum radius of curvature
vesselSeg.maxCurvFrac float      maximum length of segment based on  curvature (fraction of pi radians)
vesselSeg.minEndRad   float      min and max end radius as fraction of start radius
vesselSeg.maxEndRad   float      min and max end radius as fraction of start radius
vesselSeg.angleWt     float      cost function preferential angle weighting
vesselSeg.densityWt   float      cost function density weighting
vesselSeg.numTry      integer    number of trial segments to generate
vesselSeg.maxTry      integer    maximum number of segments to generate before giving up and reducing length
vesselSeg.absMaxTry   integer    total number of segment tries before completely giving up
vesselSeg.roiStep     float (mm) step size for checking segment is valid
===================== ========== ===========================================================================
