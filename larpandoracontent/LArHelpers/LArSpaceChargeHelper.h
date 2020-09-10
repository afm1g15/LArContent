#ifndef LAR_SPACE_CHARGE_HELPER_H
#define LAR_SPACE_CHARGE_HELPER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

//#include "TGraph.h"
//#include "TFile.h"
//#include "TF1.h"

#include "Objects/CartesianVector.h"
 
namespace lar_content
{

    class LArSpaceChargeHelper
    {
    public:

        static void Configure(std::string fileName);
        static pandora::CartesianVector GetSpaceChargeCorrectedPosition(pandora::CartesianVector &inputPoistion);
        static pandora::CartesianVector GetPositionOffset(pandora::CartesianVector &positionOffset);

    protected:

        static float GetParametricPositionOffset(pandora::CartesianVector &transformedPosition, std::string axisLabel);
        static pandora::CartesianVector TransformPosition(pandora::CartesianVector &inputPosition);
        static float TransformX(float xPosition);
        static float TransformY(float yPosition);
        static float TransformZ(float zPosition);

    };
} 
 
#endif // LAR_SPACE_CHARGE_HELPER_H
