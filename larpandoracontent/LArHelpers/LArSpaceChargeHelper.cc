/**
 *  @file   larpandoracontent/LArHelpers/LArSpaceChargeHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "LArSpaceChargeHelper.h"

using namespace pandora;

namespace lar_content
{
  //TGraph **g1_x = new TGraph*[7];
  //TGraph **g2_x = new TGraph*[7];
  //TGraph **g3_x = new TGraph*[7];
  //TGraph **g4_x = new TGraph*[7];
  // TGraph **g5_x = new TGraph*[7];

  //TGraph **g1_y = new TGraph*[7];
  // TGraph **g2_y = new TGraph*[7];
  // TGraph **g3_y = new TGraph*[7];
  // TGraph **g4_y = new TGraph*[7];
  // TGraph **g5_y = new TGraph*[7];
  // TGraph **g6_y = new TGraph*[7];

  // TGraph **g1_z = new TGraph*[7];
  // TGraph **g2_z = new TGraph*[7];
  // TGraph **g3_z = new TGraph*[7];
  // TGraph **g4_z = new TGraph*[7];

  // TF1 *f1_x = new TF1("f1_x","pol6");
  // TF1 *f2_x = new TF1("f2_x","pol6");
  // TF1 *f3_x = new TF1("f3_x","pol6");
  // TF1 *f4_x = new TF1("f4_x","pol6");
  // TF1 *f5_x = new TF1("f5_x","pol6");
  // TF1 *fFinal_x = new TF1("fFinal_x","pol4");

  // TF1 *f1_y = new TF1("f1_y","pol5");
  // TF1 *f2_y = new TF1("f2_y","pol5");
  // TF1 *f3_y = new TF1("f3_y","pol5");
  //TF1 *f4_y = new TF1("f4_y","pol5");
  // TF1 *f5_y = new TF1("f5_y","pol5");
  // TF1 *f6_y = new TF1("f6_y","pol5");
  // TF1 *fFinal_y = new TF1("fFinal_y","pol5");

  // TF1 *f1_z = new TF1("f1_z","pol4");
  // TF1 *f2_z = new TF1("f2_z","pol4");
  // TF1 *f3_z = new TF1("f3_z","pol4");
  // TF1 *f4_z = new TF1("f4_z","pol4");
  //  TF1 *fFinal_z = new TF1("fFinal_z","pol3");

//------------------------------------------------------------------------------------------------------------------------------------------

void LArSpaceChargeHelper::Configure(std::string fileName)
{
  //TFile infile(fileName.c_str(), "READ");
  std::cout << fileName << std::endl;

  //if(!infile.IsOpen())
  //{
  //    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
  //    std::cout << "Could not find the space charge effect file '" << fileName << "'!\n";
  // }

    for(int i = 0; i < 5; i++)
    {
      std::cout << i << "spacecharge" << std::endl;
      //g1_x[i] = (TGraph*)infile.Get(Form("deltaX/g1_%d",i));
      // g2_x[i] = (TGraph*)infile.Get(Form("deltaX/g2_%d",i));
      //g3_x[i] = (TGraph*)infile.Get(Form("deltaX/g3_%d",i));   
      //g4_x[i] = (TGraph*)infile.Get(Form("deltaX/g4_%d",i));
      //g5_x[i] = (TGraph*)infile.Get(Form("deltaX/g5_%d",i));

      // g1_y[i] = (TGraph*)infile.Get(Form("deltaY/g1_%d",i));
      // g2_y[i] = (TGraph*)infile.Get(Form("deltaY/g2_%d",i));
      // g3_y[i] = (TGraph*)infile.Get(Form("deltaY/g3_%d",i));   
      // g4_y[i] = (TGraph*)infile.Get(Form("deltaY/g4_%d",i));
      // g5_y[i] = (TGraph*)infile.Get(Form("deltaY/g5_%d",i));
      // g6_y[i] = (TGraph*)infile.Get(Form("deltaY/g6_%d",i));

      // g1_z[i] = (TGraph*)infile.Get(Form("deltaZ/g1_%d",i));
      // g2_z[i] = (TGraph*)infile.Get(Form("deltaZ/g2_%d",i));
      // g3_z[i] = (TGraph*)infile.Get(Form("deltaZ/g3_%d",i));   
      // g4_z[i] = (TGraph*)infile.Get(Form("deltaZ/g4_%d",i));
    }

    //g1_x[5] = (TGraph*)infile.Get("deltaX/g1_5");
    //g2_x[5] = (TGraph*)infile.Get("deltaX/g2_5");
    //g3_x[5] = (TGraph*)infile.Get("deltaX/g3_5");   
    //g4_x[5] = (TGraph*)infile.Get("deltaX/g4_5");
    //g5_x[5] = (TGraph*)infile.Get("deltaX/g5_5");

    //g1_y[5] = (TGraph*)infile.Get("deltaY/g1_5");
    //g2_y[5] = (TGraph*)infile.Get("deltaY/g2_5");
    //g3_y[5] = (TGraph*)infile.Get("deltaY/g3_5");   
    //g4_y[5] = (TGraph*)infile.Get("deltaY/g4_5");
    // g5_y[5] = (TGraph*)infile.Get("deltaY/g5_5");
    //g6_y[5] = (TGraph*)infile.Get("deltaY/g6_5");

    // g1_x[6] = (TGraph*)infile.Get("deltaX/g1_6");
    //g2_x[6] = (TGraph*)infile.Get("deltaX/g2_6");
    //g3_x[6] = (TGraph*)infile.Get("deltaX/g3_6");   
    //g4_x[6] = (TGraph*)infile.Get("deltaX/g4_6");
    //g5_x[6] = (TGraph*)infile.Get("deltaX/g5_6");

    //infile.Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(CartesianVector &positionVector)
{
    return positionVector - GetPositionOffset(positionVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::GetPositionOffset(CartesianVector &positionVector)
{
    CartesianVector transformedPosition(TransformPosition(positionVector));
    CartesianVector positionOffsets(-GetParametricPositionOffset(transformedPosition, "X"), GetParametricPositionOffset(transformedPosition, "Y"), GetParametricPositionOffset(transformedPosition, "Z")); 

    

    return positionOffsets;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::GetParametricPositionOffset(CartesianVector &transformedPosition, std::string axisLabel)
{
  //double parA[6][7], parB[6];
  std::cout<< transformedPosition.GetZ()<<"transformed z" <<std::endl;
     
    for (int j = 0; j < 6; j++)
    {
        for(int i = 0; i < 7; i++)
	  //parA[j][i] = 0.0;
	  std::cout << i << "primary" << std::endl;
	  //parB[j] = 0.0;
    }
     
    if (axisLabel == "X")
    {
        for(int j = 0; j < 7; j++)
        {
	  //parA[0][j] = g1_x[j]->Eval(transformedPosition.GetZ());
	  //parA[1][j] = g2_x[j]->Eval(transformedPosition.GetZ());
	  //parA[2][j] = g3_x[j]->Eval(transformedPosition.GetZ());
	  //parA[3][j] = g4_x[j]->Eval(transformedPosition.GetZ());
	  //parA[4][j] = g5_x[j]->Eval(transformedPosition.GetZ());
        }
     
	// f1_x->SetParameters(parA[0]);
        //f2_x->SetParameters(parA[1]);
        //f3_x->SetParameters(parA[2]);
        //f4_x->SetParameters(parA[3]);
        //f5_x->SetParameters(parA[4]);
    }
    else if (axisLabel == "Y")
    {
        for(int j = 0; j < 6; j++)
        {
	  //parA[0][j] = g1_y[j]->Eval(transformedPosition.GetZ());
	  //parA[1][j] = g2_y[j]->Eval(transformedPosition.GetZ());
	  //parA[2][j] = g3_y[j]->Eval(transformedPosition.GetZ());
	  //parA[3][j] = g4_y[j]->Eval(transformedPosition.GetZ());
	  //parA[4][j] = g5_y[j]->Eval(transformedPosition.GetZ());
	  //parA[5][j] = g6_y[j]->Eval(transformedPosition.GetZ());
        }
    
        //f1_y->SetParameters(parA[0]);
        //f2_y->SetParameters(parA[1]);
        //f3_y->SetParameters(parA[2]);
        //f4_y->SetParameters(parA[3]);
        //f5_y->SetParameters(parA[4]);
        //f6_y->SetParameters(parA[5]);
    }
    else if (axisLabel == "Z")
    {
        for (int j = 0; j < 5; j++)
        {
	  //parA[0][j] = g1_z[j]->Eval(transformedPosition.GetZ());
	  //parA[1][j] = g2_z[j]->Eval(transformedPosition.GetZ());
	  ///parA[2][j] = g3_z[j]->Eval(transformedPosition.GetZ());
	  //parA[3][j] = g4_z[j]->Eval(transformedPosition.GetZ());
        }
     
        //f1_z->SetParameters(parA[0]);
        //f2_z->SetParameters(parA[1]);
        //f3_z->SetParameters(parA[2]);
        //f4_z->SetParameters(parA[3]);
    }
     
    //float aValNew, bValNew;
     
    if (axisLabel == "Y")
    {
      //aValNew = transformedPosition.GetX();
      //bValNew = transformedPosition.GetY();
    }
    else
    {
      //aValNew = transformedPosition.GetY();
      //bValNew = transformedPosition.GetX();
    }
     
    float oneDimensionalOffset(0.f);

    if (axisLabel == "X")
    {
      //parB[0] = f1_x->Eval(aValNew);
      //parB[1] = f2_x->Eval(aValNew);
        //parB[2] = f3_x->Eval(aValNew);
        //parB[3] = f4_x->Eval(aValNew);
        //parB[4] = f5_x->Eval(aValNew);
     
        //fFinal_x->SetParameters(parB);
        //oneDimensionalOffset = 100.0*fFinal_x->Eval(bValNew);
    }
    else if (axisLabel == "Y")
    {
      //parB[0] = f1_y->Eval(aValNew);
      //parB[1] = f2_y->Eval(aValNew);
      //parB[2] = f3_y->Eval(aValNew);
      //parB[3] = f4_y->Eval(aValNew);
      //parB[4] = f5_y->Eval(aValNew);
      //parB[5] = f6_y->Eval(aValNew);
     
        //fFinal_y->SetParameters(parB);
        //oneDimensionalOffset = 100.0*fFinal_y->Eval(bValNew);
    }

    else if(axisLabel == "Z")
    {
      //parB[0] = f1_z->Eval(aValNew);
      //parB[1] = f2_z->Eval(aValNew);
      //parB[2] = f3_z->Eval(aValNew);
      //parB[3] = f4_z->Eval(aValNew);
        
        //fFinal_z->SetParameters(parB);
        //oneDimensionalOffset = 100.0*fFinal_z->Eval(bValNew);
    }
  
    return oneDimensionalOffset;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArSpaceChargeHelper::TransformPosition(CartesianVector &positionVector)
{
    CartesianVector transformedPosition(TransformX(positionVector.GetX()), TransformY(positionVector.GetY()), TransformZ(positionVector.GetZ()));
    return transformedPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformX(float xPosition)
{
    return (1.25 - (2.50/2.56)*(xPosition/100.0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformY(float yPosition)
{
    return ((2.50/2.33)*((yPosition/100.0)+1.165) - 1.25);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArSpaceChargeHelper::TransformZ(float zPosition)
{
    return (10.0/10.37)*(zPosition/100.0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
