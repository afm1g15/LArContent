//Get reconstructed pfos and pull out neutrino type and energy

#include "larpandoracontent/LArControlFlow/NDQualitiesAlg.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

  NDQualitiesAlg::NDQualitiesAlg() {
  }

  NDQualitiesAlg::~NDQualitiesAlg() {
    try
      {
	//	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree", "output.root", "UPDATE"));
	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreeb", "outputb.root", "UPDATE"));
	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreec", "outputc.root", "UPDATE"));
	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreed", "outputd.root", "UPDATE"));
	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreee", "outpute.root", "UPDATE"));
	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreef", "outputf.root", "UPDATE"));
      }
    catch (const StatusCodeException &)
      {
	std::cout << " Unable to write tree  to file " << std::endl;
      }
    
  }
  //-----------------------------------------------------------------------

  StatusCode NDQualitiesAlg::Run() {
    std::cout << "///////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;
    std::cout << "/////////////////////START OF ALG RUN\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;
    std::cout << "///////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    std::cout << "parentMinX " << parentMinX << std::endl;
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    std::cout << "parentMaxX " << parentMaxX << std::endl;
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());


    const PfoList *pPfoList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));  //RecreatedPfos
    HitChargeVector forwardsFitPoints, backwardsFitPoints;
    float forwardsChiSquared = 0;
    float backwardsChiSquared = 0;
    float minchisquaredval = 0;
    int numberHitsToConsider = 0;
    double neutrinoenergy = 0;
    double andyenergytotal = 0;
    double andyenergytotal_f_L = 0;
    double andyenergytotal_b_L = 0;
    double energybylengthtotal = 0;
    double energybylengthtotalb = 0;

    // double energybylengthtotal_end = 0;
    // double energybylengthtotalb_end = 0;
    double fdeltaEfromminchi = 0;
    double bdeltaEfromminchi = 0;
    double front_b_dE = 0;
    double front_f_dE = 0;
    double prednuenergy = 0.0;
    for (const ParticleFlowObject *const pPfo : *pPfoList)
      {
	std::cout << "            " << std::endl;
	std::cout << "/////////////////////NOW LOOKING AT A NEW PFO\\\\\\\\\\\\\\\\\\\\\\\\" << std::endl;

	std::cout << "NDQualitiesAlg :: pPfo " << pPfo << std::endl;
	std::cout << "NDQualitiesAlg ::  PDG code: " << pPfo->GetParticleId()  << " Daughters:  " << pPfo->GetNDaughterPfos() << " Parents:  " << pPfo->GetNParentPfos() << std::endl;
	std::cout << "NDQualitiesAlg ::  mass: " << pPfo->GetMass()*1000 << " MeV" << std::endl; //given in GeV

	double energybylengthtotal_end = 0;
	double energybylengthtotalb_end = 0;
	double fdeltaEfromminchi = 0;
	double bdeltaEfromminchi = 0;
	for (const ParticleFlowObject *const pPfo : *pPfoList)
	  {

	std::cout << "NDQualitiesAlg: pPfo " << pPfo << std::endl;
	std::cout << "NDQualitiesAlg: PDG code: " << pPfo->GetParticleId()  << " Daughters:  " << pPfo->GetNDaughterPfos() << " Parents:  " << pPfo->GetNParentPfos() << std::endl;
	std::cout << "NDQualitiesAlg: mass: " << pPfo->GetMass()*1000 << " MeV" << std::endl; //given in GeV



	double energypfo = 0;
	double energycount = 0;
	double predenergy = 0.0;
	CaloHitList caloHitList;
	ClusterList clusterList = pPfo->GetClusterList();
	const Cluster* pCluster = clusterList.back();

	//-------Find the type of track/shower presented------------------------------------------------------------------
	int w = 0;
	for (auto c :clusterList) {
	  if (w == 3) {
	    pCluster = c;
	  }
	  w++;
	}
	
	bool isContained = false;
	bool isContainedTrack = false;
	bool isUnContainedTrack = false;
	bool isContainedShower = false;
	bool isUnContainedShower = false;

        if(clusterList.size() !=0) {
	  pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
	  HitChargeVector hitChargeVector;
	  CartesianVector endPoint1(0.0, 0.0, 0.0);
	  CartesianVector endPoint2(0.0, 0.0, 0.0);
	  float trackLength = 0.f;
	  this->AddToSlidingFitCache(pCluster);
	  this->FillHitChargeVector(pCluster, hitChargeVector, endPoint1, endPoint2);
	  this->GetTrackLength(hitChargeVector, trackLength);
	  std::cout << "trackLength (type): " << trackLength << " cm" << std::endl;

	  const float upperY((endPoint1.GetY() > endPoint2.GetY()) ? endPoint1.GetY() : endPoint2.GetY());
	  const float lowerY((endPoint1.GetY() < endPoint2.GetY()) ? endPoint1.GetY() : endPoint2.GetY());

	  const float zAtUpperY((endPoint1.GetY() > endPoint2.GetY()) ? endPoint1.GetZ() : endPoint2.GetZ());
	  const float zAtLowerY((endPoint1.GetY() < endPoint2.GetY()) ? endPoint1.GetZ() : endPoint2.GetZ());

	  const float xAtUpperY((endPoint1.GetY() > endPoint2.GetY()) ? endPoint1.GetX() : endPoint2.GetX());
	  const float xAtLowerY((endPoint1.GetY() < endPoint2.GetY()) ? endPoint1.GetX() : endPoint2.GetX());
	  

	  isContained = ((upperY < parentMaxY - 1) && (upperY > parentMinY + 1) &&
		      (lowerY < parentMaxY - 1) && (lowerY > parentMinY + 1) &&
		      (zAtUpperY < parentMaxZ - 1) && (zAtUpperY > parentMinZ + 1) &&
		      (zAtLowerY < parentMaxZ - 1) && (zAtLowerY > parentMinZ + 1) &&
		      (xAtUpperY < parentMaxX - 1) && (xAtUpperY > parentMinX + 1) &&
		      (xAtLowerY < parentMaxX - 1) && (xAtLowerY > parentMinX + 1 ));


	  if (isContained == true && pPfo->GetParticleId() == 13) {
	    isContainedTrack = true;
	  }
	  else if (isContained == false && pPfo->GetParticleId() == 13) {
	    isUnContainedTrack = true;
	  }
	  else if (isContained == true && pPfo->GetParticleId() == 11) {
	    isContainedShower = true;
	  }
	  else if (isContained == false && pPfo->GetParticleId() == 11) {
	    isUnContainedShower = true;
	  } 

	  std::cout << "   " << std::endl;
	  std::cout << "isContained track? : " << isContainedTrack << std::endl;
	  std::cout << "isUnContained track? : " << isUnContainedTrack << std::endl;
	  std::cout << "isContained shower? : " << isContainedShower << std::endl;
	  std::cout << "isUnContained shower? : " << isUnContainedShower << std::endl;
	  std::cout << "   " << std::endl;
	}

	//-------Get the W cluster---------------------------------------------------------------------------------------

	std::cout << "clusterList.size()  " << clusterList.size() << std::endl;
	int r = 0;
	for (auto c :clusterList) {
	  if (r == 2 && clusterList.size() == 4) {
	    pCluster = c;
	  }
	  r++;
	}


	//-------Send it through the TDT -------------------------------------------------------------------------------
	//	if ((isContainedTrack == true  || isUnContainedTrack == true) && clusterList.size() == 4) {
	if (clusterList.size() == 4) {
	  if ((isContainedTrack == true  || isUnContainedTrack == true) && clusterList.size() == 4) {

	    if(clusterList.size() !=0) {
	      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

	      HitChargeVector hitChargeVector;
	      CartesianVector endPoint1(0.0, 0.0, 0.0);
	      CartesianVector endPoint2(0.0, 0.0, 0.0);
	      float trackLength = 0.f;
	      this->AddToSlidingFitCache(pCluster);
	      this->FillHitChargeVector(pCluster, hitChargeVector, endPoint1, endPoint2);
	      this->GetTrackLength(hitChargeVector, trackLength);
	      std::cout << "trackLength: (TDT) " << trackLength << " cm" << std::endl;


	      //------set lookup tables and fit parameters----------------------------------------------------------------------------
	      lar_content::NDQualitiesAlg::LookupTable lookupTableMuon;
	      lar_content::NDQualitiesAlg::LookupTable lookupTableProton;
	      lar_content::NDQualitiesAlg::LookupTable lookupTableElectron;
	      double InitialEnergy = 5000;
	      double BinWidth = 1;
	      double Mm = 105.66;               //now this is MeV 
	      double Mp = 938.28;
	      double Me = 0.511;

	    if (lookupTableMuon.GetMap().empty()){
	      lookupTableMuon.SetInitialEnergy(InitialEnergy); 
	      lookupTableMuon.SetBinWidth(BinWidth);

	      this->FillLookupTable(lookupTableMuon, Mm); 
	    }

	    if (lookupTableProton.GetMap().empty()){
	      lookupTableProton.SetInitialEnergy(InitialEnergy); 
	      lookupTableProton.SetBinWidth(BinWidth);

	      this->FillLookupTable(lookupTableProton, Mp); 
	    }


	    if (lookupTableElectron.GetMap().empty()){
	      lookupTableElectron.SetInitialEnergy(InitialEnergy); 
	      lookupTableElectron.SetBinWidth(BinWidth);

	      this->FillLookupTable(lookupTableElectron, Me); 
	    }
	    
	    double TotalCharge = (0.f);
	    double TotalHitWidth = (0.f);
	    double TotalHitEnergy = (0.f);
	    int hitsize = hitChargeVector.size();

	    std::cout << "  " << std::endl;
	    std::cout << "unbinned hit charge vector size : " << hitsize << std::endl;
	    for (HitCharge &hitCharge : hitChargeVector)
	      {
		TotalHitWidth += hitCharge.GetHitWidth();
		TotalCharge += hitCharge.GetCharge();
		TotalHitEnergy += hitCharge.GetCaloHit()->GetInputEnergy();
	      }

	    double maxScaleM(lookupTableMuon.GetMaxRange()/trackLength);   //maxScale = maxrange/tracklength
	    if (maxScaleM < 1.0) {
	      maxScaleM = 1.1;
	    }
	    double maxScaleP(lookupTableProton.GetMaxRange()/trackLength);
	    if (maxScaleP < 1.0) {
	      maxScaleP = 1.1;
	    }
	    double maxScaleE(lookupTableElectron.GetMaxRange()/trackLength);
	    if (maxScaleE < 1.0) {
	      maxScaleE = 1.1;
	    }
	    std::cout << "lookupTableM.GetMaxRange() " << lookupTableMuon.GetMaxRange() << std::endl;
	    std::cout << "lookupTableP.GetMaxRange() " << lookupTableProton.GetMaxRange() << std::endl;
	    std::cout << "lookupTableE.GetMaxRange() " << lookupTableElectron.GetMaxRange() << std::endl;

	    lar_content::NDQualitiesAlg::HitChargeVector binnedHitChargeVector;
	    BinHitChargeVector(hitChargeVector, binnedHitChargeVector);
	    if ( binnedHitChargeVector.size() == 0) {
	      binnedHitChargeVector = hitChargeVector;
	    }  

	    int binnedsize = binnedHitChargeVector.size();
	    std::cout << "binned hit charge vector size : " << binnedsize << std::endl;
	    std::cout << "  " << std::endl;

	    //--------------------forwards fit--------------------------------------------------------------------------------
	    const int nParameters = 3;
	    const std::string parName[nParameters]   = {"ENDENERGY", "SCALE", "EXTRA"};
	    const double vstart[nParameters] = {2.1, 1.0, 1.0};
	    std::list<double> chisquaredlist;
	    std::vector<double> dEdx_2D_av_list;
	    std::vector<double> dEdx_2D_raw_av_list;
	    std::vector<double> dEdx_2D_tot_list;
	    std::vector<double> dEdx_2D_raw_tot_list;
	    //  std::vector<double> diff_raw_list;
	    std::vector<double> p0list;
	    std::vector<double> p1list;
	    std::vector<double> p2list;
	    std::vector<int> countlist;

	    double dEdx_2D_av_f = 0.0;
	    double dEdx_2D_raw_av_f = 0.0;
	    // double dEdx_2D_raw_av_f_diff = 0.0;
	    double dEdx_dE_total_f = 0.0;
	    double dEdx_dE_total_raw_f = 0.0;
	    int forwardscount = 0;
	    // double changeval = 0.0

	    if (isContainedTrack == true) {                                            //if a contained track
	      const double step[nParameters] = {1, 100, 1};
	      const double highphysbound[nParameters] = {25, maxScaleP, 5};

	      for (double p0=vstart[0]; p0 < highphysbound[0];p0 = p0 + step[0]){
		for (double p1=vstart[1]; p1 < highphysbound[1];p1 = p1 + step[1]){
		  for (double p2=vstart[2]; p2 < highphysbound[2];p2 = p2 + step[2]){
		    double Ee(p0), L(p1 * trackLength); //energy, length
		    double Le(GetLengthfromEnergy(lookupTableProton, Ee));  //length
		    double Ls(Le - L);
		    double Es(GetEnergyfromLength(lookupTableProton, Ls));    //energy
		    double alpha((Es - Ee)/TotalCharge), beta(L/TotalHitWidth);
    
		    double chisquared(0.0);
		    double dEdx_2D_av = 0.0;
		    double dEdx_2D_tot = 0.0;
		    double dEdx_2D_raw_av = 0.0;
		    double dEdx_2D_raw_tot = 0.0;

		    int count = 0;
		    //  double diff_raw_av = 0.0;
		    //  double previous = 0.0;
		    //minimise this bit-----
		    // std::cout << "  " << std::endl;
		    for (lar_content::NDQualitiesAlg::HitCharge hitCharge : binnedHitChargeVector) //for each hit charge in the binned vector
		      {
			double L_i(Ls + (p1 * hitCharge.GetLongitudinalPosition())); //length
			double E_i(GetEnergyfromLength(lookupTableProton, L_i));       //energy

			double dEdx_2D(p2 * (beta/alpha) * BetheBloch(E_i, Mp));   //energy per length, calculated
			double dEdx_2D_raw((beta/alpha) * BetheBloch(E_i, Mp));
			double ChargeOverWidth(hitCharge.GetChargeOverWidth());  //energy per length, experimental
			//	std::cout << "charge over width f" << ChargeOverWidth << "  dEdx_2D_raw " << dEdx_2D_raw <<std::endl;

			if (dEdx_2D_raw > 0.0) {
			  dEdx_2D_raw_tot = dEdx_2D_raw_tot + dEdx_2D_raw;
			}

			if (dEdx_2D > 0.0) {
			  dEdx_2D_tot = dEdx_2D_tot + dEdx_2D;
			}
			
			//	if (dEdx_2D_raw < 5.0 && dEdx_2D_raw > 0.0) {
			dEdx_2D_raw_av = dEdx_2D_raw_av + ChargeOverWidth;
			count = count + 1;
			//	}
			if (dEdx_2D < 5.0 && dEdx_2D > 0.0) {
			  dEdx_2D_av = dEdx_2D_av + dEdx_2D;
			}

			chisquared += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
		      }
		    // std::cout << "  " << std::endl;
		    //-----------
		    if (chisquaredlist.size() > 0) {
		      if (chisquared < chisquaredlist.back()){
			chisquaredlist.push_back(chisquared);
			dEdx_2D_av_list.push_back(dEdx_2D_av);
			dEdx_2D_raw_av_list.push_back(dEdx_2D_raw_av);
			dEdx_2D_tot_list.push_back(dEdx_2D_tot);
			dEdx_2D_raw_tot_list.push_back(dEdx_2D_raw_tot);
			//  diff_raw_list.push_back(diff_raw_av);
			p0list.push_back(p0);
			p1list.push_back(p1);
			p2list.push_back(p2);
			countlist.push_back(count);
			chisquaredlist.pop_front();
			dEdx_2D_av_list.erase(dEdx_2D_av_list.begin());
			dEdx_2D_raw_av_list.erase(dEdx_2D_raw_av_list.begin());
			dEdx_2D_tot_list.erase(dEdx_2D_tot_list.begin());
			dEdx_2D_raw_tot_list.erase(dEdx_2D_raw_tot_list.begin());
			p0list.erase(p0list.begin());
			p1list.erase(p1list.begin());
			p2list.erase(p2list.begin());
			countlist.erase(countlist.begin());
		      }
		    }
		    else {
		      chisquaredlist.push_back(chisquared);
		      dEdx_2D_av_list.push_back(dEdx_2D_av);
		      dEdx_2D_raw_av_list.push_back(dEdx_2D_raw_av);
		      dEdx_2D_tot_list.push_back(dEdx_2D_tot);
		      dEdx_2D_raw_tot_list.push_back(dEdx_2D_raw_tot);
		      p0list.push_back(p0);
		      p1list.push_back(p1);
		      p2list.push_back(p2);
		      countlist.push_back(count);
		    }
		  }
		}
	      }
	    }
	    else if (isUnContainedTrack == true || isUnContainedShower == true || isContainedShower == true) {                                          //if an uncontained track, find end energy
	      double useScale;
	      LookupTable useLookup;
	      double useMass;

	      if (isUnContainedTrack == true) {
		useScale = maxScaleM;
		useLookup = lookupTableMuon;
		useMass = Mm;
	      }
	      else {
		useScale = maxScaleE;
		useLookup = lookupTableElectron;
		useMass = Me;
	      }

	      const double step[nParameters] = {1, 100, 1};
	      const double highphysbound[nParameters] = {5000, useScale, 5};

	      int countb = 0;
	      for (double p0=vstart[0]; p0 < highphysbound[0];p0 = p0 + step[0]){
		for (double p1=vstart[1]; p1 < highphysbound[1];p1 = p1 + step[1]){
		  for (double p2=vstart[2]; p2 < highphysbound[2];p2 = p2 + step[2]){
		    double Ee(p0), L(p1 * trackLength); //energy, length
		    double Le(GetLengthfromEnergy(useLookup, Ee));  //length
		    double Ls(Le - L);
		    double Es(GetEnergyfromLength(useLookup, Ls));    //energy
		    double alpha((Es - Ee)/TotalCharge), beta(L/TotalHitWidth);
    
		    double chisquared(0.0);
		    double dEdx_2D_av = 0.0;
		    double dEdx_2D_raw_av = 0.0;
		    int count = 0;
		    front_f_dE = 0.0;

		    // double diff_raw_av = 0.0;
		    //   double previous = 0.0;
		    //minimise this bit-----
		    // std::cout << "  " << std::endl;
		    for (lar_content::NDQualitiesAlg::HitCharge hitCharge : binnedHitChargeVector) //for each hit charge in the binned vector
		      {
			double L_i(Ls + (p1 * hitCharge.GetLongitudinalPosition())); //length
			double E_i(GetEnergyfromLength(useLookup, L_i));       //energy

			double dEdx_2D(p2 * (beta/alpha) * BetheBloch(E_i, useMass));   //energy per length, calculated
			//	double dEdx_2D_raw((beta/alpha) * BetheBloch(E_i, useMass));
			double ChargeOverWidth(hitCharge.GetChargeOverWidth());  //energy per length, experimental
			//	std::cout << "charge over width f " << ChargeOverWidth << " dEdx_2D_raw  " << dEdx_2D_raw <<std::endl;

			//	if (dEdx_2D_raw > 0.0) {
			dEdx_2D_raw_av = dEdx_2D_raw_av + ChargeOverWidth;
			count = count + 1;
			//	}
			if (dEdx_2D > 0.0) {
			  dEdx_2D_av = dEdx_2D_av + dEdx_2D;
			}

			if (isContainedShower == true || isUnContainedShower == true) {
			  front_f_dE = front_f_dE + ChargeOverWidth;
			}
			//	else if (isUnContainedShower == true && (countb == 0 || countb == 1 || countb == 2)) { 
			//	  front_f_dE = front_f_dE + ChargeOverWidth;
			//	}

			chisquared += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
			countb++;
		      }
		    // std::cout << "  " << std::endl;
		    //-----------
		    if (chisquaredlist.size() > 0) {
		      if (chisquared < chisquaredlist.back()){
			chisquaredlist.push_back(chisquared);
			dEdx_2D_av_list.push_back(dEdx_2D_av);
			dEdx_2D_raw_av_list.push_back(dEdx_2D_raw_av);
			p0list.push_back(p0);
			p1list.push_back(p1);
			p2list.push_back(p2);
			countlist.push_back(count);
			chisquaredlist.pop_front();
			dEdx_2D_av_list.erase(dEdx_2D_av_list.begin());
			dEdx_2D_raw_av_list.erase(dEdx_2D_raw_av_list.begin());
			p0list.erase(p0list.begin());
			p1list.erase(p1list.begin());
			p2list.erase(p2list.begin());
			countlist.erase(countlist.begin());
		      }
		    }
		    else {
		      chisquaredlist.push_back(chisquared);
		      dEdx_2D_av_list.push_back(dEdx_2D_av);
		      dEdx_2D_raw_av_list.push_back(dEdx_2D_raw_av);
		      // diff_raw_list.push_back(diff_raw_av);
		      p0list.push_back(p0);
		      p1list.push_back(p1);
		      p2list.push_back(p2);
		      countlist.push_back(count);
		    }
		  }
		}
	      }
	    }

	    //-------Get the minimum chi squared--------------------------------------------------------------------------------------
	    //  std::list<double>::iterator minchisquared = std::min_element(chisquaredlist.begin(), chisquaredlist.end());
	    // auto index = std::distance(chisquaredlist.begin(), minchisquared);
	    // std::cout << "index! " << index  << std::endl;
	    //std::cout << "size! " << chisquaredlist.size()  << std::endl;

	    //  int indexvalue = index;
	    int indexvalue = 0;
	    double outpar[nParameters];

	    outpar[0] = p0list[indexvalue];
	    outpar[1] = p1list[indexvalue];
	    outpar[2] = p2list[indexvalue];


	    std::cout << "  " << std::endl;
	    std::cout << "0 = : " << outpar[0] << std::endl;
	    std::cout << "1 = : " << outpar[1] << std::endl;
	    std::cout << "2 = : " << outpar[2] << std::endl;

	    //  std::cout << "----------------------" << std::endl;
	    //  t = std::clock() - t;
	    //  std::printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	    //  time(&now2);
	    //  seconds = std::difftime(now2,now);
	    //  std::printf ("%.f seconds to run.\n", seconds);
	    //  std::cout << "----------------------" << std::endl;
	    
	    //-------trying to get the dEdx--------------------------------------------------------------------------------------------
	    // if(isUnContainedTrack == true) {
	    double forwardsp[5];

	    forwardsp[0] = dEdx_2D_av_list[indexvalue];
	    forwardsp[1] = dEdx_2D_raw_av_list[indexvalue];
	    // forwardsp[2] = diff_raw_list[indexvalue];


	    if(isContainedTrack == true) {
	      forwardsp[2] = dEdx_2D_raw_tot_list[indexvalue];
	      forwardsp[3] = dEdx_2D_tot_list[indexvalue];

	      dEdx_dE_total_f = forwardsp[3];
	      dEdx_dE_total_raw_f = forwardsp[2];
	    }
	     
	    forwardsp[4] = countlist[indexvalue];
	    forwardscount = forwardsp[4];

	    if (isContainedShower == true || isUnContainedShower == true) {
	      front_f_dE = front_f_dE/binnedHitChargeVector.size();
	    }
	      

	    dEdx_2D_av_f = forwardsp[0]/forwardscount;
	    dEdx_2D_raw_av_f = forwardsp[1]/forwardscount;

	    std::cout << "forwardscount = " << forwardscount << std::endl;
	    std::cout << "dEdx_2D_raw = " << forwardsp[1] << std::endl;
	      
	    //  }


	    //-----------------------------backwards fit------------------------------------------------------------------------------------
	      //  int t2;
	      // t2 = std::clock();
	      //   // time_t now3;
	     // time_t now4;
		// double seconds2;
	      //  time(&now3);

	    const int nParameters2 = 3;
	    const std::string parName2[nParameters2]   = {"ENDENERGY", "SCALE", "EXTRA"};
	    const double vstart2[nParameters2] = {2.1, 1.0, 0.0};
	    std::list<double> chisquaredlist2;
	    std::vector<double> dEdx_2D_av_list_2;
	    std::vector<double> dEdx_2D_raw_av_list_2;
	    std::vector<double> dEdx_2D_tot_list_2;
	    std::vector<double> dEdx_2D_raw_tot_list_2;
	    //  std::vector<double> diff_raw_list_2;
	    std::vector<double> p0list2;
	    std::vector<double> p1list2;
	    std::vector<double> p2list2;
	    std::vector<int> countlist2;


	    double dEdx_2D_av_b = 0.0;
	    double dEdx_2D_raw_av_b = 0.0;
	    // double dEdx_2D_raw_av_b_diff = 0.0;
	    double dEdx_dE_total_b = 0.0;
	    double dEdx_dE_total_raw_b = 0.0;
	    int backwardscount = 0;


	    if (isContainedTrack == true) {                                                              //if a contained track
	      const double step2[nParameters2] = {1, 100, 1};
	      const double highphysbound2[nParameters2] = {25, maxScaleP, 5};
	      for (double p02=vstart2[0];  p02 < highphysbound2[0];p02 = p02 + step2[0]){
		for (double p12=vstart2[1]; p12 < highphysbound2[1];p12 = p12 + step2[1]){
		  for (double p22=vstart2[2];  p22 < highphysbound2[2];p22 = p22 + step2[2]){
		    double Ee(p02), L(p12 * trackLength); //energy, length
		    double Le(GetLengthfromEnergy(lookupTableProton, Ee));  //length
		    double Ls(Le - L);
		    double Es(GetEnergyfromLength(lookupTableProton, Ls));    //energy

		    double alpha((Es - Ee)/TotalCharge), beta(L/TotalHitWidth);
	     
		    double chisquared2(0.0);
		    double dEdx_2D_av_2 = 0.0;
		    double dEdx_2D_raw_av_2 = 0.0;
		    double dEdx_2D_tot_2 = 0.0;
		    double dEdx_2D_raw_tot_2 = 0.0;
		    int count2 = 0;

		    //minimise this bit-----
		    // std::cout << "   " << std::endl;
		    for (lar_content::NDQualitiesAlg::HitCharge hitCharge : binnedHitChargeVector) //for each hit charge in the binned vector
		      {
			double L_i(Ls + (p12 * hitCharge.GetLongitudinalPosition())); //length
			double E_i(GetEnergyfromLength(lookupTableProton, L_i));       //energy

			double dEdx_2D(p22 * (beta/alpha) * BetheBloch(E_i, Mp));   //energy per length, calculated
			double dEdx_2D_raw((beta/alpha) * BetheBloch(E_i, Mp));
			double ChargeOverWidth(hitCharge.GetChargeOverWidth());  //energy per length, experimental
			//	std::cout << "charge over width b " << ChargeOverWidth << " dEdx_2D_raw  " << dEdx_2D_raw <<std::endl;

			if (dEdx_2D_raw > 0.0) {
			  dEdx_2D_raw_tot_2 = dEdx_2D_raw_tot_2 + dEdx_2D_raw;
			}
			if (dEdx_2D > 0.0) {
			  dEdx_2D_tot_2 = dEdx_2D_tot_2 + dEdx_2D;
			}
			
			//	if (dEdx_2D_raw < 5.0 && dEdx_2D_raw > 0.0) {
			dEdx_2D_raw_av_2 = dEdx_2D_raw_av_2 + ChargeOverWidth;
			count2 = count2 + 1;
			//	}
			if (dEdx_2D < 5.0 && dEdx_2D > 0.0) {
			  dEdx_2D_av_2 = dEdx_2D_av_2 + dEdx_2D;
			}


			chisquared2 += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
		      }
		    // std::cout << "  " << std::endl;
		    //-----------
		    if (chisquaredlist2.size() > 0) {
		      if (chisquared2 < chisquaredlist2.back()){
			chisquaredlist2.push_back(chisquared2);
			dEdx_2D_av_list_2.push_back(dEdx_2D_av_2);
			dEdx_2D_raw_av_list_2.push_back(dEdx_2D_raw_av_2);
			dEdx_2D_tot_list_2.push_back(dEdx_2D_tot_2);
			dEdx_2D_raw_tot_list_2.push_back(dEdx_2D_raw_tot_2);
			p0list2.push_back(p02);
			p1list2.push_back(p12);
			p2list2.push_back(p22);
			countlist2.push_back(count2);
			chisquaredlist2.pop_front();
			dEdx_2D_av_list_2.erase(dEdx_2D_av_list_2.begin());
			dEdx_2D_raw_av_list_2.erase(dEdx_2D_raw_av_list_2.begin());
			dEdx_2D_tot_list_2.erase(dEdx_2D_tot_list_2.begin());
			dEdx_2D_raw_tot_list_2.erase(dEdx_2D_raw_tot_list_2.begin());
			p0list2.erase(p0list2.begin());
			p1list2.erase(p1list2.begin());
			p2list2.erase(p2list2.begin());
			countlist2.erase(countlist2.begin());
		      }
		    }
		    else {
		      chisquaredlist2.push_back(chisquared2);
		      dEdx_2D_av_list_2.push_back(dEdx_2D_av_2);
		      dEdx_2D_raw_av_list_2.push_back(dEdx_2D_raw_av_2);
		      dEdx_2D_tot_list_2.push_back(dEdx_2D_tot_2);
		      dEdx_2D_raw_tot_list_2.push_back(dEdx_2D_raw_tot_2);
		      p0list2.push_back(p02);
		      p1list2.push_back(p12);
		      p2list2.push_back(p22);
		      countlist2.push_back(count2);
		    }
		  }
		}
	      }

	    }
	    else if (isUnContainedTrack == true || isUnContainedShower == true || isContainedShower == true ) {                                 //if it's an uncontained track
	      double useScale;
	      LookupTable useLookup;
	      double useMass;

	      if (isUnContainedTrack == true) {
		useScale = maxScaleM;
		useLookup = lookupTableMuon;
		useMass = Mm;
	      }
	      else {
		useScale = maxScaleE;
		useLookup = lookupTableElectron;
		useMass = Me;
	      }


	      const double  step2[nParameters2] = {1, 100, 1};
	      const double highphysbound2[nParameters2] = {5000, useScale, 5};

	      int countb = 0;
	      for (double p02=vstart2[0];  p02 < highphysbound2[0];p02 = p02 + step2[0]){
		for (double p12=vstart2[1]; p12 < highphysbound2[1];p12 = p12 + step2[1]){
		  for (double p22=vstart2[2];  p22 < highphysbound2[2];p22 = p22 + step2[2]){
		    double Ee(p02), L(p12 * trackLength); //energy, length
		    double Le(GetLengthfromEnergy(useLookup, Ee));  //length
		    double Ls(Le - L);
		    double Es(GetEnergyfromLength(useLookup, Ls));    //energy

		    double alpha((Es - Ee)/TotalCharge), beta(L/TotalHitWidth);
	     
		    double chisquared2(0.0);
		    double dEdx_2D_av_2 = 0.0;
		    double dEdx_2D_raw_av_2 = 0.0;
		    int count2 = 0;
		    front_b_dE = 0.0;
		    //minimise this bit-----
		    // std::cout << "    " << std::endl;
		    for (lar_content::NDQualitiesAlg::HitCharge hitCharge : binnedHitChargeVector) //for each hit charge in the binned vector
		      {
			double L_i(Ls + (p12 * hitCharge.GetLongitudinalPosition())); //length
			double E_i(GetEnergyfromLength(useLookup, L_i));       //energy

			double dEdx_2D(p22 * (beta/alpha) * BetheBloch(E_i, useMass));   //energy per length, calculated
			//	double dEdx_2D_raw((beta/alpha) * BetheBloch(E_i, useMass));
			double ChargeOverWidth(hitCharge.GetChargeOverWidth());  //energy per length, experimental
			//	std::cout << "charge over width f " << ChargeOverWidth << " dEdx_2D_raw  " << dEdx_2D_raw << std::endl;
		       
			//	if (dEdx_2D_raw > 0.0) {
			dEdx_2D_raw_av_2 = dEdx_2D_raw_av_2 + ChargeOverWidth;
			count2 = count2 + 1;
			//	}
			if (dEdx_2D > 0.0) {
			  dEdx_2D_av_2 = dEdx_2D_av_2 + dEdx_2D;
			}

        
			if (isContainedShower == true || isUnContainedShower == true) {
			  front_b_dE = front_b_dE + ChargeOverWidth;
			}
			//	//	else if (isUnContainedShower == true && (countb == 0 || countb == 1 || countb == 2)) { 
			  //	  front_b_dE = front_b_dE + ChargeOverWidth;
			  //	}
			
			chisquared2 += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
			countb++;
		      }
		    // std::cout << "    " << std::endl;
		    //-----------
		    if (chisquaredlist2.size() > 0) {
		      if (chisquared2 < chisquaredlist2.back()){
			chisquaredlist2.push_back(chisquared2);
			dEdx_2D_av_list_2.push_back(dEdx_2D_av_2);
			dEdx_2D_raw_av_list_2.push_back(dEdx_2D_raw_av_2);
			//  diff_raw_list_2.push_back(diff_raw_av_2);
			p0list2.push_back(p02);
			p1list2.push_back(p12);
			p2list2.push_back(p22);
			countlist2.push_back(count2);
			chisquaredlist2.pop_front();
			dEdx_2D_av_list_2.erase(dEdx_2D_av_list_2.begin());
			dEdx_2D_raw_av_list_2.erase(dEdx_2D_raw_av_list_2.begin());
			p0list2.erase(p0list2.begin());
			p1list2.erase(p1list2.begin());
			p2list2.erase(p2list2.begin());
			countlist2.erase(countlist2.begin());
		      }
		    }
		    else {
		      chisquaredlist2.push_back(chisquared2);
		      dEdx_2D_av_list_2.push_back(dEdx_2D_av_2);
		      dEdx_2D_raw_av_list_2.push_back(dEdx_2D_raw_av_2);
		      p0list2.push_back(p02);
		      p1list2.push_back(p12);
		      p2list2.push_back(p22);
		      countlist2.push_back(count2);
		    }
		  }
		}
	      }
	    }


	    //--------------find best backwards chi sqauared------------------------------------------------------------------------
	    // std::list<double>::iterator minchisquared2 = std::min_element(chisquaredlist2.begin(), chisquaredlist2.end());
	    // auto index2 = std::distance(chisquaredlist2.begin(), minchisquared2);   

	    // int indexvalue2 = index2;
	    int indexvalue2 = 0;
	    // std::cout << "index2 " << index2 << std::endl;

	    double outpar2[nParameters];

	    outpar2[0] = p0list2[indexvalue2];
	    outpar2[1] = p1list2[indexvalue2];
	    outpar2[2] = p2list2[indexvalue2];

	    std::cout << "0b = : " << outpar2[0] << std::endl;
	    std::cout << "1b = : " << outpar2[1] << std::endl;
	    std::cout << "2b = : " << outpar2[2] << std::endl;
	    std::cout << "  " << std::endl;

	    //  std::cout << "----------------------" << std::endl;
	    //  t2 = std::clock() - t2;
	    // std::printf ("It took me %d clicks (%f seconds).\n",t2,((float)t2)/CLOCKS_PER_SEC);
	    //time(&now4);
	    //seconds2 = std::difftime(now4,now3);
	    //std::printf ("%.f seconds to run.\n", seconds2);
	    //std::cout << "----------------------" << std::endl;


	    //---------------finding backwards dEdx---------------------------------------------------------------------------------
	    // if(isUnContainedTrack == true) {
	    double backwardsp[5];

	    backwardsp[0] = dEdx_2D_av_list_2[indexvalue2];
	    backwardsp[1] = dEdx_2D_raw_av_list_2[indexvalue2];
	    //  backwardsp[2] = diff_raw_list_2[indexvalue];

	    if(isContainedTrack == true) {
	      backwardsp[2] = dEdx_2D_tot_list_2[indexvalue2];
	      backwardsp[3] = dEdx_2D_raw_tot_list_2[indexvalue2];

	      dEdx_dE_total_b = backwardsp[2];
	      dEdx_dE_total_raw_b = backwardsp[3];
	    }
	  
	    forwardsp[4] = countlist2[indexvalue];
	    backwardscount = forwardsp[4];
	    
	    if (isContainedShower == true || isUnContainedShower == true) {
	      front_b_dE = front_b_dE/binnedHitChargeVector.size();
	    }

	    dEdx_2D_av_b = backwardsp[0]/backwardscount;
	    dEdx_2D_raw_av_b = backwardsp[1]/backwardscount;


	    std::cout << "backwardscount = " << backwardscount  << std::endl;
	    std::cout << "dEdx_2D_raw = " << backwardsp[1]  << std::endl;


	    //-----------get forawrds output details--------------------------------------------------------------------
	    double f_Ee(outpar[0]), f_L(outpar[1] * trackLength);
	    double andyenergy = 2.99817*0.0001*pow(trackLength, 1.26749) + 1.20364*0.001*pow(trackLength, 5.45196*0.1);
	    // std::cout << "andyenergy = " << andyenergy << std::endl;
	  
	    double f_Le = 0;
	    if (isContainedTrack == true) {
	      f_Le = (GetLengthfromEnergy(lookupTableProton, f_Ee));
	    }
	    else if (isUnContainedTrack == true) {
	      f_Le = (GetLengthfromEnergy(lookupTableMuon, f_Ee));
	      //  std::cout << "f_Le " << f_Le << std::endl;
	    }
	    else if (isUnContainedShower == true || isContainedShower == true) {
	      f_Le = (GetLengthfromEnergy(lookupTableElectron, f_Ee));
	      //  std::cout << "f_Le " << f_Le << std::endl;
	    }

	    double f_Ls = f_Le - f_L;                                   //this is track *start* point, f_L will be track Length
	    // std::cout << "f_Ee: (ending energy) " << f_Ee << std::endl;               //this is *end* energy
	    double andyenergy_f_L = pow(2.99817*0.0001*f_L, 1.26749) + pow(1.20364*0.001*f_L, 5.45196*0.1);
	    // std::cout << "andyenergy_f_L = " << andyenergy_f_L << std::endl;

	    double f_Es = 0;
	    if (isContainedTrack == true) {
	      f_Es = GetEnergyfromLength(lookupTableProton, f_Ls);
	    }
	    else if (isUnContainedTrack == true) {
	      f_Es = GetEnergyfromLength(lookupTableMuon, f_Ls);
	      std::cout << "f_Es " << f_Es << std::endl;
	    }
	    else if (isUnContainedShower == true || isContainedShower == true) {
	      f_Es = GetEnergyfromLength(lookupTableElectron, f_Ls);
	      //  std::cout << "f_Le " << f_Le << std::endl;
	    }
	    //e = end, i = start, s = ?, maybe here s = i 

	    double f_deltaE = f_Es - f_Ee;
	    // std::cout << "f_deltaE " << f_deltaE << std::endl;

	    double f_alpha = f_deltaE/TotalCharge;
	    double f_beta = f_L/TotalHitWidth;


	    //--------get backwards output details-----------------------------------------------------------------------------------
	    double b_Ee(outpar2[0]), b_L(outpar2[1] * trackLength);
	    double andyenergy_b_L = 2.99817*0.0001*pow(f_L, 1.26749) + 1.20364*0.001*pow(f_L, 5.45196*0.1);
	    //  std::cout << "andyenergy_b_L = " << andyenergy_b_L << std::endl;

	    double b_Le = 0;
	    if (isContainedTrack == true) {
	      b_Le = (GetLengthfromEnergy(lookupTableProton, b_Ee));
	    }
	    else if (isUnContainedTrack == true) {
	      b_Le = (GetLengthfromEnergy(lookupTableMuon, b_Ee));
	      //  std::cout << "b_Le " << b_Le << std::endl;
	    }
	    else if (isUnContainedShower == true || isContainedShower == true) {
	      b_Le = (GetLengthfromEnergy(lookupTableElectron, b_Ee));
	      //  std::cout << "f_Le " << f_Le << std::endl;
	    }
	    double b_Ls = b_Le - b_L;
	    //  std::cout << "b_Ee: (ending energy) " << b_Ee << std::endl;

	    double b_Es = 0;
	    if (isContainedTrack == true) {
	      b_Es = GetEnergyfromLength(lookupTableProton, b_Ls);
	    }
	    else if (isUnContainedTrack == true) {
	      b_Es = GetEnergyfromLength(lookupTableMuon, b_Ls);
	      //  std::cout << "b_Es " << b_Es << std::endl;
	    }
	    else if (isUnContainedShower == true || isContainedShower == true) {
	      b_Es = GetEnergyfromLength(lookupTableElectron, b_Ls);
	      //  std::cout << "f_Le " << f_Le << std::endl;
	    }
	    double b_deltaE = b_Es - b_Ee;
	    // std::cout << "b_deltaE " << b_deltaE << std::endl;
	  

	    double b_alpha = b_deltaE/TotalCharge;
	    double b_beta = b_L/TotalHitWidth;

	    //------------------------------grab variables from min chi squared----------------------------------------------

	    fdeltaEfromminchi = fdeltaEfromminchi +  f_deltaE;
	    bdeltaEfromminchi = bdeltaEfromminchi +  b_deltaE;



	    //--------------------------go though hit charge vector to get variables at each dQdx--------------------------
	    float i = 0.0;
	    double energybylength = 0;
	    double energybylengthb = 0;
	    double energybylength_end = 0;
	    double energybylengthb_end = 0;
	    double addingenergybylength = 0;

	    int nHitsConsidered(0);

	    double t_dEdx_total = 0.0;
	    double t_dEdx_ave = 0.0;
	    double f_dEdx_2D = 0.0;
	    double b_dEdx_2D = 0.0;
	    double f_dEdx_2D_average = 0.0;
	    double b_dEdx_2D_average = 0.0;
	    if (1 == 2) {
	      std::cout << t_dEdx_ave << std::endl;
	    }

	    for (HitCharge &hitCharge : hitChargeVector)
	      {

		if(isContainedTrack == true) {                                               //contained track
		  double f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());         //sum these? Gte the first one only
		  double f_E_i = GetEnergyfromLength(lookupTableProton, f_L_i);
		  f_dEdx_2D = outpar[2] * (f_beta/f_alpha) * BetheBloch(f_E_i, Mp);
		  

		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "f_dEdx_2D", f_dEdx_2D ));
		  if (i == 0) {
		    energybylength = (f_dEdx_2D*(f_L/1))/1000;
		  }
		  addingenergybylength = f_dEdx_2D*hitCharge.GetHitWidth() + addingenergybylength;

		  double b_L_i = b_Ls + (outpar2[1] * (trackLength - hitCharge.GetLongitudinalPosition()));
		  double b_E_i = GetEnergyfromLength(lookupTableProton, b_L_i);
		  b_dEdx_2D = outpar2[2] * (b_beta/b_alpha) * BetheBloch(b_E_i, Mp);

		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "b_dEdx_2D", b_dEdx_2D ));
		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "i", i));
		  // PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttree"));
		}

		if(isUnContainedTrack == true || isUnContainedShower == true || isContainedShower == true) {                                                  //uncontained track
		  LookupTable useLookup;
		  if (isUnContainedTrack == true) {
		    useLookup = lookupTableMuon;
		  }
		  else {
		    useLookup = lookupTableElectron;
		  }
		  
		  double f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());         //sum these? Gte the first one only
		  double f_E_i = GetEnergyfromLength(useLookup, f_L_i);
		  //  std::cout << "F E I " << f_E_i  << std::endl;
		  f_dEdx_2D = outpar[2] * (f_beta/f_alpha) * BetheBloch(f_E_i, Mm);
		  //  std::cout << "f_dEdx_2D " << f_dEdx_2D << std::endl;
		  f_dEdx_2D_average = f_dEdx_2D_average + std::abs(f_dEdx_2D);
		 

		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "f_dEdx_2D", f_dEdx_2D ));
		  if (i == 0) {
		    energybylength = (f_dEdx_2D*(f_L/1))/1000;
		  }
		  addingenergybylength = f_dEdx_2D*hitCharge.GetHitWidth() + addingenergybylength;

		  double b_L_i = b_Ls + (outpar2[1] * (trackLength - hitCharge.GetLongitudinalPosition()));
		  double b_E_i = GetEnergyfromLength(useLookup, b_L_i);

		  b_dEdx_2D = outpar2[2] * (b_beta/b_alpha) * BetheBloch(b_E_i, Mm);
		  b_dEdx_2D_average = b_dEdx_2D_average + std::abs(b_dEdx_2D);

		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "b_dEdx_2D", b_dEdx_2D ));
		  // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "i", i));
		  //PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttree"));
		}


		//--------------------alternate calulation of dEdx-----------------------------------------------
		double t_dQdx = hitCharge.GetChargeOverWidth();
		double a = 0.93;
		double b = 0.562;
		double wion = 23.6;
		double e = 1.0;

		double t_dEdx = ((exp(b*(wion/e)*t_dQdx)) - a)/(b);
		//	std::cout << "t_dEdx_t" << t_dEdx << std::endl;
		t_dEdx_total = t_dEdx_total + t_dEdx;
		//	std::cout << "t_dEdx_total running" << t_dEdx_total << std::endl;


		if (i == 0) {
		  energybylengthb = (b_dEdx_2D*(b_L/0.5))/1000;                                                  
		}

		energybylength_end = (f_dEdx_2D*f_L*0.5)/1000;
		energybylengthb_end = (b_dEdx_2D*b_L*0.5)/1000;


		double Q_fit_f(f_dEdx_2D * hitCharge.GetHitWidth());
		double Q_fit_b(b_dEdx_2D * hitCharge.GetHitWidth());

		float forwardsDelta(hitCharge.GetChargeOverWidth() - f_dEdx_2D), backwardsDelta(hitCharge.GetChargeOverWidth() - b_dEdx_2D);

		float f_sigma(std::sqrt((0.00164585 * f_dEdx_2D * f_dEdx_2D) + (0.0201838 * std::abs(f_dEdx_2D)))); //80%
		float b_sigma(std::sqrt((0.00164585 * b_dEdx_2D * b_dEdx_2D) + (0.0201838 * std::abs(b_dEdx_2D)))); //80%
		//std::cout << "F DEDX " << f_dEdx_2D << std::endl;

		float lp(hitCharge.GetLongitudinalPosition()), hw(hitCharge.GetHitWidth());
		float f_Q_fit_f(Q_fit_f), f_Q_fit_b(Q_fit_b);
		HitCharge forwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_f, f_sigma);
		forwardsFitPoints.push_back(forwardsRecoHitCharge);
		HitCharge backwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_b, b_sigma);
		backwardsFitPoints.push_back(backwardsRecoHitCharge);

		float forwardsHitChisquared((forwardsDelta * forwardsDelta)/(f_sigma * f_sigma));
		float backwardsHitChisquared((backwardsDelta * backwardsDelta)/(b_sigma * b_sigma));
		//	std::cout << "F DELTA " << forwardsDelta << std::endl;
		//	std::cout << "F SIGMA " << f_sigma << std::endl;

		float Q_fit_forwards(Q_fit_f), Q_fit_backwards(Q_fit_b); 

		hitCharge.SetForwardsFitCharge(Q_fit_forwards); 
		hitCharge.SetForwardsSigma(f_sigma);
		hitCharge.SetForwardsDelta(forwardsDelta);
		hitCharge.SetForwardsChiSquared(forwardsHitChisquared);
		//std::cout << "F HIT CHI SQAURE" << forwardsChiSquared << std::endl;

		hitCharge.SetBackwardsFitCharge(Q_fit_backwards); 
		hitCharge.SetBackwardsSigma(b_sigma);
		hitCharge.SetBackwardsDelta(backwardsDelta);
		hitCharge.SetBackwardsChiSquared(backwardsHitChisquared);

		if (!((hitChargeVector.size() >= 2 * numberHitsToConsider) && nHitsConsidered > numberHitsToConsider && nHitsConsidered < hitChargeVector.size() - numberHitsToConsider))
		  {
		    forwardsChiSquared += forwardsHitChisquared;
		    backwardsChiSquared += backwardsHitChisquared;
		    //  std::cout << "in addition loop" << std::endl;
		    //  std::cout << "F CHI SQAURE loop" << forwardsChiSquared << std::endl;
		  }

		nHitsConsidered++;


		i++;
	      }


	    //-----------------calculated pfo energies--------------------------------------------------
	    std::cout << "   " << std::endl;
	    std::cout << "t_dEdx_total " << t_dEdx_total/i << std::endl;
	    t_dEdx_ave = t_dEdx_total/i;
	    // std::cout <<"width*i  " << BinWidth*i << "   i " << i <<  std::endl;
	    // std::cout <<"binwidth  " << BinWidth << std::endl;

	    f_dEdx_2D_average = f_dEdx_2D_average/i;
	    b_dEdx_2D_average = b_dEdx_2D_average/i;
	    std::cout << "f_dEdx_2D_average " << f_dEdx_2D_average << std::endl;
	    std::cout << "b_dEdx_2D_average " << b_dEdx_2D_average << std::endl;
	    //std::cout << "F CHI SQAURE " << forwardsChiSquared << std::endl;
	  

	    energypfo = energypfo + (f_Es/1000);
	    energycount = energycount + 1;
	    andyenergytotal = andyenergytotal + andyenergy;
	    andyenergytotal_f_L = andyenergytotal_f_L + andyenergy_f_L;
	    andyenergytotal_b_L = andyenergytotal_b_L + andyenergy_b_L;
	    


	    if (energybylength <= 10) {
	      energybylengthtotal = energybylengthtotal + energybylength;
	    }
	    if (energybylengthb <= 10) {
	      energybylengthtotalb = energybylengthtotalb + energybylengthb;
	    }

	    std::cout << "energy in this pfo by length end: " << energybylength_end << " GeV" << std::endl;
	    std::cout << "energy in this pfo by length back end: " << energybylengthb_end  << " GeV" << std::endl;

	    if (pPfo->GetNParentPfos() == 1) {
	      neutrinoenergy = neutrinoenergy + energypfo/energycount;
	    }

	    std::cout << "   " << std::endl;
	    //--------------------------get mc particle energy------------------------------------------------------
	    float mcenergyp = 0.0;
	    int id = 0;
	    try{
	      const MCParticle *pMCParticle = LArMCParticleHelper::GetMainMCParticle(pPfo);
	      const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
	      std::cout << " -----------------------" << std::endl;
	      std::cout << "pParentMCParticle energy  " << pParentMCParticle->GetEnergy() << std::endl;
	      std::cout << "is mc parent neutrino " << LArMCParticleHelper::IsNeutrino(pParentMCParticle) << std::endl;
	      std::cout << "mc nuance code " <<LArMCParticleHelper::GetNuanceCode(pMCParticle) << std::endl;
	      std::cout << pMCParticle << " <----- MC particle " << std::endl;
	      std::cout << pMCParticle->GetParticleId() << " <----- MC particle Id " << std::endl;
	      mcenergyp = pMCParticle->GetEnergy();
	      id = pMCParticle->GetParticleId();
	      std::cout << pMCParticle->GetEnergy() << " <----- MC particle energy " << std::endl;
	      std::cout << " -----------------------" << std::endl;
	    } catch (...){
	      std::cout << "lost one" << std::endl;
	    }


	    //-----------------------what is current pfo-------------------------------------------------------------
	    std::cout << "pfoId = " << pPfo->GetParticleId() << std::endl;
	    int pdg = pPfo->GetParticleId();
	    float ori = 0.0;
	    float raw = 0.0;
	    float front = 0.0;
	    int countuse = 0;

	    //  float diff = 0.0;
	    

	    // if (isContainedTrack == true && id == 2212){
	    if (isContainedTrack == true){
	      float foundenergy = 0;
	      float endenergy = 0;
	      float average_dEdx = 0;
	      float total_dEdx = 0;
	      float total_dEdx_raw = 0;

	      if(forwardsChiSquared > backwardsChiSquared) {
		minchisquaredval = backwardsChiSquared;
		//std::cout << "minchisquaredval forwards  " << minchisquaredval << std::endl;


		endenergy = b_Ee;
		//foundenergy = endenergy + andyenergy;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = b_dEdx_2D_average;

		total_dEdx = dEdx_dE_total_b;
		total_dEdx_raw = dEdx_dE_total_raw_b;

		ori = dEdx_2D_av_b;
		raw = dEdx_2D_raw_av_b;

		countuse = backwardscount;

		predenergy = 0.00359282*(std::pow(raw,-1.5)) + 1.04223;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
		//	diff = dEdx_2D_raw_av_b_diff;
	      }
	      else if (forwardsChiSquared <= backwardsChiSquared) {
		minchisquaredval = forwardsChiSquared;
		//	std::cout << "minchisquaredval backwards   " << minchisquaredval << std::endl;

		endenergy = f_Ee;
		//	foundenergy = endenergy + andyenergy;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = f_dEdx_2D_average;

		total_dEdx = dEdx_dE_total_f;
		total_dEdx_raw = dEdx_dE_total_raw_f;

		ori = dEdx_2D_av_f;
		raw = dEdx_2D_raw_av_f;
		//	diff = dEdx_2D_raw_av_f_diff;
		countuse = forwardscount;

		predenergy = 0.00359282*(std::pow(raw,-1.5)) + 1.04223;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}	
	      }
	      std::cout << "   " << std::endl;
	      std::cout << "minchisquaredval  " << minchisquaredval << std::endl;
	      std::cout << "foundenergy " << foundenergy << std::endl;
	      std::cout << "average_dEdx " << average_dEdx  << std::endl;
	      std::cout << "raw dEdx " << raw << std::endl;
	      std::cout << "   " << std::endl;
	      std::cout << "predenergy " << predenergy  << std::endl;

	      // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "foundenergy", andyenergy)); 
	      //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "foundenergyb", t_dEdx_total )); 
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "trueenergy", mcenergyp));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "Id", id));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "IdRecon", pdg));
	      // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "tracklength", f_L));
	      //  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "tracklengthoriginal", trackLength));
	      //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "foundenergyoriginal", andyenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "endenergy", endenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "ori_dEdx", ori ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "raw_dEdx", raw ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "binnedvector", binnedsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "hitvector", hitsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "countbins", countuse));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "total_dEdx", total_dEdx));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "total_dEdx_raw", total_dEdx_raw));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreeb", "predenergy", predenergy));
	      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreeb"));
	      //  PANDORA_MONITORING_API(ScanTree(this->GetPandora(), "ttreeb"));
	    }

	    //  if (isUnContainedTrack == true && id==13 && binnedsize > 10){
	    if (isUnContainedTrack == true){
	      float foundenergy = 0;
	      float endenergy = 0;
	      float average_dEdx = 0;

	      if(forwardsChiSquared > backwardsChiSquared) {
		minchisquaredval = backwardsChiSquared;
		//	std::cout << "minchisquaredval back  " << minchisquaredval << std::endl;
		endenergy = b_Ee;
		//foundenergy = endenergy + andyenergy;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = b_dEdx_2D_average;

		ori = dEdx_2D_av_b;
		raw = dEdx_2D_raw_av_b;
		std::cout << "dEdx_2D_av_b " << raw << std::endl;
		//	diff = dEdx_2D_raw_av_b_diff;
		countuse = backwardscount;
		
		predenergy = 0.000165797*(std::pow(raw,-2.5)) + 0.2;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
		
	      }
	      else if (forwardsChiSquared <= backwardsChiSquared) {
		minchisquaredval = forwardsChiSquared;
		//	std::cout << "minchisquaredval forwards   " << minchisquaredval << std::endl;
		endenergy = f_Ee;
		//	foundenergy = endenergy + andyenergy;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = f_dEdx_2D_average;

		ori = dEdx_2D_av_f;
		raw = dEdx_2D_raw_av_f;

		std::cout << "dEdx_2D_av_b " << raw << std::endl;
		//	diff = dEdx_2D_raw_av_f_diff;
		countuse = forwardscount;

		predenergy = 0.000165797*(std::pow(raw,-2.5)) + 0.2;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
		
	      }
	      std::cout << "   " << std::endl;
	      std::cout << "minchisquaredval  " << minchisquaredval << std::endl;
	      std::cout << "foundenergy " << foundenergy << std::endl;
	      std::cout << "average_dEdx " << average_dEdx  << std::endl;
	      std::cout << "raw dEdx " << raw << std::endl;
	      std::cout << "   " << std::endl;
	      std::cout << "predenergy " << predenergy  << std::endl;
	      
	      // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "foundenergy", foundenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "trueenergy", mcenergyp));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "Id", id));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "IdRecon", pdg));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "endenergy", endenergy));
	      //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "average_dEdx_c", t_dEdx_ave ));
	      //  PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "diff_dEdx", diff ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "ori_dEdx", ori ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "raw_dEdx", raw ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "binnedvector", binnedsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "hitvector", hitsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "countbins", countuse));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "predenergy", predenergy));
	      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreec"));
	      // PANDORA_MONITORING_API(ScanTree(this->GetPandora(), "ttreec"));
	    }


	    if (isContainedShower == true){
	      float foundenergy = 0;
	      float endenergy = 0;
	      float average_dEdx = 0;
	      float total_dEdx = 0;
	      float total_dEdx_raw = 0;

	      if(forwardsChiSquared > backwardsChiSquared) {
		minchisquaredval = backwardsChiSquared;
		//	std::cout << "minchisquaredval forwards  " << minchisquaredval << std::endl;

		endenergy = b_Ee;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = b_dEdx_2D_average;

		total_dEdx = dEdx_dE_total_b;
		total_dEdx_raw = dEdx_dE_total_raw_b;

		ori = dEdx_2D_av_b;
		raw = dEdx_2D_raw_av_b;

		front = front_b_dE;
		countuse = forwardscount;

		//	predenergy = 50.7615*raw + 0.0368726 - 1.842;
		predenergy = (raw - 0.0270035)/0.0137571;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
	      }
	      else if (forwardsChiSquared <= backwardsChiSquared) {
		minchisquaredval = forwardsChiSquared;
		//	std::cout << "minchisquaredval backwards   " << minchisquaredval << std::endl;
		endenergy = f_Ee;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = f_dEdx_2D_average;

		total_dEdx = dEdx_dE_total_f;
		total_dEdx_raw = dEdx_dE_total_raw_f;

		ori = dEdx_2D_av_f;
		raw = dEdx_2D_raw_av_f;

		front = front_f_dE;
		countuse = backwardscount;

		predenergy = (raw - 0.0270035)/0.0137571;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
	      }
	      std::cout << "   " << std::endl;
	      std::cout << "minchisquaredval  " << minchisquaredval << std::endl;
	      std::cout << "foundenergy " << foundenergy << std::endl;
	      std::cout << "average_dEdx " << average_dEdx  << std::endl;
	      std::cout << "raw dEdx " << raw << std::endl;
	      std::cout << "   " << std::endl;
	      std::cout << "predenergy " << predenergy  << std::endl;

	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "foundenergy", andyenergy)); 
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "foundenergyb", t_dEdx_total )); 
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "trueenergy", mcenergyp));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "Id", id));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "IdRecon", pdg));
	      //   PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "foundenergyoriginal", andyenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "endenergy", endenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "ori_dEdx", ori ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "raw_dEdx", raw ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "binnedvector", binnedsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "hitvector", hitsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "countbins", countuse));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "total_dEdx", total_dEdx));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "total_dEdx_raw", total_dEdx_raw));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "front", front));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "predenergy", predenergy));
	      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreed"));
	      //  PANDORA_MONITORING_API(ScanTree(this->GetPandora(), "ttreed"));
	    }

	    //  if (isUnContainedTrack == true && id==13 && binnedsize > 10){
	    if (isUnContainedShower == true){
	      float foundenergy = 0;
	      float endenergy = 0;
	      float average_dEdx = 0;

	      if(forwardsChiSquared > backwardsChiSquared) {
		minchisquaredval = backwardsChiSquared;
		//	std::cout << "minchisquaredval back  " << minchisquaredval << std::endl;
		endenergy = b_Ee;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = b_dEdx_2D_average;

		ori = dEdx_2D_av_b;
		raw = dEdx_2D_raw_av_b;

		front = front_b_dE;
		countuse = forwardscount;

		//	predenergy = 104.239*raw + 0.0649419 - 3.064;
		predenergy = (raw + 0.058042)/0.0780376;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
	      }
	      else if (forwardsChiSquared <= backwardsChiSquared) {
		minchisquaredval = forwardsChiSquared;
		//	std::cout << "minchisquaredval forwards   " << minchisquaredval << std::endl;
		endenergy = f_Ee;
		foundenergy = t_dEdx_total + endenergy;
		average_dEdx = f_dEdx_2D_average;

		ori = dEdx_2D_av_f;
		raw = dEdx_2D_raw_av_f;

		front = front_f_dE;
		countuse = backwardscount;

		predenergy = (raw + 0.058042)/0.0780376;
		if (predenergy < 0.0) {
		  predenergy = 0.0;
		}
	      }
	      std::cout << "   " << std::endl;
	      std::cout << "minchisquaredval  " << minchisquaredval << std::endl;
	      std::cout << "foundenergy " << foundenergy << std::endl;
	      std::cout << "average_dEdx " << average_dEdx  << std::endl;
	      std::cout << "raw dEdx " << raw << std::endl;
	      std::cout << "   " << std::endl;
	      std::cout << "predenergy " << predenergy  << std::endl;
	      
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "foundenergy", foundenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "trueenergy", mcenergyp));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "Id", id));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "IdRecon", pdg));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "endenergy", endenergy));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "ori_dEdx", ori ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "raw_dEdx", raw ));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "binnedvector", binnedsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "hitvector", hitsize));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "countbins", countuse));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "front", front));
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreee", "predenergy", predenergy));
	      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreee"));
	      //  PANDORA_MONITORING_API(ScanTree(this->GetPandora(), "ttreee"));
	    }


	  }
	}
	else {
	  std::cout << "wrong cluster list size" << std::endl;
	}
	if (predenergy > 0.0) {
	  prednuenergy = prednuenergy + predenergy;	  
	}
      }
  

    std::cout << "   " << std::endl;
    std::cout << " ------------------  " << std::endl;
    std::cout << "prednuenergy: " << prednuenergy << " GeV " << std::endl;
    std::cout << " ------------------  " << std::endl;
    

    //-----------------------------------------------------------------------------------------------------------------------------------
     int Size = 0;
    std::list<std::pair<const pandora::MCParticle*, int>> mchitslist;
    std::list<std::pair<const pandora::MCParticle*, int>> mcpairslist;
    MCParticleList mcparticlelist;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it2;
    CaloHitList totalcalohits;
    // double mcenergy = 0;

    float truenuenergy = 0.0;

    for (const Pfo *const pPPPPfo : *pPfoList ) {
      try{
	const MCParticle *pMCParticle = LArMCParticleHelper::GetMainMCParticle(pPPPPfo);
	mcparticlelist.push_back(pMCParticle);
	const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
	std::cout << " -----------------------" << std::endl;
	std::cout << "TheMC Parent: " << std::endl;
	std::cout << "pParentMCParticle energy  " << pParentMCParticle->GetEnergy() << std::endl;
	std::cout << "pParentMCParticle Id  " << pParentMCParticle->GetParticleId() << std::endl;
	if (pParentMCParticle->GetParticleId() == 14) {
	  truenuenergy = pParentMCParticle->GetEnergy();
	}
	//	mcenergy = pParentMCParticle->GetEnergy();
	std::cout << "is mc parent neutrino " << LArMCParticleHelper::IsNeutrino(pParentMCParticle) << std::endl;
	std::cout << "   " << std::endl;
	std::cout << "TheMC Particle: " << std::endl;
	std::cout << pMCParticle << " <----- MC particle " << std::endl;
	std::cout << pMCParticle->GetEnergy() << " <----- MC particle energy " << std::endl;
	std::cout << "mc nuance code of particle " <<LArMCParticleHelper::GetNuanceCode(pMCParticle) << std::endl;
	std::cout << " -----------------------" << std::endl;
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_U, totalcalohits);
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_V, totalcalohits);
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_W, totalcalohits);
      } catch (...){
	std::cout << "lost one of the mc particles" << std::endl;
      }

    }
     
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreef", "truenunergy", truenuenergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreef", "prednuenergy", prednuenergy));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreef"));
  
    // std::cout << "mc energy before tree " << mcenergy << std::endl;

     //get the MC particles to hits map
    LArMCParticleHelper::MCContributionMap mcMCParticlesToGoodHitsMap;
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 0;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_minPrimaryGoodViews = 0;
    parameters.m_minHitSharingFraction = 0;
    parameters.m_selectInputHits = true;

    LArMCParticleHelper::SelectReconstructableMCParticles(&mcparticlelist, &totalcalohits, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcMCParticlesToGoodHitsMap); 
    //LArMCParticleHelper::SelectReconstructableMCParticles(&mcparticlelist, &totalcalohits, parameters, 1, mcMCParticlesToGoodHitsMap); 

  
    //get the pfos to hits map
    LArMCParticleHelper::PfoContributionMap PfosToGoodHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, {mcMCParticlesToGoodHitsMap}, PfosToGoodHitsMap, false);  //changed


    //find hits shared between Pfos and MCParticles
    LArMCParticleHelper::PfoToMCParticleHitSharingMap PfotoMCParticleMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap ParticletoPfoMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(PfosToGoodHitsMap, {mcMCParticlesToGoodHitsMap}, PfotoMCParticleMap, ParticletoPfoMap);
    
    MCParticleVector mcParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcMCParticlesToGoodHitsMap}, mcParticleVector);
    PfoVector pfoVector;
    LArMonitoringHelper::GetOrderedPfoVector({PfosToGoodHitsMap}, pfoVector);



    unsigned int nMatches = std::numeric_limits<unsigned int>::max();
    LArMonitoringHelper::PrintMatchingTable(pfoVector, mcParticleVector, ParticletoPfoMap, nMatches);


    for (const auto &pMCParticle : mcParticleVector)
      {
	const auto &caloHitList2 = mcMCParticlesToGoodHitsMap.at(pMCParticle);
	std::cout << "  " << std::endl;
	std::cout << "Primary MCParticle " << pMCParticle << std::endl;
	std::cout << "  - PDG : " << pMCParticle->GetParticleId() << std::endl;
	std::cout << "  - NHits : " << caloHitList2.size() << std::endl; 
	std::cout << "  - Energy : " << pMCParticle->GetEnergy() << std::endl; 
	Size = caloHitList2.size();
	mchitslist.push_back(std::make_pair(pMCParticle, Size));
	std::cout << "  " << std::endl;
      
      }


     return STATUS_CODE_SUCCESS;
  }

  StatusCode NDQualitiesAlg::ReadSettings(const TiXmlHandle xmlHandle)
  {
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
													 "pfoListName", m_pfoListName));

    
    return STATUS_CODE_SUCCESS;
  }
  //------------------------------------------------------------------------------------------------------------------
  void NDQualitiesAlg::FillHitChargeVector(const Cluster *const pCluster, HitChargeVector &hitChargeVector,  CartesianVector &endPoint1, CartesianVector &endPoint2 )
  {
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFit(this->GetCachedSlidingFit(pCluster));

    endPoint1 = slidingFit.GetGlobalMinLayerPosition();
    endPoint2 = slidingFit.GetGlobalMaxLayerPosition();
     
    std::cout <<" end points (Fill Hit Charge Vector) :  " << endPoint1 <<"    " << endPoint2 << std::endl;


    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
      {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float hitWidth(pCaloHit->GetCellSize1());

        float caloHitEnergy(pCaloHit->GetInputEnergy());
        caloHitEnergy *= 273.5; //ADC to electron
        caloHitEnergy *= 23.6/1000000; //ionisation energy per electron in MeV
        caloHitEnergy /= 0.62;

        float rL(0.f), rT(0.f);
        slidingFit.GetLocalPosition(caloHitPosition, rL, rT);
        if (rL == 0.)
	  continue;

        float calibratedUncertainty(std::sqrt((0.00419133 * (caloHitEnergy/hitWidth) * (caloHitEnergy/hitWidth)) + (0.00967141 * (caloHitEnergy/hitWidth)))); //70%
        HitCharge hitCharge(pCaloHit, rL, hitWidth, caloHitEnergy, calibratedUncertainty);
        hitChargeVector.push_back(hitCharge);
      }

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);

  }

  //------------------------------------------------------------------------------------------------------------
  
  void NDQualitiesAlg::GetTrackLength(HitChargeVector &hitChargeVector, float &trackLength)
  {
    trackLength = 0.f;

    for (HitCharge &hitCharge : hitChargeVector)
      {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
	  trackLength = hitCharge.GetLongitudinalPosition();
      }
  }
  //----------------------------------------------------------------------------
  void NDQualitiesAlg::BinHitChargeVector(lar_content::NDQualitiesAlg::HitChargeVector &hitChargeVector, lar_content::NDQualitiesAlg::HitChargeVector &binnedHitChargeVector)
  {
    //This is always commented out in TrackDirectionTool.cc  // Jan 2020: may need it!
    float binSize = (hitChargeVector.size() > 50 ? (0.5 + (hitChargeVector.size() - 50) * 2.5/300) : 0.0);

    if (binSize <= 0.0)
      {
        binnedHitChargeVector = hitChargeVector;
        return;
      }
    else if (binSize > 4.0)
      binSize = 4.0;


    float trackLength(0.f);

    for (lar_content::NDQualitiesAlg::HitCharge &hitCharge : hitChargeVector)
      {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
	  trackLength = hitCharge.GetLongitudinalPosition();
      }

    for (float i = binSize; i <= trackLength; i += binSize)
      {
        int nHitsBin(0);
        float meanBinPosition(0.f), meanBinWidth(0.f);
        float meanBinCharge(0.f), sumSquaredSigmas(0.f);

        for (lar_content::NDQualitiesAlg::HitCharge &hitCharge : hitChargeVector)
	  {
            if (hitCharge.GetLongitudinalPosition() > i)
	      break;

            if (hitCharge.GetLongitudinalPosition() < i && hitCharge.GetLongitudinalPosition() >= (i - binSize))
	      {
                nHitsBin++;
                meanBinPosition += hitCharge.GetLongitudinalPosition();
                meanBinWidth += hitCharge.GetHitWidth();

                meanBinCharge += hitCharge.GetCharge();
                sumSquaredSigmas += (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
	      }
	  }

        if (nHitsBin == 0)
	  continue;


        meanBinPosition /= nHitsBin;
        meanBinWidth /= nHitsBin;
        meanBinCharge /= nHitsBin;
        sumSquaredSigmas /= (nHitsBin * nHitsBin);

        float binUncertainty(std::sqrt(sumSquaredSigmas));
        lar_content::NDQualitiesAlg::HitCharge binnedHitCharge(NULL, meanBinPosition, meanBinWidth, meanBinCharge, binUncertainty);

	if (binnedHitCharge.GetCharge() > 0.1)                  // why?
	  binnedHitChargeVector.push_back(binnedHitCharge);
      }
  }
  
  //----------------------------------------------------------------------------------------------------------

  double NDQualitiesAlg::GetEnergyfromLength(lar_content::NDQualitiesAlg::LookupTable &lookupTable, double &trackLength)
  {
    std::map<int, double> lookupMap(lookupTable.GetMap());
    double binWidth(lookupTable.GetBinWidth());

    if (trackLength >= lookupTable.GetMaxRange())
      return 0.5; //0 energy means infinite dE/dx
    else if (trackLength <= 1.0)
      return lookupTable.GetInitialEnergy();

    int n(std::floor(trackLength/binWidth));
    //double sampleLength = lookupTable.GetMaxRange() - trackLength;
    //int n(std::floor((sampleLength)/binWidth));


    // std::cout << "n: " << n << std::endl;
    //  for (const auto& entry : lookupMap) {
    //   std::cout << entry.first << "  "  << entry.second << std::endl;
    // }
    std::map<int, double>::iterator nextEntryIterator(lookupMap.upper_bound(n)); //points to n+1
    std::map<int, double>::iterator previousEntryIterator(std::prev(nextEntryIterator , 1));
    //std::cout << std::distance(lookupMap.begin(),nextEntryIterator) << std::endl;


    double leftLength(previousEntryIterator->first * binWidth), rightLength(nextEntryIterator->first * binWidth);
    double leftEnergy(previousEntryIterator->second), rightEnergy(nextEntryIterator->second);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(leftEnergy - rightEnergy);

    //double finalEnergy(leftEnergy - (((sampleLength - leftLength)/(lengthDifference)) * (energyDifference)));
    double finalEnergy(leftEnergy - (((rightLength - leftLength)/(lengthDifference)) * (energyDifference)));

    //very small energy values leadt to huge dE/dx values: truncate
    if (finalEnergy <= 2.0)
      return 2.0;
    else
      return finalEnergy;
  }
  //----------------------------------------------------------------------------
  double NDQualitiesAlg::GetLengthfromEnergy(lar_content::NDQualitiesAlg::LookupTable &lookupTable, double &currentEnergy)
  {
    std::map<double, int> reverseLookupMap(lookupTable.GetReverseMap());
    double binWidth(lookupTable.GetBinWidth());

    if (currentEnergy <= 0.0)
      return lookupTable.GetMaxRange();
    else if (currentEnergy >= lookupTable.GetInitialEnergy())
      return 0.0;

    std::map<double, int>::iterator nextEntryIterator(reverseLookupMap.upper_bound(currentEnergy));
    std::map<double, int>::iterator previousEntryIterator(std::prev(nextEntryIterator , 1));

    double upperEnergy(nextEntryIterator->first), lowerEnergy(previousEntryIterator->first);
    double leftLength(nextEntryIterator->second * binWidth), rightLength(previousEntryIterator->second * binWidth);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(upperEnergy - lowerEnergy);

    return leftLength + (((upperEnergy - currentEnergy)/energyDifference) * lengthDifference);
  }
  
  
  //----------------------------------------------------------------------------
  const TwoDSlidingFitResult &NDQualitiesAlg::GetCachedSlidingFit(const Cluster *const pCluster) const
  {
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
      {
        std::cout << "Sliding fit retrieval failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
      }

    return iter->second;
  }
  //------------------------------------------------------------------------
  bool NDQualitiesAlg::SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2)
  {
    return hitCharge1.GetLongitudinalPosition() < hitCharge2.GetLongitudinalPosition();
  }

  //----------------------------------------------------------------------------------------------------------------------------------
  void NDQualitiesAlg::AddToSlidingFitCache(const Cluster *const pCluster)
  {

    TwoDSlidingFitResultMap slidingFitResultMap;
    if (slidingFitResultMap.size() != 0) {
      std::unordered_map<const pandora::Cluster*,lar_content::TwoDSlidingFitResult>::const_iterator got = slidingFitResultMap.find(pCluster);
      if ( got !=  slidingFitResultMap.end() ) {
	return;
      }
    }
 
    if (slidingFitResultMap.find(pCluster) != slidingFitResultMap.end()){
      return;
    }
  
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, 20, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFit)).second)
      {
        std::cout << "Sliding fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
      }
  }

  //----------------------------------------------------------------------------------------------------------------------------------

  double NDQualitiesAlg::DensityCorrection(double &T, double &M)
  {
    const double C = 5.2146;
    const double a = 0.19559;
    const double m = 3.0;
    const double X1 = 3.0;
    const double X0 = 0.2000;
    const double delta0 = 0.0;

    double p = std::sqrt((T*T) + 2*T*M);
    double gamma = std::sqrt(1 + ((p/M) * (p/M)));
    double beta = std::sqrt(1 - 1 / (gamma*gamma));
    double X = std::log10(beta*gamma);

    if (X < X0)
      return delta0;
    else if ((X > X0) && (X < X1))
      return 2 * X * std::log(10) - C + (a * (std::pow((X1 - X), m)));
    else
      return 2 * X * std::log(10) + C;
  }

  
  double NDQualitiesAlg::BetheBloch(double &T, double &M)
  {
    const double K = 0.307075; // constant K in MeV cm mol^-1
    const double z = 1; // charge in e
    const double Z = 18; // Atomic number Z
    const double A = 39.948; // Atomic mass in g mol-1
    //double M = 105.7; // Mass of heavy particle in MeV
    const double m_e = 0.511; // Mass of electron in MeV
    const double rho = 1.396; // Density of material in g cm^-3 (here: argon density)
    const double I = 0.000188; // Ionisation energy in MeV


    double p(std::sqrt((T*T) + 2*T*M));
    double gamma(std::sqrt(1 + ((p/M) * (p/M))));
    double beta(std::sqrt(1 - 1 / (gamma*gamma)));

    double T_max(2 * m_e * (beta*gamma*beta*gamma) / (1 + 2 * gamma * m_e / M + ((m_e/M) * (m_e/M))));

    return rho * ((K * z * z * Z) / A) * (0.5*std::log(2 * m_e * T_max * (beta*gamma*beta*gamma) / (I*I) ) - (beta*beta) - (0.5*DensityCorrection(p, M))) / (beta*beta); //in MeV/cm

  }

   //----------------------------------------------------------------------------------------------------------------------------------

  void NDQualitiesAlg::FillLookupTable(lar_content::NDQualitiesAlg::LookupTable& lookupTable, double M)
  {
    std::map<int, double> lookupMap;
    std::map<double, int> reverseLookupMap;

    double currentEnergy(lookupTable.GetInitialEnergy()), binWidth(lookupTable.GetBinWidth());
    int maxBin(0);

    std::cout << "binwidth FLT " << binWidth << std::endl;
    

    for (double n = 0; n < 100000; ++n)  //lost a 0
// for (double n = 0; n < 10000; ++n)
      {
        double currentdEdx = BetheBloch(currentEnergy, M);

	//	std::cout << "current E " << currentEnergy << std::endl;
	//	std::cout << "currentdEdx " << currentdEdx << std::endl;
	//	std::cout << " ---------------------------- " << std::endl;

        if ((currentdEdx * binWidth) >= currentEnergy)
	  {
            double maxRange = (n * binWidth) + (currentEnergy/currentdEdx);
            lookupTable.SetMaxRange(maxRange);
            maxBin = n;

            lookupMap.insert(std::pair<int, double>(n, 0.0));
            reverseLookupMap.insert(std::pair<double, int>(0.0, n));
            break;
	  }
        else
	  {
            lookupMap.insert(std::pair<int, double>(n, currentEnergy));
            reverseLookupMap.insert(std::pair<double, int>(currentEnergy, n));
	  }
	///std::cout << "current energy: " << currentEnergy << std::endl;
	//std::cout << "n: " << n << std::endl;
	//	std::cout << "currentdEdx: " << currentdEdx << std::endl;
	
        currentEnergy -= (currentdEdx * binWidth);
      }

    // std::cout << "maxBin " << maxBin << std::endl;

    //double maxRange(lookupTable.GetMaxRange());

    //remove redundant entries to make lookup much faster
    
    for (std::map<int, double>::iterator it = lookupMap.begin(); it != lookupMap.end(); it++)
      {
        double n(it->first);
        double val(it->second);
        double distanceFromMaxRange((maxBin - n) * binWidth);

        if (n == 0 || n == maxBin)
	  continue;
        else
	  {
            double distanceMagnitude(floor(distanceFromMaxRange/2.0));
            double samplingDistance((1 + distanceMagnitude) * binWidth);


            if (!(remainder(distanceFromMaxRange, samplingDistance) == 0.0))
	      {
                lookupMap.erase(n);
                reverseLookupMap.erase(val);
	      }
	  }
      }
    

    //std::cout << "max range = " << lookupTable.GetMaxRange() << "    bin " << maxBin << std::endl;
    lookupTable.SetMap(lookupMap);
    lookupTable.SetReverseMap(reverseLookupMap);

    if (lookupTable.GetMaxRange() == 0.f)
      std::cout << "Warning: the lookup table max range has not been correctly set." << std::endl;
  }

// }
}
