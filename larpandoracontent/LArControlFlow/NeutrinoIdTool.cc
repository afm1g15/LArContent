/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/NeutrinoIdTool.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

  //typedef std::unordered_map<const pandora::ParticleFlowObject*, float> PfoToFloatMap;

NeutrinoIdTool::NeutrinoIdTool() :
    m_useTrainingMode(false),
    m_selectNuanceCode(false),
    m_nuance(-std::numeric_limits<int>::max()),
    m_minPurity(0.9f),
    m_minCompleteness(0.9f),
    m_minProbability(0.0f),
    m_maxNeutrinos(1),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

 /*
NeutrinoIdTool::~NeutrinoIdTool() {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree6", "output6.root", "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << " Unable to write tree  to file " << std::endl;
        }
  }
  */

//------------------------------------------------------------------------------------------------------------------------------------------

  void NeutrinoIdTool::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos, const PfoToFloatMap &pfotoprobabilitymapb, const SliceVector &sliceVector)
{
  std::cout << "slice size neutrino Id " << sliceVector.size() << std::endl;
  std::cout << "TRAINING MODE:  " << m_useTrainingMode << std::endl;
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nSlices(nuSliceHypotheses.size());
    if (nSlices == 0) return;

    SliceFeaturesVector sliceFeaturesVector;
    this->GetSliceFeatures(this, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, pfotoprobabilitymapb);

    if (m_useTrainingMode)
    {
      std::cout << "Training Mode!" << std::endl;
        // ATTN in training mode, just return everything as a cosmic-ray
        this->SelectAllPfos(pAlgorithm, crSliceHypotheses, selectedPfos);

	std::cout << "Getting the bestSliceIndex!" << std::endl;
        unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
        if (!this->GetBestMCSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex)) return;

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
	  std::cout << "Looking at the SliceFeatures!" << std::endl;
            const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
            if (!features.IsFeatureVectorAvailable()) continue;

            LArMvaHelper::MvaFeatureVector featureVector;
	    std::cout << "Getting the featureVector!" << std::endl;
            features.GetFeatureVector(featureVector);
	    std::cout << "Producing the TrainingExample!" << std::endl;
            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, sliceIndex == bestSliceIndex, featureVector);
        }

        return;
    }

    this->SelectPfosByProbability(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos, pfotoprobabilitymapb);

}

//------------------------------------------------------------------------------------------------------------------------------------------

  void NeutrinoIdTool::GetSliceFeatures(const NeutrinoIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector, const PfoToFloatMap &pfotoprobabilitymapb) const
{
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
      sliceFeaturesVector.push_back(SliceFeatures(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), pTool, pfotoprobabilitymapb));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetBestMCSliceIndex(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex) const
{
    unsigned int nHitsInBestSlice(0), nNuHitsInBestSlice(0);
    std::cout << "Inside the bestsliceindex" << std::endl;

    // Get all hits in all slices to find true number of mc hits
    const CaloHitList *pAllReconstructedCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pAllReconstructedCaloHitList));
    std::cout << "NeutrinoTool: CaloHitList: " << pAllReconstructedCaloHitList->size() << std::endl;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));
    std::cout << "NeutrinoTool: MCHitList: " << pMCParticleList->size() << std::endl;

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList reconstructableCaloHitList;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectCaloHits(pAllReconstructedCaloHitList, mcToPrimaryMCMap, reconstructableCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);
    std::cout << "pAllReconstructedCaloHitList : " << pAllReconstructedCaloHitList->size() << std::endl;

    const int nuNHitsTotal(this->CountNeutrinoInducedHits(reconstructableCaloHitList));
    std::cout << "nuNHitsTotal " << nuNHitsTotal << std::endl;

    const CaloHitSet reconstructableCaloHitSet(reconstructableCaloHitList.begin(), reconstructableCaloHitList.end());
    std::cout << "reconstructableCaloHitSet size " << reconstructableCaloHitSet.size() << std::endl;


    CaloHitList parentCaloHitList;
    for (const CaloHit  *const pCaloHit : reconstructableCaloHitList)
      {
	//	std::cout << "compare calo: " << pCaloHit << std::endl;
	const CaloHit *const pParentHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
	parentCaloHitList.push_back(pParentHit);
	//	std::cout << "compare parent: " << pParentHit << std::endl;
      }

    std::cout << "looping through" << std::endl;
    std::cout << "-------------------------------T-----------------" << std::endl;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
      std::cout << "//////////////////New loop\\\\\\\\\\\\\\\\\\\\ " << std::endl;

        CaloHitList reconstructedCaloHitList;
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedCaloHitList, reconstructableCaloHitSet, reconstructableCaloHitList);
	
        for (const ParticleFlowObject *const pNeutrino : nuSliceHypotheses.at(sliceIndex))
        {
            const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());
            this->Collect2DHits(nuFinalStates, reconstructedCaloHitList, reconstructableCaloHitSet, reconstructableCaloHitList);
	    std::cout << "nuFinalStates.size() " << nuFinalStates.size() << std::endl;
        }
	
	std::cout <<"reconstructedCaloHitList.size() " << reconstructedCaloHitList.size() << std::endl;
        const unsigned int nNuHits(this->CountNeutrinoInducedHits(reconstructedCaloHitList));
	std::cout << "nuHits " << nNuHits << std::endl;
	std::cout << "nNuHitsInBestSlice " << nNuHitsInBestSlice << std::endl;

        if ((nNuHits > nNuHitsInBestSlice) == true)
        {
	  nNuHitsInBestSlice = nNuHits;
	  nHitsInBestSlice = reconstructedCaloHitList.size();
	  bestSliceIndex = sliceIndex;
        }
    }
    std::cout << "-------------------------------B-----------------" << std::endl;

    // ATTN for events with no neutrino induced hits, default neutrino purity and completeness to zero
    const float purity(nHitsInBestSlice > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nHitsInBestSlice) : 0.f);
    const float completeness(nuNHitsTotal > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nuNHitsTotal) : 0.f);

    std::cout << "------" << std::endl;
    std::cout << "nHitsInBestSlice " << nHitsInBestSlice << std::endl;
    std::cout << "nuNHitsTotal " << nuNHitsTotal << std::endl;
    std::cout << "------" << std::endl;
    return this->PassesQualityCuts(pAlgorithm, purity, completeness);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::PassesQualityCuts(const Algorithm *const pAlgorithm, const float purity, const float completeness) const
{
  std::cout << "1st cut " << (purity < m_minPurity || completeness < m_minCompleteness) << std::endl;
  std::cout << "purity: " << purity << "  m_minPurity " << m_minPurity << std::endl;
  std::cout << "completeness: " << completeness << " m_minCompleteness  " << m_minCompleteness  << std::endl;
  std::cout << "2nd cut " << (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance)) << std::endl;

    if (purity < m_minPurity || completeness < m_minCompleteness) return false;
    if (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance)) return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

  void NeutrinoIdTool::Collect2DHits(const PfoList &pfos, CaloHitList &reconstructedCaloHitList, const CaloHitSet &reconstructableCaloHitSet,  CaloHitList &matchCaloHitList) const
{
  if (1==2) {
    std::cout << reconstructableCaloHitSet.size() << std::endl;
  }
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    std::cout << "collectedHits.size() " << collectedHits.size() << std::endl;

    for (const CaloHit *const pCaloHit : collectedHits)
      { //see if it's in container
	// if (!reconstructableCaloHitSet.count(pCaloHit)) {
	//	std::cout << "continued." << std::endl;
	//	continue;
	//  }
	const CaloHit *const pParentHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
	if (std::find(matchCaloHitList.begin(), matchCaloHitList.end(), pParentHit) == matchCaloHitList.end()){
	  continue;
	}


        // Ensure no hits have been double counted
        if (std::find(reconstructedCaloHitList.begin(), reconstructedCaloHitList.end(), pCaloHit) == reconstructedCaloHitList.end())
            reconstructedCaloHitList.push_back(pParentHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::CountNeutrinoInducedHits(const CaloHitList &caloHitList) const
{
    unsigned int  nNuHits(0);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
	  if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit)))){ 
	    nNuHits =  nNuHits + 1;
	  }
        }
        catch (const StatusCodeException &)
        {
        }
    }

    std::cout << "counting...nuHits: " << nNuHits << std::endl;
    return nNuHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int NeutrinoIdTool::GetNuanceCode(const Algorithm *const pAlgorithm) const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    // std::cout << "True neutrino : " << trueNeutrinos.front() << std::endl;

    if (trueNeutrinos.size() != 1)
    {
        std::cout << "NeutrinoIdTool::GetNuanceCode - Error: number of true neutrinos in event must be exactly one" << std::endl;
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    return LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectAllPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &hypotheses, PfoList &selectedPfos) const
{
    for (const PfoList &pfos : hypotheses)
    {
        for (const ParticleFlowObject *const pPfo : pfos)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        this->SelectPfos(pfos, selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

  void NeutrinoIdTool::SelectPfosByProbability(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, PfoList &selectedPfos, const PfoToFloatMap &pfotoprobabilitymapb) const
  {
    if (1 == 2) {
      std::cout <<pfotoprobabilitymapb.size() << std::endl; 
    }

    std::cout << nuSliceHypotheses.size() << " <----- nuSliceHypotheses.size() (nSlices)" << std::endl;
    float maxprobcr_f = -1;
    float maxprobnu_f = -1;
    float minprobcr_f = -1;
    float minprobnu_f = -1;
    std::cout << maxprobcr_f << maxprobnu_f << std::endl;
    std::cout << minprobcr_f << minprobnu_f << std::endl;
    // Calculate the probability of each slice that passes the minimum probability cut
    std::vector<UintFloatPair> sliceIndexProbabilityPairs;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
      {
	std::cout << "-----------------start of new slice loop---------------------------------" << std::endl;
	std::vector<float> downprobcr;
	std::vector<float> downprobnu;
        //const float nuProbability(sliceFeaturesVector.at(sliceIndex).GetNeutrinoProbability(m_supportVectorMachine));
	//std::cout << "nuProbability : " <<  nuProbability << std::endl;

        for (const ParticleFlowObject *const pPfo : crSliceHypotheses.at(sliceIndex))
	  {
	    object_creation::ParticleFlowObject::Metadata metadata;
	    std::cout << "Cosmic: Pfo : " << pPfo << std::endl;
	    // std::cout << "Cosmic: Pfo Id : " << pPfo->GetParticleId() << std::endl;
	    // std::cout << "   " << std::endl;
	    
	    
	    auto search = pfotoprobabilitymapb.find(pPfo);
	    if (search !=  pfotoprobabilitymapb.end()) {
	      std::cout << "Found " << search->first << " " << search->second << '\n';
	      metadata.m_propertiesToAdd["downProb"] = search->second;
	      if(search->second != -1 && search->second != -2  && search->second != -3 && search->second != -4) {
		downprobcr.push_back(search->second);
	      }
	    }
	    
	    //  metadata.m_propertiesToAdd["NuScore"] = nuProbability;
	    //  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
	  }

        for (const ParticleFlowObject *const pPfo : nuSliceHypotheses.at(sliceIndex))
	  {
	    object_creation::ParticleFlowObject::Metadata metadata;
	    std::cout << "Neutrino: Pfo : " << pPfo << std::endl;
	    std::cout << "Neutrino: Pfo Id : " << pPfo->GetParticleId() << std::endl;
	    CaloHitList collectedHits;
	    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, collectedHits);
	    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, collectedHits);
	    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, collectedHits);
	    std::cout << "Neutrino: Hits : " << collectedHits.size() << std::endl;
	    //  std::cout << "   " << std::endl;

	    PfoList daughterPfos = pPfo->GetDaughterPfoList();
	    
	    for (const ParticleFlowObject *const pPPfo : daughterPfos) {
	      std::cout << "Daughter: Pfo : " << pPPfo << std::endl;
	      std::cout << "Daughter: Pfo Id : " << pPPfo->GetParticleId() << std::endl;
	      CaloHitList collectedHitsD;
	      LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, collectedHitsD);
	      LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, collectedHitsD);
	      LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, collectedHitsD);
	      std::cout << "Daughter: Hits : " << collectedHitsD.size() << std::endl;

	      auto search = pfotoprobabilitymapb.find(pPPfo);
	      if (search !=  pfotoprobabilitymapb.end()) {
		std::cout << "Found " << search->first << " " << search->second << '\n';
		metadata.m_propertiesToAdd["downProb"] = search->second;
		if(search->second != -1 && search->second != -2 && search->second != -3 && search->second != -4) {
		  downprobnu.push_back(search->second);
		}
	      }
	    }

	    //   metadata.m_propertiesToAdd["NuScore"] = nuProbability;
	    //  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
	  }

	//could add another cut - alex
	//track needs to be long. cut

	std::cout << "downprobcr size = " << downprobcr.size() << std::endl;
	std::cout << "downprobnu size = " << downprobnu.size() << std::endl;

	if(downprobcr.size() > 0) {
	  std::vector<float>::iterator minprobcr = std::min_element(downprobcr.begin(), downprobcr.end());
	  auto index = std::distance(downprobcr.begin(), minprobcr);   
	  int indexvalue = index;
	  std::cout << " --- " << std::endl;
	  std::cout <<" minprobcr "   << downprobcr[indexvalue] << std::endl;
	  minprobcr_f = downprobcr[indexvalue];
	  std::vector<float>::iterator maxprobcr = std::max_element(downprobcr.begin(), downprobcr.end());
	  auto index2 = std::distance(downprobcr.begin(), maxprobcr);   
	  int indexvalue2 = index2;
	  maxprobcr_f = downprobcr[indexvalue2];
	  std::cout <<" maxprobcr "   << downprobcr[indexvalue2] << std::endl;
	  std::cout << " --- " << std::endl;
	}

	if(downprobnu.size() > 0) {
	  std::vector<float>::iterator minprobnu = std::min_element(downprobnu.begin(), downprobnu.end());
	  auto indexn = std::distance(downprobnu.begin(), minprobnu);   
	  int indexvaluen = indexn;
	  std::cout <<" minprobnu "   << downprobnu[indexvaluen] << std::endl;
	  minprobnu_f = downprobnu[indexvaluen];
	  std::vector<float>::iterator maxprobnu = std::max_element(downprobnu.begin(), downprobnu.end());
	  auto index2n = std::distance(downprobnu.begin(), maxprobnu);   
	  int indexvalue2n = index2n;
	  maxprobnu_f = downprobnu[indexvalue2n];
	  std::cout <<" maxprobnu "   << downprobnu[indexvalue2n] << std::endl;
	  std::cout << " ------- " << std::endl;
	}
	
	//------------------------------------------------------------
	std::cout <<"probs " << maxprobcr_f << maxprobnu_f << std::endl;
	std::cout << minprobcr_f << minprobnu_f << std::endl;

	std::cout << "-----------------nuProbCalc---------------------------------" << std::endl;
	//	std::cout << "m_supportVectorMachine " << m_supportVectorMachine << std::endl;
	std::cout << "sliceindex " <<  sliceIndex << std::endl;
	std::cout << "SFV size " << sliceFeaturesVector.size() << std::endl;
	const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
	if (!features.IsFeatureVectorAvailable()){
	  std::cout << "not available...." << std::endl;
	}
	const float nuProbability(sliceFeaturesVector.at(sliceIndex).GetNeutrinoProbability(m_supportVectorMachine));
	std::cout << "nuProbability : " <<  nuProbability << std::endl;
	std::cout << "SVM Ran Correctly!" << std::endl;

	for (const ParticleFlowObject *const pPfo : crSliceHypotheses.at(sliceIndex))
	  {
	    object_creation::ParticleFlowObject::Metadata metadata;
	    
            metadata.m_propertiesToAdd["NuScore"] = nuProbability;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
	  }

        for (const ParticleFlowObject *const pPfo : nuSliceHypotheses.at(sliceIndex))
	  {
            object_creation::ParticleFlowObject::Metadata metadata;
	      metadata.m_propertiesToAdd["NuScore"] = nuProbability;
	      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
	  }
	      

	//--------------------------------------------------------------
	if (nuProbability < m_minProbability) //if ((nuProbability < 0.7) && ( > 0.96) && (tracklength > 50))  //minprobability is 0 atm...
	  //	if (downprobnu.size() != 0 && downprobcr.size() !=0) {
	  //  if (((nuProbability < 0.5) && ( minprobcr_f > 0.75)) || ((nuProbability < 0.5) && ( minprobnu_f > 0.75)) )
	  {
	    //if below, it's a cosmic ray
	    std::cout << " ^^ taken out as a cosmic ray " << std::endl;
	    std::cout << "------------------------------" << std::endl;
	    this->SelectPfos(crSliceHypotheses.at(sliceIndex), selectedPfos);
	    continue;
	    //  }
	  }

	//	if (minprobnu_f == 0.5 || minprobnu_f == -1) {
	//	if (downprobnu.size() != 0 || downprobcr.size() !=0) {
	// if (((nuProbability < 0.2) && ( minprobcr_f > 0.90)) || ((nuProbability < 0.2) && ( minprobnu_f > 0.90)) ) {
	// sliceIndexProbabilityPairs.push_back(UintFloatPair(sliceIndex, 0.0));  ///multiply all together? add?
	// std::cout << "Edited by prob - 0.2, 0.90" << std::endl;
	    //  }
	    // else {
	    sliceIndexProbabilityPairs.push_back(UintFloatPair(sliceIndex, nuProbability));
	    //  }
	    //	}
	    //	else {
	    //	  sliceIndexProbabilityPairs.push_back(UintFloatPair(sliceIndex, nuProbability));
	    //	}
	    //  }
	    /*
	      if (downprobnu.size() != 0) {
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree6", "nu_f", minprobnu_f));
	      }
	      else {
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree6", "nu_f", -1));
	      }

	      if (downprobcr.size() !=0) {
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree6", "cr_f", minprobcr_f));
	      }
	      else {
	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree6", "cr_f", -1));
	      }

	      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree6", "nuProb", nuProbability));
	      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttree6"));
	    */
      }
    // Sort the slices by probability
    std::sort(sliceIndexProbabilityPairs.begin(), sliceIndexProbabilityPairs.end(), [] (const UintFloatPair &a, const UintFloatPair &b)
	      {
		return (a.second > b.second);
	      });

    // Select the first m_maxNeutrinos as neutrinos, and the rest as cosmic
    std::cout << "m_maxNeutrinos " << m_maxNeutrinos << std::endl;
    unsigned int nNuSlices(0);
    for (const UintFloatPair &slice : sliceIndexProbabilityPairs)
      {
	if (nNuSlices < m_maxNeutrinos) // && maxprobcr_f < 0.55)
	  {
	    std::cout << "Calling this a neutrino!  " << slice.first << std::endl;

	    this->SelectPfos(nuSliceHypotheses.at(slice.first), selectedPfos);
	    nNuSlices++;
	    continue;
	  }

	std::cout << "Calling this a cosmic!" << std::endl;
	this->SelectPfos(crSliceHypotheses.at(slice.first), selectedPfos);
      }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectPfos(const PfoList &pfos, PfoList &selectedPfos) const
{
    selectedPfos.insert(selectedPfos.end(), pfos.begin(), pfos.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

// TODO think about how to make this function cleaner when features are more established
  NeutrinoIdTool::SliceFeatures::SliceFeatures(const PfoList &nuPfos, const PfoList &crPfos, const NeutrinoIdTool *const pTool, const PfoToFloatMap &pfotoprobabilitymapb) :
    m_isAvailable(false),
    m_pTool(pTool)
{
    try
    {
        const ParticleFlowObject *const pNeutrino(this->GetNeutrino(nuPfos));
        const CartesianVector &nuVertex(LArPfoHelper::GetVertex(pNeutrino)->GetPosition());
        const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());

        // Neutrino features
        CartesianVector nuWeightedDirTotal(0.f, 0.f, 0.f);
        unsigned int nuNHitsUsedTotal(0);
        unsigned int nuNHitsTotal(0);
        CartesianPointVector nuAllSpacePoints;
        for (const ParticleFlowObject *const pPfo : nuFinalStates)
        {
	  //add probability feature
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nuAllSpacePoints.insert(nuAllSpacePoints.end(), spacePoints.begin(), spacePoints.end());
            nuNHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5) continue;

            const CartesianVector dir(this->GetDirectionFromVertex(spacePoints, nuVertex));
            nuWeightedDirTotal += dir * static_cast<float>(spacePoints.size());
            nuNHitsUsedTotal += spacePoints.size();
        }

        if (nuNHitsUsedTotal == 0) return;
        const CartesianVector nuWeightedDir(nuWeightedDirTotal * (1.f / static_cast<float>(nuNHitsUsedTotal)));

        CartesianPointVector pointsInSphere;
        this->GetPointsInSphere(nuAllSpacePoints, nuVertex, 10, pointsInSphere);

        CartesianVector centroid(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenValues eigenValues(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenVectors eigenVectors;
        LArPcaHelper::RunPca(pointsInSphere, centroid, eigenValues, eigenVectors);


        const float nuNFinalStatePfos(static_cast<float>(nuFinalStates.size()));
        const float nuVertexY(nuVertex.GetY());
        const float nuWeightedDirZ(nuWeightedDir.GetZ());
        const float nuNSpacePointsInSphere(static_cast<float>(pointsInSphere.size()));

        if (eigenValues.GetX() <= std::numeric_limits<float>::epsilon()) return;
        const float nuEigenRatioInSphere(eigenValues.GetY() / eigenValues.GetX());

	std::vector<float> downprobnu;
	float minprobnu_f = -1;
	float maxprobnu_f = -1;
	for (const ParticleFlowObject *const pPfo : nuPfos)
	  {

	    PfoList daughterPfos = pPfo->GetDaughterPfoList();
	    
	    for (const ParticleFlowObject *const pPPfo : daughterPfos) {
	      auto search = pfotoprobabilitymapb.find(pPPfo);
	      if (search !=  pfotoprobabilitymapb.end()) {
		std::cout << "Found " << search->first << " " << search->second << '\n';
		if(search->second != -1 && search->second != -2 && search->second != -3 && search->second != -4) {
		  downprobnu.push_back(search->second);
		}
	      }
	    }

	  }


	if(downprobnu.size() > 0) {
	  std::vector<float>::iterator minprobnu = std::min_element(downprobnu.begin(), downprobnu.end());
	  auto indexn = std::distance(downprobnu.begin(), minprobnu);   
	  int indexvaluen = indexn;
	  minprobnu_f = downprobnu[indexvaluen];
	  std::vector<float>::iterator maxprobnu = std::max_element(downprobnu.begin(), downprobnu.end());
	  auto index2n = std::distance(downprobnu.begin(), maxprobnu);   
	  int indexvalue2n = index2n;
	  maxprobnu_f = downprobnu[indexvalue2n];
	}

        // Cosmic-ray features
        unsigned int nCRHitsMax(0);
        unsigned int nCRHitsTotal(0);
        float crLongestTrackDirY(std::numeric_limits<float>::max());
        float crLongestTrackDeflection(-std::numeric_limits<float>::max());

        for (const ParticleFlowObject *const pPfo : crPfos)
        {
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nCRHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5) continue;

            if (spacePoints.size() > nCRHitsMax)
            {
                nCRHitsMax = spacePoints.size();
                const CartesianVector upperDir(this->GetUpperDirection(spacePoints));
                const CartesianVector lowerDir(this->GetLowerDirection(spacePoints));

                crLongestTrackDirY = upperDir.GetY();
                crLongestTrackDeflection = 1.f - upperDir.GetDotProduct(lowerDir * (-1.f));
            }
        }

        if (nCRHitsMax == 0) return;
        if (nCRHitsTotal == 0) return;

        const float crFracHitsInLongestTrack = static_cast<float>(nCRHitsMax)/static_cast<float>(nCRHitsTotal);

	std::vector<float> downprobcr;
	float minprobcr_f = -1;
	float maxprobcr_f = -1;
	for (const ParticleFlowObject *const pPfo : crPfos)
	  {
	    	    
	    auto search = pfotoprobabilitymapb.find(pPfo);
	    if (search !=  pfotoprobabilitymapb.end()) {
	      std::cout << "Found " << search->first << " " << search->second << '\n';
	      if(search->second != -1 && search->second != -2  && search->second != -3 && search->second != -4) {
		downprobcr.push_back(search->second);
	      }
	    }
	   
	  }

	if(downprobcr.size() > 0) {
	  std::vector<float>::iterator minprobcr = std::min_element(downprobcr.begin(), downprobcr.end());
	  auto index = std::distance(downprobcr.begin(), minprobcr);   
	  int indexvalue = index;
	  minprobcr_f = downprobcr[indexvalue];
	  std::vector<float>::iterator maxprobcr = std::max_element(downprobcr.begin(), downprobcr.end());
	  auto index2 = std::distance(downprobcr.begin(), maxprobcr);   
	  int indexvalue2 = index2;
	  maxprobcr_f = downprobcr[indexvalue2];
	}

	// Push the features to the feature vector
	m_featureVector.push_back(nuNFinalStatePfos);
	m_featureVector.push_back(nuNHitsTotal);
	m_featureVector.push_back(nuVertexY);
	m_featureVector.push_back(nuWeightedDirZ);
	m_featureVector.push_back(nuNSpacePointsInSphere);
	m_featureVector.push_back(nuEigenRatioInSphere);
	m_featureVector.push_back(crLongestTrackDirY);
	m_featureVector.push_back(crLongestTrackDeflection);
	m_featureVector.push_back(crFracHitsInLongestTrack);
        m_featureVector.push_back(nCRHitsMax);
	m_featureVector.push_back(minprobnu_f);
	m_featureVector.push_back(maxprobnu_f);
	m_featureVector.push_back(minprobcr_f);
	m_featureVector.push_back(maxprobcr_f);

        m_isAvailable = true;
    }
    catch (StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    featureVector.insert(featureVector.end(), m_featureVector.begin(), m_featureVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float NeutrinoIdTool::SliceFeatures::GetNeutrinoProbability(const SupportVectorMachine &supportVectorMachine) const
{
    // ATTN if one or more of the features can not be calculated, then default to calling the slice a cosmic ray
  std::cout << "GetNeutrinoProb in NeutrinoIDTool" << std::endl;
    if (!this->IsFeatureVectorAvailable()) return 0.f;

    LArMvaHelper::MvaFeatureVector featureVector;
    std::cout << "GetFeatureVector" << std::endl;
    this->GetFeatureVector(featureVector);
    std::cout << "Go To LArMvaHelper" << std::endl;
    return LArMvaHelper::CalculateProbability(supportVectorMachine, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *NeutrinoIdTool::SliceFeatures::GetNeutrino(const PfoList &nuPfos) const
{
    // ATTN we should only ever have one neutrino reconstructed per slice
    if (nuPfos.size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return nuPfos.front();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetSpacePoints(const ParticleFlowObject *const pPfo, CartesianPointVector &spacePoints) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    if (clusters3D.size() > 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    if (clusters3D.empty()) return;

    CaloHitList caloHits;
    clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHits);

    for (const CaloHit *const pCaloHit : caloHits)
        spacePoints.push_back(pCaloHit->GetPositionVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetDirectionFromVertex(const CartesianPointVector &spacePoints, const CartesianVector &vertex) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return ((pointA - vertex).GetMagnitude() < (pointB - vertex).GetMagnitude());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetUpperDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return (pointA.GetY() > pointB.GetY());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetLowerDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return (pointA.GetY() < pointB.GetY());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetDirection(const CartesianPointVector &spacePoints, std::function<bool(const CartesianVector &pointA, const CartesianVector &pointB)> fShouldChooseA) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(m_pTool->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    const ThreeDSlidingFitResult fit(&spacePoints, 5, layerPitch);
    const CartesianVector endMin(fit.GetGlobalMinLayerPosition());
    const CartesianVector endMax(fit.GetGlobalMaxLayerPosition());
    const CartesianVector dirMin(fit.GetGlobalMinLayerDirection());
    const CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());

    const bool isMinStart(fShouldChooseA(endMin, endMax));
    const CartesianVector startPoint(isMinStart ? endMin : endMax);
    const CartesianVector endPoint(isMinStart ? endMax : endMin);
    const CartesianVector startDir(isMinStart ? dirMin : dirMax);

    const bool shouldFlip((endPoint - startPoint).GetUnitVector().GetDotProduct(startDir) < 0.f);
    return (shouldFlip ? startDir*(-1.f) : startDir);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetPointsInSphere(const CartesianPointVector &spacePoints, const CartesianVector &vertex, const float radius, CartesianPointVector &spacePointsInSphere) const
{
    for (const CartesianVector &point : spacePoints)
    {
        if ((point - vertex).GetMagnitudeSquared() <= radius*radius)
            spacePointsInSphere.push_back(point);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectNuanceCode", m_selectNuanceCode));

    if (m_selectNuanceCode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "NuanceCode", m_nuance));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumNeutrinoProbability", m_minProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaximumNeutrinos", m_maxNeutrinos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    if (!m_useTrainingMode)
    {
        std::string svmName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmName", svmName));

        std::string svmFileName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmFileName", svmFileName));

        const std::string fullSvmFileName(LArFileHelper::FindFileInPath(svmFileName, m_filePathEnvironmentVariable));
        m_supportVectorMachine.Initialize(fullSvmFileName, svmName);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
