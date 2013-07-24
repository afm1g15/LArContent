/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "LArClusterAssociation/ClusterExtensionAlgorithm.h"
#include "LArClusterAssociation/IsolatedHitMergingAlgorithm.h"
#include "LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArClustering/ClusterCreationAlgorithm.h"
#include "LArClustering/ClusteringParentAlgorithm.h"
#include "LArClustering/RemnantClusteringAlgorithm.h"
#include "LArClusterSeedAssociation/BoundedClusterMergingAlgorithm.h"
#include "LArClusterSeedAssociation/ConeBasedMergingAlgorithm.h"
#include "LArClusterSeedAssociation/ParallelClusterMergingAlgorithm.h"
#include "LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/NtupleWritingAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"
#include "LArReclustering/ShowerMipSeparationAlgorithm.h"
#include "LArReclustering/TrackSplittingAlgorithm.h"
#include "LArThreeDSeed/ThreeDShowersAlgorithm.h"
#include "LArThreeDSeed/ThreeDTracksAlgorithm.h"
#include "LArTwoDSeed/SeedBranchGrowingAlgorithm.h"
#include "LArTwoDSeed/SeedConsolidationAlgorithm.h"
#include "LArTwoDSeed/SeedFindingAlgorithm.h"
#include "LArTwoDSeed/SeedLengthGrowingAlgorithm.h"
#include "LArTwoDSeed/SeedRelegationAlgorithm.h"
#include "LArTwoDSeed/VertexSeedFindingAlgorithm.h"
#include "LArUtility/EventPreparationAlgorithm.h"
#include "LArUtility/TwoDPreparationAlgorithm.h"
#include "LArVertex/VertexFindingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArBFieldCalculator.h"
#include "LArPseudoLayerCalculator.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArClusterExtension",                    lar::ClusterExtensionAlgorithm::Factory)                                    \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArClusterCreation",                     lar::ClusterCreationAlgorithm::Factory)                                     \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArRemnantClustering",                   lar::RemnantClusteringAlgorithm::Factory)                                   \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArParallelClusterMerging",              lar::ParallelClusterMergingAlgorithm::Factory)                              \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArLongitudinalAssociation",             lar::LongitudinalAssociationAlgorithm::Factory)                             \
        d("LArEventDisplay",                        lar::EventDisplayAlgorithm::Factory)                                        \
        d("LArNtupleWriting",                       lar::NtupleWritingAlgorithm::Factory)                                       \
        d("LArVisualMonitoring",                    lar::VisualMonitoringAlgorithm::Factory)                                    \
        d("LArShowerMipSeparation",                 lar::ShowerMipSeparationAlgorithm::Factory)                                 \
        d("LArTrackSplitting",                      lar::TrackSplittingAlgorithm::Factory)                                      \
        d("LArThreeDShowers",                       lar::ThreeDShowersAlgorithm::Factory)                                       \
        d("LArThreeDTracks",                        lar::ThreeDTracksAlgorithm::Factory)                                        \
        d("LArSeedBranchGrowing",                   lar::SeedBranchGrowingAlgorithm::Factory)                                   \
        d("LArSeedConsolidation",                   lar::SeedConsolidationAlgorithm::Factory)                                   \
        d("LArSeedFinding",                         lar::SeedFindingAlgorithm::Factory)                                         \
        d("LArSeedLengthGrowing",                   lar::SeedLengthGrowingAlgorithm::Factory)                                   \
        d("LArSeedRelegation",                      lar::SeedRelegationAlgorithm::Factory)                                      \
        d("LArTransverseAssociation",               lar::TransverseAssociationAlgorithm::Factory)                               \
        d("LArVertexSeedFinding",                   lar::VertexSeedFindingAlgorithm::Factory)                                   \
        d("LArEventPreparation",                    lar::EventPreparationAlgorithm::Factory)                                    \
        d("LArTwoDPreparation",                     lar::TwoDPreparationAlgorithm::Factory)                                     \
        d("LArVertexFinding",                       lar::VertexFindingAlgorithm::Factory)                                       \
        d("LArVertexSplitting",                     lar::VertexSplittingAlgorithm::Factory)

    #define LAR_PARTICLE_ID_LIST(d)                                                                                             \
        d("LArEmShowerId",                          &lar::LArParticleIdHelper::LArEmShowerId)                                   \
        d("LArPhotonId",                            &lar::LArParticleIdHelper::LArPhotonId)                                     \
        d("LArElectronId",                          &lar::LArParticleIdHelper::LArElectronId)                                   \
        d("LArMuonId",                              &lar::LArParticleIdHelper::LArMuonId)

    #define LAR_SETTINGS_LIST(d)                                                                                                \
        d("LArClusterHelper",                       &lar::LArClusterHelper::ReadSettings)                                       \
        d("LArGeometryHelper",                      &lar::LArGeometryHelper::ReadSettings)                                      \
        d("LArParticleIdHelper",                    &lar::LArParticleIdHelper::ReadSettings)                                    \
        d("LArPointingClusterHelper",               &lar::LArPointingClusterHelper::ReadSettings)                               \
        d("LArVertexHelper",                        &lar::LArVertexHelper::ReadSettings)

    /**
     *  @brief  Register all the fine granularity algorithms with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(pandora::Pandora &pandora);

    /**
     *  @brief  Register all the fine granularity helper functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterHelperFunctions(pandora::Pandora &pandora);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterAlgorithms(pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(PANDORA_REGISTER_ALGORITHM);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterHelperFunctions(pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);
    LAR_SETTINGS_LIST(PANDORA_REGISTER_SETTINGS);

    return pandora::STATUS_CODE_SUCCESS;
}

#endif // #ifndef LAR_CONTENT_H
