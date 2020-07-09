/**
 *  @file   larpandoracontent/LArControlFlow/NDQUalitiesAlg.h
 *
 *  @brief  Header file for the neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ND_QUALITIES_ALG_H
#define LAR_ND_QUALITIES_ALG_H 1

#include "Pandora/ExternallyConfiguredAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

//#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "Pandora/Algorithm.h"
#include <string>
#include <vector>
#include "Objects/CaloHit.h"

namespace lar_content
{

  class NDQualitiesAlg : public pandora::Algorithm
  {
  public:
    NDQualitiesAlg();
    virtual ~NDQualitiesAlg();

     class HitCharge
    {
    public:

        HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty);
        HitCharge();

        const pandora::CaloHit* GetCaloHit() const;
        float GetLongitudinalPosition() const;
        float GetHitWidth() const;
        float GetCharge() const;
        float GetChargeOverWidth() const;
        float GetUncertainty() const;

        void SetDistanceToNN(float &distance);
        float GetDistanceToNN() const;

        void SetForwardsFitCharge(float &Q_fit_f); 
        void SetForwardsSigma(float &f_sigma);
        void SetForwardsDelta(float &forwardsDelta);
        void SetForwardsChiSquared(float &forwardsHitChisquared);

        void SetBackwardsFitCharge(float &Q_fit_b); 
        void SetBackwardsSigma(float &b_sigma);
        void SetBackwardsDelta(float &backwardsDelta);
        void SetBackwardsChiSquared(float &backwardsHitChisquared);

        float GetForwardsFitCharge(); 
        float GetForwardsSigma();
        float GetForwardsDelta();
        float GetForwardsChiSquared();

        float GetBackwardsFitCharge(); 
        float GetBackwardsSigma();
        float GetBackwardsDelta();
        float GetBackwardsChiSquared();

        bool                                       m_intails;

    private:
        const pandora::CaloHit*                    m_calohit;
        float                                      m_longitudinalposition;    ///<
        float                                      m_hitwidth;                ///<
        float                                      m_charge;                  ///<
        float                                      m_qoverx;                  ///<
        float                                      m_uncertainty;          ///<
        float                                      m_distancetonearestneighbour;          ///<

        float                                      m_forwardsfitcharge;
        float                                      m_forwardssigma;
        float                                      m_forwardsdelta;
        float                                      m_forwardschisquared;

        float                                      m_backwardsfitcharge;
        float                                      m_backwardssigma;
        float                                      m_backwardsdelta;
        float                                      m_backwardschisquared;
    };

    typedef std::vector<HitCharge> HitChargeVector;

     class LookupTable
    {
    public:

        LookupTable();

        LookupTable(double &initialEnergy, double &binWidth);

        std::map<int, double> GetMap();

        void SetMap(std::map<int, double> &map);

        std::map<double, int> GetReverseMap();

        void SetReverseMap(std::map<double, int> &map);

        double GetInitialEnergy();

        void SetInitialEnergy(double &initialEnergy);

        double GetBinWidth();

        void SetBinWidth(double &binWidth);

        void SetMaxRange(double &maxRange);

        double GetMaxRange();

    private:
        std::map<int, double>                       m_map;
        std::map<double, int>                       m_reversemap;
        double                                      m_binwidth;               ///<
        double                                      m_initialenergy;
        double                                      m_maxrange;
    };



   
  private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    std::string                     m_pfoListName;
    unsigned int            m_slidingFitWindow; 
    TwoDSlidingFitResultMap m_slidingFitResultMap;              ///< The sliding fit result map
    unsigned int            m_minClusterCaloHits;               ///< The min number of hits in base cluster selection method
    float                   m_minClusterLength;          ///< The min length (squared) in base cluster selection method

    int                     m_targetParticlePDG;
    int                     m_numberTrackEndHits;
    bool                    m_enableFragmentRemoval;
    bool                    m_enableSplitting;

    double                  m_tableInitialEnergy;
    double                  m_tableStepSize;

    bool                    m_writeTable;

    std::string             m_lookupTableFileName;
    std::string             m_probabilityFileName;
    std::string             m_treeName;
    

    void FillHitChargeVector(const pandora::Cluster *const pCluster, HitChargeVector &hitChargeVector, pandora::CartesianVector&, pandora::CartesianVector&);

    void GetTrackLength(HitChargeVector &hitChargeVector, float &trackLength);
    
    const TwoDSlidingFitResult &GetCachedSlidingFit(const pandora::Cluster *const pCluster) const;

    static bool SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2);

    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    void BinHitChargeVector(lar_content::NDQualitiesAlg::HitChargeVector &hitChargeVector, lar_content::NDQualitiesAlg::HitChargeVector &binnedHitChargeVector);

    double GetEnergyfromLength(lar_content::NDQualitiesAlg::LookupTable &lookupTable, double &trackLength);

    double GetLengthfromEnergy(lar_content::NDQualitiesAlg::LookupTable &lookupTable, double &currentEnergy);

    double DensityCorrection(double &T, double &M);

    double BetheBloch(double &T, double &M);

    void FillLookupTable(lar_content::NDQualitiesAlg::LookupTable& lookupTable, double M);

    
    
  };



inline NDQualitiesAlg::HitCharge::HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty) :
    m_intails((hitCharge/hitWidth) <= 1.4 ? true : false),
    m_calohit(caloHit),
    m_longitudinalposition(longitudinalPosition),
    m_hitwidth(hitWidth),
    m_charge(hitCharge),
    m_qoverx(hitCharge/hitWidth),
    m_uncertainty(uncertainty),
    m_forwardsfitcharge(0.f),
    m_forwardssigma(0.f),
    m_forwardsdelta(0.f),
    m_forwardschisquared(0.f),
    m_backwardsfitcharge(0.f),
    m_backwardssigma(0.f),
    m_backwardsdelta(0.f),
    m_backwardschisquared(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDQualitiesAlg::HitCharge::HitCharge() :
    m_intails(false),
    m_calohit(NULL),
    m_longitudinalposition(0.f),
    m_hitwidth(0.f),
    m_charge(0.f),
    m_qoverx(0.f),
    m_uncertainty(0.f),
    m_forwardsfitcharge(0.f),
    m_forwardssigma(0.f),
    m_forwardsdelta(0.f),
    m_forwardschisquared(0.f),
    m_backwardsfitcharge(0.f),
    m_backwardssigma(0.f),
    m_backwardsdelta(0.f),
    m_backwardschisquared(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit* NDQualitiesAlg::HitCharge::GetCaloHit() const
{
    return m_calohit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetLongitudinalPosition() const
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetHitWidth() const
{
    return m_hitwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetCharge() const
{
    return m_charge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetChargeOverWidth() const
{
    return m_qoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetUncertainty() const
{
    return m_uncertainty;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::HitCharge::SetDistanceToNN(float &distance)
{
    m_distancetonearestneighbour = distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetDistanceToNN() const
{
    return m_distancetonearestneighbour;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetForwardsFitCharge(float &Q_fit_f) 
{
    m_forwardsfitcharge = Q_fit_f;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetForwardsSigma(float &f_sigma)
{
    m_forwardssigma = f_sigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetForwardsDelta(float &forwardsDelta)
{
    m_forwardsdelta = forwardsDelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetForwardsChiSquared(float &forwardsHitChisquared)
{
    m_forwardschisquared = forwardsHitChisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::HitCharge::SetBackwardsFitCharge(float &Q_fit_b) 
{
    m_backwardsfitcharge = Q_fit_b;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetBackwardsSigma(float &b_sigma)
{
    m_backwardssigma= b_sigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetBackwardsDelta(float &backwardsDelta)
{
    m_backwardsdelta = backwardsDelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void NDQualitiesAlg::HitCharge::SetBackwardsChiSquared(float &backwardsHitChisquared)
{
    m_backwardschisquared = backwardsHitChisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetForwardsFitCharge() 
{
    return m_forwardsfitcharge;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetForwardsSigma()
{
    return m_forwardssigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetForwardsDelta()
{
    return m_forwardsdelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetForwardsChiSquared()
{
    return m_forwardschisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float NDQualitiesAlg::HitCharge::GetBackwardsFitCharge() 
{
    return m_backwardsfitcharge;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetBackwardsSigma()
{
    return m_backwardssigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetBackwardsDelta()
{
    return m_backwardsdelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float NDQualitiesAlg::HitCharge::GetBackwardsChiSquared()
{
    return m_backwardschisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline NDQualitiesAlg::LookupTable::LookupTable()
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = 0.f;
    m_initialenergy = 0.f,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDQualitiesAlg::LookupTable::LookupTable(double &initialEnergy, double &binWidth)
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = binWidth;
    m_initialenergy = initialEnergy,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<int, double> NDQualitiesAlg::LookupTable::GetMap()
{
    return m_map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::LookupTable::SetMap(std::map<int, double> &map)
{
    m_map = map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<double, int> NDQualitiesAlg::LookupTable::GetReverseMap()
{
    return m_reversemap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::LookupTable::SetReverseMap(std::map<double, int> &reverseMap)
{
    m_reversemap = reverseMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double NDQualitiesAlg::LookupTable::GetInitialEnergy()
{
    return m_initialenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::LookupTable::SetInitialEnergy(double &initialEnergy)
{
    m_initialenergy = initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double NDQualitiesAlg::LookupTable::GetBinWidth()
{
    return m_binwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::LookupTable::SetBinWidth(double &binWidth)
{
    m_binwidth = binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDQualitiesAlg::LookupTable::SetMaxRange(double &maxRange)
{
    m_maxrange = maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double NDQualitiesAlg::LookupTable::GetMaxRange()
{
    return m_maxrange;
}

//------------------------------------------------------------------------------------------------------------------------------------------



}
#endif //POSTMASTERANALYSIS_H
