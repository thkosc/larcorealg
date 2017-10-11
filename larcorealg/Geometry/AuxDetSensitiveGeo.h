////////////////////////////////////////////////////////////////////////
/// \file  AuxDetSensitiveGeo.h
/// \brief Encapsulate the geometry of the sensitive portion of an auxiliary detector
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GEO_AuxDetSensitiveGeo_H
#define GEO_AuxDetSensitiveGeo_H
#include <vector>
#include "TGeoVolume.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::array(), ...

class TGeoNode;
class TGeoHMatrix;
class TGeoMatrix;

namespace geo {
  
  class AuxDetSensitiveGeo {
  public:
    AuxDetSensitiveGeo(std::vector<const TGeoNode*>& path, 
		       int                           depth);
    AuxDetSensitiveGeo(const TGeoVolume* volume, 
		       TGeoHMatrix*       rotation);
    ~AuxDetSensitiveGeo();
    
    void GetCenter(double* xyz, 
		   double localz=0.0)    const;
    void GetNormalVector(double* xyzDir) const;

    //box geometry
    double Length()                  const { return fLength;      }
    double HalfWidth1()      	     const { return fHalfWidth1;  }
    double HalfWidth2()      	     const { return fHalfWidth2;  }
    double HalfHeight()      	     const { return fHalfHeight;  }
    const TGeoVolume* TotalVolume()  const { return fTotalVolume; }
    
    double DistanceToPoint(double * xyz) const;
    
    void LocalToWorld    (const double* local, double* world) const;
    void LocalToWorldVect(const double* local, double* world) const;
    void WorldToLocal    (const double* world, double* local) const;
    void WorldToLocalVect(const double* world, double* local) const;
    
    /**
     * @brief Prints information about this auxiliary sensitive detector.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     * 
     * Note that the first line out the output is _not_ indented.
     * 
     * Verbosity levels
     * -----------------
     * 
     * * 0 _(default)_: only center
     * * 1: also size
     * * 2: also normal direction
     * 
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintAuxDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;
    
    /// Maximum verbosity supported by `PrintAuxDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;
    
  private:

    TGeoHMatrix*                  fGeoMatrix;   ///< Transformation matrix to world frame		      
    const TGeoVolume*     	  fTotalVolume; ///< Total volume of AuxDet, called vol*		      
    double                	  fLength;      ///< length of volume, along z direction in local	      
    double                	  fHalfWidth1;  ///< 1st half width of volume, at -z/2 in local coordinates 
    double                	  fHalfWidth2;  ///< 2nd half width (width1==width2 for boxes), at +z/2     
    double                	  fHalfHeight;  ///< half height of volume                                  
  };
}


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::AuxDetSensitiveGeo::PrintAuxDetInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 0 */
) const {
  
  //----------------------------------------------------------------------------
  std::array<double, 3U> center;
  GetCenter(center.data());
  out << "centered at " << lar::dump::array<3U>(center) << " cm";
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  lar::util::RealComparisons<double> coordIs(1e-4);
  out << ", size ( " << (2.0 * HalfWidth1());
  if (coordIs.nonEqual(HalfWidth1(), HalfWidth2()))
    out << "/" << (2.0 * HalfWidth2());
  out << " x " << (2.0 * HalfHeight()) << " x " << Length() << " ) cm";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  std::array<double, 3U> normal;
  GetNormalVector(normal.data());
  out << ", normal facing " << lar::dump::array<3U>(normal);
  
//  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  
} // geo::AuxDetSensitiveGeo::PrintAuxDetInfo()


//------------------------------------------------------------------------------

#endif
