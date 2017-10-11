////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeo.h
/// \brief Encapsulate the geometry of an auxiliary detector
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GEO_AuxDetGeo_H
#define GEO_AuxDetGeo_H
#include <vector>
#include "TGeoVolume.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::array(), ...

class TGeoNode;
class TGeoHMatrix;
class TGeoMatrix;

namespace geo {
  
  class AuxDetGeo {
  public:
    AuxDetGeo(std::vector<const TGeoNode*>& path, 
              int                           depth);
    ~AuxDetGeo();
    
    void GetCenter(double* xyz, 
                   double localz=0.0)    const;
    void GetNormalVector(double* xyzDir) const;

    //box geometry
    double Length()                  const { return fLength;      }
    double HalfWidth1()      	       const { return fHalfWidth1;  }
    double HalfWidth2()      	       const { return fHalfWidth2;  }
    double HalfHeight()      	       const { return fHalfHeight;  }
    const TGeoVolume* TotalVolume()  const { return fTotalVolume; }
    
    double DistanceToPoint(double * xyz) const;
    
    void LocalToWorld    (const double* local, double* world) const;
    void LocalToWorldVect(const double* local, double* world) const;
    void WorldToLocal    (const double* world, double* local) const;
    void WorldToLocalVect(const double* world, double* local) const;

    std::string Name() const { return fTotalVolume->GetName(); }

    // methods for the sensitive volumes in the aux det
    
    size_t             const  FindSensitiveVolume      (double const worldLoc[3])  const;
    AuxDetSensitiveGeo const& PositionToSensitiveVolume(double const worldLoc[3],
                                                        size_t     & sv)           const;
    AuxDetSensitiveGeo const& SensitiveVolume(size_t sv) const { return *fSensitive[sv];   }
    size_t             const  NSensitiveVolume()         const { return fSensitive.size(); }
    
    void                      SortSubVolumes(geo::GeoObjectSorter const& sorter);

    
    /**
     * @brief Prints information about this auxiliary detector.
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
     * * 0: only detector name
     * * 1 _(default)_: also center
     * * 2: also size
     * * 3: also number of sensitive detectors
     * * 4: also normal direction
     * 
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintAuxDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;
    
    /// Maximum verbosity supported by `PrintAuxDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 4;
    
    
  private:

    void FindAuxDetSensitive(std::vector<const TGeoNode*>& path,
                             unsigned int                  depth);
    void MakeAuxDetSensitive(std::vector<const TGeoNode*>& path,
                             int                           depth);

    const TGeoVolume*     	         fTotalVolume; ///< Total volume of AuxDet, called vol*
    TGeoHMatrix*                     fGeoMatrix;   ///< Transformation matrix to world frame
    double                	         fLength;      ///< length of volume, along z direction in local
    double                	         fHalfWidth1;  ///< 1st half width of volume, at -z/2 in local coordinates
    double                	         fHalfWidth2;  ///< 2nd half width (width1==width2 for boxes), at +z/2
    double                	         fHalfHeight;  ///< half height of volume
    std::vector<AuxDetSensitiveGeo*> fSensitive;   ///< sensitive volumes in the detector
  };
}


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::AuxDetGeo::PrintAuxDetInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {
  
  //----------------------------------------------------------------------------
  out << "\"" << Name() << "\"";
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  std::array<double, 3U> center;
  GetCenter(center.data());
  out << " centered at " << lar::dump::array<3U>(center) << " cm";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  lar::util::RealComparisons<double> coordIs(1e-4);
  out << ", size ( " << (2.0 * HalfWidth1());
  if (coordIs.nonEqual(HalfWidth1(), HalfWidth2()))
    out << "/" << (2.0 * HalfWidth2());
  out << " x " << (2.0 * HalfHeight()) << " x " << Length() << " ) cm";
  
  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  out << "\n" << indent
    << "with ";
  switch (NSensitiveVolume()) {
    case 0:  out << "no sensitive volume"; break;
    case 1:  out << "1 sensitive volume"; break;
    default: out << NSensitiveVolume() << " sensitive volumes"; break;
  } // switch
  
  if (verbosity-- <= 0) return; // 3
  
  //----------------------------------------------------------------------------
  std::array<double, 3U> normal;
  GetNormalVector(normal.data());
  out << ", normal facing " << lar::dump::array<3U>(normal);
  
//  if (verbosity-- <= 0) return; // 4
  
  //----------------------------------------------------------------------------
  
} // geo::AuxDetGeo::PrintAuxDetInfo()


//------------------------------------------------------------------------------

#endif
