////////////////////////////////////////////////////////////////////////
/// \file  OpDetGeo.h
/// \brief Encapsulate the geometry of an optical detector
///
/// \author  bjpjones@mit.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_OPDETGEO_H
#define LARCOREALG_GEOMETRY_OPDETGEO_H

#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::array(), ...

// ROOT libraries
#include "TGeoMatrix.h" // TGeoHMatrix

// C/C++ standard libraries
#include <vector>
#include <array>


// forward declarations
class TGeoNode;

namespace geo {

  class OpDetGeo {
  public:
    OpDetGeo(std::vector<const TGeoNode*>& path, 
	    int depth);

    void   GetCenter(double* xyz, double localz=0.0) const;
    double RMax() const;
    double HalfL() const;
    double RMin() const;
    double ThetaZ(bool degrees = false) const;  ///< returns angle of detector
                                                ///< with respect to z axis 
                                                ///< in the Y-Z plane, in 
                                                ///< radians by default
    double CosThetaFromNormal(double const* xyz) const;
    double DistanceToPoint(double const* xyz) const;


    void LocalToWorld(const double* local, double* world)     const;
    void LocalToWorldVect(const double* local, double* world) const;
    void WorldToLocal(const double* world, double* local)     const;
    void WorldToLocalVect(const double* world, double* local) const;

    const TGeoNode*     Node() const { return fOpDetNode; }

    /**
     * @brief Prints information about this optical detector.
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
     * * 2: also angle from z axis
     * 
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintOpDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;
    
    /// Maximum verbosity supported by `PrintOpDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;
    
  private:
    const TGeoNode* fOpDetNode;  ///< Pointer to theopdet node
    TGeoHMatrix     fGeoMatrix; ///< Transformation matrix to world frame
  };
}


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::OpDetGeo::PrintOpDetInfo(
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
  out << ", radius: " << RMax() << " cm";
  if (RMin() != 0.0) out << " (inner: " << RMin() << " cm)";
  out << ", length: " << (2.0*HalfL()) << " cm";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  out << ", theta(z): " << ThetaZ() << " rad";
  
//  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  
} // geo::OpDetGeo::PrintOpDetInfo()


#endif // LARCOREALG_GEOMETRY_OPDETGEO_H
