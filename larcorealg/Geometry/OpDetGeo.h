////////////////////////////////////////////////////////////////////////
/// \file  OpDetGeo.h
/// \brief Encapsulate the geometry of an optical detector
///
/// \author  bjpjones@mit.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_OPDETGEO_H
#define LARCOREALG_GEOMETRY_OPDETGEO_H

// ROOT libraries
#include "TGeoMatrix.h" // TGeoHMatrix

// C/C++ standard libraries
#include <vector>


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

  private:
    const TGeoNode* fOpDetNode;  ///< Pointer to theopdet node
    TGeoHMatrix     fGeoMatrix; ///< Transformation matrix to world frame
  };
}


#endif // LARCOREALG_GEOMETRY_OPDETGEO_H
