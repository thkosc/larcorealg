/**
 * @file   larcorealg/Geometry/LocalTransformation.tcc
 * @brief  Class containing local-to-world transformations
 *         (template implementation)
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 30, 2016
 * @see    LocalTransformation.h
 *
 * This file is expected to be included directly in the header.
 *
 */

#ifndef LARCOREALG_GEOMETRY_LOCALTRANSFORMATION_TCC
#define LARCOREALG_GEOMETRY_LOCALTRANSFORMATION_TCC

// framework libraries
#include "cetlib_except/exception.h"

// ROOT
#include "TGeoNode.h"
#include "TGeoMatrix.h"

// CLHEP
#include "CLHEP/Geometry/Transform3D.h"

// C standard library
#include <cassert>


//------------------------------------------------------------------------------
namespace geo {
  namespace details {
    
    template <typename T, std::size_t SrcN = 3, std::size_t DestN = SrcN>
    bool doBuffersOverlap(T const* src, T const* dest)
      { return (dest < (src + SrcN)) && (src < (dest + DestN)); }
    
    template <typename T, std::size_t SrcN = 3, std::size_t DestN = SrcN>
    void checkVectorBufferOverlap(T const* src, T const* dest) {
      if (doBuffersOverlap<T, SrcN, DestN>(src, dest)) {
        throw cet::exception("LocalTransformation")
          << "source " << SrcN << "@[" << ((void*) src) << "]"
          << " and destination " << DestN << "@[" << ((void*) dest) << "]"
          << " buffers overlap!\n";
      }
      assert(!doBuffersOverlap(src, dest));
    } // checkVectorBufferOverlap()
    
  } // namespace details
} // namespace geo


//------------------------------------------------------------------------------
template <typename Matrix>
void geo::LocalTransformation<Matrix>::LocalToWorld
  (double const* local, double* world) const
{
  details::checkVectorBufferOverlap(local, world);
  fGeoMatrix.LocalToMaster(local, world); 
} // geo::LocalTransformation::LocalToWorld()


//------------------------------------------------------------------------------
template <typename Matrix>
void geo::LocalTransformation<Matrix>::LocalToWorldVect
  (double const* local, double* world) const
{
  details::checkVectorBufferOverlap(local, world);
  fGeoMatrix.LocalToMasterVect(local, world);
} // geo::LocalTransformation::LocalToWorldVect()


//------------------------------------------------------------------------------
template <typename Matrix>
void geo::LocalTransformation<Matrix>::WorldToLocal
  (double const* world, double* local) const
{
  details::checkVectorBufferOverlap(local, world);
  fGeoMatrix.MasterToLocal(world, local);
} // geo::LocalTransformation::WorldToLocal()


//------------------------------------------------------------------------------
template <typename Matrix>
void geo::LocalTransformation<Matrix>::WorldToLocalVect
  (const double* world, double* local) const
{
  details::checkVectorBufferOverlap(local, world);
  fGeoMatrix.MasterToLocalVect(world, local);
} // geo::LocalTransformation::WorldToLocalVect()


//------------------------------------------------------------------------------
template <typename Matrix>
template <typename DestPoint, typename SrcPoint>
DestPoint geo::LocalTransformation<Matrix>::WorldToLocalImpl
  (SrcPoint const& world) const
{
  double const worldArray[3] = { world.X(), world.Y(), world.Z() };
  double localArray[3];
  WorldToLocal(worldArray, localArray);
  return { localArray[0], localArray[1], localArray[2] };
} // geo::LocalTransformation::WorldToLocal()


//......................................................................
template <typename Matrix>
template <typename DestVector, typename SrcVector>
DestVector geo::LocalTransformation<Matrix>::WorldToLocalVectImpl
  (SrcVector const& world) const
{
  double const worldArray[3] = { world.X(), world.Y(), world.Z() };
  double localArray[3];
  WorldToLocalVect(worldArray, localArray);
  return { localArray[0], localArray[1], localArray[2] };
} // geo::LocalTransformation::WorldToLocalVect()


//......................................................................
template <typename Matrix>
template <typename DestPoint, typename SrcPoint>
DestPoint geo::LocalTransformation<Matrix>::LocalToWorldImpl
  (SrcPoint const& local) const
{
  double const localArray[3] = { local.X(), local.Y(), local.Z() };
  double worldArray[3];
  LocalToWorld(localArray, worldArray);
  return { worldArray[0], worldArray[1], worldArray[2] };
} // geo::LocalTransformation::LocalToWorld()

  
//......................................................................
template <typename Matrix>
template <typename DestVector, typename SrcVector>
DestVector geo::LocalTransformation<Matrix>::LocalToWorldVectImpl
  (SrcVector const& local) const
{
  double const localArray[3] = { local.X(), local.Y(), local.Z() };
  double worldArray[3];
  LocalToWorldVect(localArray, worldArray);
  return { worldArray[0], worldArray[1], worldArray[2] };
} // geo::LocalTransformation::LocalToWorldVect()


//------------------------------------------------------------------------------
// specialisations (forward declarations)
// 
namespace geo {
  
  template <>
  TGeoHMatrix LocalTransformation<TGeoHMatrix>::transformationFromPath
    (std::vector<TGeoNode const*> const& path, size_t depth);

  template <>
  HepGeom::Transform3D
  LocalTransformation<HepGeom::Transform3D>::transformationFromPath
    (std::vector<TGeoNode const*> const& path, size_t depth);

} // namespace geo

//------------------------------------------------------------------------------
  
#endif // LARCOREALG_GEOMETRY_LOCALTRANSFORMATION_TCC
