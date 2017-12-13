////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/CryostatGeo.h
/// \brief Encapsulate the construction of a single cyostat
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_CRYOSTATGEO_H
#define LARCOREALG_GEOMETRY_CRYOSTATGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// ROOT libraries
#include "TGeoVolume.h"
#include "TGeoMatrix.h" // TGeoHMatrix

// C/C++ standard libraries
#include <vector>
#include <string>

// forward declarations
class TGeoNode;

namespace geo {

  //......................................................................
  /// Geometry information for a single cryostat
  class CryostatGeo: public geo::BoxBoundedGeo {
    
      public:
    
    using GeoNodePath_t = geo::WireGeo::GeoNodePath_t;
    
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     * 
     * These types represents points and displacement vectors in the reference
     * frame defined in the cryostat geometry box from the GDML geometry
     * description.
     * 
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     * 
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::CryostatGeo` have the same type but are not
     * compatible.
     */
    
    /// Tag for vectors in the "local" GDML coordinate frame of the cryostat.
    struct CryostatGeoCoordinatesTag {};
    
    /// Type of points in the local GDML cryostat frame.
    using LocalPoint_t = geo::Point3DBase_t<CryostatGeoCoordinatesTag>;
    
    /// Type of displacement vectors in the local GDML cryostat frame.
    using LocalVector_t = geo::Vector3DBase_t<CryostatGeoCoordinatesTag>;
    
    ///@}
    
    
    /// Construct a representation of a single cryostat of the detector.
    CryostatGeo(std::vector<const TGeoNode*>& path, int depth);
    
    
    /// @{
    /// @name Cryostat geometry information
    
    /// Half width of the cryostat [cm]
    double            HalfWidth()                               const;
    /// Half height of the cryostat [cm]
    double            HalfHeight()                              const;
    /// Half height of the cryostat [cm]
    double            HalfLength()                              const;
    /// Full width of the cryostat [cm]
    double            Width()                                   const { return 2. * HalfWidth(); }
    /// Full height of the cryostat [cm]
    double            Height()                                  const { return 2. * HalfHeight(); }
    /// Length of the cryostat [cm]
    double            Length()                                  const { return 2. * HalfLength(); }
    /// Mass of the cryostat
    double            Mass()                                    const { return fVolume->Weight(); }
    /// Pointer to ROOT's volume descriptor
    const TGeoVolume* Volume()                                  const { return fVolume;           }
    
    /// @brief Returns boundaries of the cryostat (in centimetres).
    /// @return boundaries in a geo::BoxBoundedGeo
    geo::BoxBoundedGeo const& Boundaries() const
      { return BoundingBox(); }
    
    /// @brief Fills boundaries of the cryostat (in centimetres).
    /// @param boundaries filled as: [0] -x [1] +x [2] -y [3] +y [4] -z [5] +z
    void Boundaries(double* boundaries) const;
    
    
    /// Returns the geometrical center of the cryostat.
    geo::Point_t GetCenter() const
      { return Boundaries().Center(); }
    
    /// Returns the bounding box of this cryostat.
    geo::BoxBoundedGeo const& BoundingBox() const
      { return static_cast<geo::BoxBoundedGeo const&>(*this); }
    
    /// Returns the identifier of this cryostat.
    geo::CryostatID const& ID() const { return fID; }
    
    
    /**
     * @brief Prints information about this cryostat.
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
     * * 0: only cryostat ID
     * * 1 _(default)_: also center and size
     * * 2: also number of TPCs, optical detectors, and maximum wires per plane
     * *    and of planes for TPC
     * * 3: also information on bounding box
     * 
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintCryostatInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;
    
    /// Maximum verbosity supported by `PrintCryostatInfo()`.
    static constexpr unsigned int MaxVerbosity = 3;
    
    /// @}
    
    
    /// @{
    /// @name TPC access
    
    //@{
    /// Number of TPCs in this cryostat
    unsigned int      NTPC()                                    const { return fTPCs.size();      }
    unsigned int      NElements()                               const { return fTPCs.size();      }
    //@}
    
    //@{
    /**
     * @brief Returns whether a TPC with index itpc is present in this cryostat
     * @param itpc index of TPC in this cryostat
     * @return whether the TPC with index itpc is present in this cryostat
     */
    bool HasTPC(unsigned int itpc) const { return itpc < NTPC(); }
    bool HasElement(unsigned int itpc) const { return HasTPC(itpc); }
    //@}
    
    //@{
    /**
     * @brief Returns whether the TPC in tpcid is present in this cryostat
     * @param tpcid full TPC ID
     * @return whether the TPC in tpcid is present in this cryostat
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    bool HasTPC(geo::TPCID const& tpcid) const { return HasTPC(tpcid.TPC); }
    bool HasElement(geo::TPCID const& tpcid) const { return HasTPC(tpcid); }
    //@}
    
    /// Return the itpc'th TPC in the cryostat.
    /// @throws cet::exception (category "TPCOutOfRange") if no such TPC
    const TPCGeo&     TPC(unsigned int itpc)                    const;
    
    //@{
    /**
     * @brief Returns the TPC in tpcid from this cryostat
     * @param tpcid full TPC ID
     * @return a constant reference to the TPC in tpcid
     * @throws cet::exception (category "TPCOutOfRange") if no such TPC
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    const TPCGeo&     TPC(TPCID const& tpcid)                   const
      { return TPC(tpcid.TPC); }
    const TPCGeo&     GetElement(TPCID const& tpcid)            const
      { return TPC(tpcid); }
    //@}
    
    /**
     * @brief Returns the TPC number itpc from this cryostat
     * @param itpc the number of local TPC
     * @return a constant pointer to the TPC, or nullptr if it does not exist
     */
    TPCGeo const*     TPCPtr(unsigned int itpc)                    const
      { return HasTPC(itpc)? &(fTPCs[itpc]): nullptr; }
    
    //@{
    /**
     * @brief Returns the TPC in tpcid from this cryostat
     * @param tpcid full TPC ID
     * @return a constant pointer to the TPC, or nullptr if it does not exist
     *
     * The cryostat number in tpcid is ignored, as it is ignored whether tpcid
     * is invalid.
     */
    TPCGeo const*     TPCPtr(TPCID const& tpcid)                   const
      { return TPCPtr(tpcid.TPC); }
    TPCGeo const*     GetElementPtr(TPCID const& tpcid)            const
      { return TPCPtr(tpcid); }
    //@}
    
    /**
     * @brief Returns the index of the TPC at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param wiggle a small factor (like 1+epsilon) to avoid rounding errors
     * @return the TPC index, or UINT_MAX if no TPC is there
     */
    unsigned int FindTPCAtPosition(double const worldLoc[3],
                                   double const wiggle) const;
    
    /// Return the TPCGeo object containing the world position worldLoc
    /// @todo What if there is none?
    const TPCGeo&     PositionToTPC(double const  worldLoc[3],
				    unsigned int &tpc,
				    double const &wiggle)       const;
    
    /// Returns the largest number of planes among the TPCs in this cryostat
    unsigned int MaxPlanes() const;
    
    /// Returns the largest number of wires among the TPCs in this cryostat
    unsigned int MaxWires() const;
    
    /// @}
    
    
    /// @{
    /// @name Optical detector access
    
    /// Number of optical detectors in this TPC
    unsigned int      NOpDet()                                  const { return fOpDets.size();    }
    
    /// Return the iopdet'th optical detector in the cryostat
    const OpDetGeo&   OpDet(unsigned int iopdet)                const;

    /// Find the nearest opdet to point in this cryostat
    unsigned int GetClosestOpDet(double const* xyz)             const;
    
    /// Get name of opdet geometry element
    std::string  OpDetGeoName()                                 const { return fOpDetGeoName; }

    /// @}

    /// @{
    /// @name Coordinate transformation
    
    /// @{
    /// @name Coordinate transformation
    
    /// Transform point from local cryostat frame to world frame.
    void LocalToWorld(const double* cryo, double* world) const
      { fTrans.LocalToWorld(cryo, world); }
      
    /// Transform point from local cryostat frame to world frame.
    /// @deprecated This method breaks the distinction between local and global
    ///             vectors, since input and output vectors share the same type;
    ///             use the "official" vector types `geo::Point_t`,
    ///             `geo::CryostatGeo::LocalPoint_t`, `geo::Vector_t` and
    ///             `geo::CryostatGeo::LocalVector_t`, and then use the method
    ///             `geo::CryostatGeo::toWorldCoords()` instead.
    template <typename Point>
    [[deprecated("use toWorldCoords() instead")]]
    Point LocalToWorld(Point const& local) const
      { return fTrans.LocalToWorld(local); }
    
    /// Transform point from local cryostat frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* cryo, double* world) const
      { fTrans.LocalToWorldVect(cryo, world); }
    
    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform point from world frame to local cryostat frame.
    void WorldToLocal(const double* world, double* cryo) const
      { fTrans.WorldToLocal(world, cryo); }
    
    /// Transform point from world frame to local cryostat frame.
    /// @deprecated This method breaks the distinction between local and global
    ///             vectors, since input and output vectors share the same type;
    ///             use the "official" vector types `geo::Point_t`,
    ///             `geo::CryostatGeo::LocalPoint_t`, `geo::Vector_t` and
    ///             `geo::CryostatGeo::LocalVector_t`, and then use the method
    ///             `geo::CryostatGeo::toLocalCoords()` instead.
    template <typename Point>
    [[deprecated("use toLocalCoords() instead")]]
    Point WorldToLocal(Point const& world) const
      { return fTrans.WorldToLocal(world); }
    
    /// Transform point from world frame to local cryostat frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* cryo) const
      { fTrans.WorldToLocalVect(world, cryo); }
    
    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// @}
    
    
    /// Method to sort TPCGeo objects
    void              SortSubVolumes(geo::GeoObjectSorter const& sorter);

    
    /// Performs all needed updates after geometry has sorted the cryostats
    void UpdateAfterSorting(geo::CryostatID cryoid);
    
  private:

    void FindTPC(std::vector<const TGeoNode*>& path,
		 unsigned int depth);
    void MakeTPC(std::vector<const TGeoNode*>& path,
		 int depth);
    
    void FindOpDet(std::vector<const TGeoNode*>& path,
		   unsigned int depth);
    void MakeOpDet(std::vector<const TGeoNode*>& path,
		   int depth);

    /// Fill the boundary information of the cryostat
    void InitCryoBoundaries();

  private:
    
    using LocalTransformation_t
      = geo::LocalTransformationGeo<TGeoHMatrix, LocalPoint_t, LocalVector_t>;
    
    LocalTransformation_t  fTrans;          ///< Cryostat-to-world transformation.
    std::vector<TPCGeo>    fTPCs;           ///< List of tpcs in this cryostat
    std::vector<OpDetGeo>  fOpDets;         ///< List of opdets in this cryostat
    TGeoVolume*            fVolume;         ///< Total volume of cryostat, called volCryostat in GDML file
    std::string            fOpDetGeoName;   ///< Name of opdet geometry elements in gdml
    geo::CryostatID        fID;             ///< ID of this cryostat
    
  };
}


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::CryostatGeo::PrintCryostatInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {
  
  //----------------------------------------------------------------------------
  out << "Cryostat " << std::string(ID());
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  out
    << " (" << Width() << " x " << Height() << " x " << Length() << ") cm^3 at "
      << GetCenter();
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  
  out << "\n" << indent
    << "hosts " << NTPC() << " TPCs (largest number of planes: " << MaxPlanes()
      << ", of wires: " << MaxWires() << ") and "
      << NOpDet() << " optical detectors"
    ;
  
  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  // print also the containing box
  geo::BoxBoundedGeo const& box = BoundingBox();
  out << "\n" << indent
    << "bounding box: " << box.Min() << " -- " << box.Max();
  
//  if (verbosity-- <= 0) return; // 3
  
  //----------------------------------------------------------------------------
} // geo::CryostatGeo::PrintCryostatInfo()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_CRYOSTATGEO_H
////////////////////////////////////////////////////////////////////////
