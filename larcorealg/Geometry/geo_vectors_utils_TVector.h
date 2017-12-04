/**
 * @file   larcorealg/Geometry/geo_vectors_utils_TVector.h
 * @brief  Specializations of `geo_vectors_utils.h` for ROOT old vector types.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 2, 2017
 * @see    `larcorealg/Geometry/geo_vectors_utils.h`
 * 
 * This library provides facilities that can be used for both LArSoft geometry
 * vectors (`geo_vectors.h`) and ROOT `TVector3` and related, with the same
 * interface. 
 * 
 * This library depends on ROOT GenVector.
 * In the CET link list in `CMakeLists.txt`, link to `${ROOT_GENVECTOR}`.
 * 
 */

#ifndef LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_VECTORS_UTILS_TVECTOR_H
#define LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_VECTORS_UTILS_TVECTOR_H

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"

// C/C++ standard library
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"


namespace geo {
  namespace vect {
    
    /// @{
    /// @name `TVector3` conversions.
    
    //--------------------------------------------------------------------------
    /// Converts a vector into a `TVector3`.
    template <typename Vector>
    TVector3 toTVector3(Vector const& v) { return convertTo<TVector3>(v); }
    
    /// @}
    //--------------------------------------------------------------------------
    
    /// Utilities to print vector types.
    namespace dump {
      /// @{
      /// @name Output of old-style ROOT vectors (`TVector3` etc.).
      
      /// Print a `TVector2` to an output stream.
      template <typename Stream>
      Stream& operator<< (Stream&& out, TVector2 const& v) {
        out << "( " << v.X() << ", " << v.Y() << " )";
        return out;
      } // operator<< (TVector2)
      
      /// Print a `TVector3` to an output stream.
      template <typename Stream>
      Stream& operator<< (Stream&& out, TVector3 const& v) {
        out << "( " << v.X() << ", " << v.Y() << ", " << v.Z() << " )";
        return out;
      } // operator<< (TVector3)
      
      /// Print a `TLorentzVector` to an output stream.
      template <typename Stream>
      Stream& operator<< (Stream&& out, TLorentzVector const& v) {
        out 
          << "( " << v.X() << ", " << v.Y() << ", " << v.Z() << "; " << v.T() << " )";
        return out;
      } // operator<< (TLorentzVector)
      
      /// @}
      
    } // namespace dump
    //--------------------------------------------------------------------------
    
  } // namespace vect
} // namespace geo


// The only way some generic code has to see the operator<< is for them to be
// exposed in the same namespace as the vectors they dump are; in TVector case,
// that's the global namespace... (Boost unit test checks still fail though)
using namespace geo::vect::dump;


//------------------------------------------------------------------------------
//--- Specialisations
//---
namespace geo {
  namespace vect {
    
    //--------------------------------------------------------------------------
    // Specialisations for: TVector2
    template <>
    inline auto mag2<TVector2>(TVector2 const& v) { return v.Mod2(); }
    
    //--------------------------------------------------------------------------
    
  } // namespace vect
} // namespace geo
//------------------------------------------------------------------------------


#endif // LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_VECTORS_UTILS_TVECTOR_H
