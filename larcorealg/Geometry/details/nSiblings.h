/**
 * The purpose of this file is to enable the separation of the
 * geometry iterators from the GeometryCore header.  The iterators are
 * defined before the GeometryCore class is defined, and because the
 * iterators require the GeometryCore object as part of the
 * implementation, there is a circularity that must be broken.
 *
 * To do this, we create free functions that call member functions of
 * the GeometryCore object in a separate compilation unit.  The
 * functions are resolved at link-time, avoiding the circularity.
 *
 * N.B. If it should be decided that the "N siblings" functionality
 *      should be templated (i.e. more than just GeometryCore might
 *      support the concept), then the circularity would be removed
 *      because GeometryCore woould become a dependent type.
 */

#ifndef LARCOREALG_GEOMETRY_DETAILS_NSIBLINGS_H
#define LARCOREALG_GEOMETRY_DETAILS_NSIBLINGS_H

// LArSoft libraries
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

namespace geo::details {
  unsigned int NSiblings(GeometryCore const* geom, CryostatID const& id);
  unsigned int NSiblings(GeometryCore const* geom, TPCID const& id);
  unsigned int NSiblings(GeometryCore const* geom, PlaneID const& id);
  unsigned int NSiblings(GeometryCore const* geom, WireID const& id);
  unsigned int NSiblings(GeometryCore const* geom, readout::TPCsetID const& id);
  unsigned int NSiblings(GeometryCore const* geom, readout::ROPID const& id);
}

#endif // LARCOREALG_GEOMETRY_DETAILS_NSIBLINGS_H
