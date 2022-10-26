#include "larcorealg/Geometry/details/nSiblings.h"
#include "larcorealg/Geometry/GeometryCore.h"

unsigned int geo::details::NSiblings(GeometryCore const* geom, CryostatID const& id)
{
  return geom->NSiblingElements(id);
}

unsigned int geo::details::NSiblings(GeometryCore const* geom, TPCID const& id)
{
  return geom->NSiblingElements(id);
}

unsigned int geo::details::NSiblings(GeometryCore const* geom, PlaneID const& id)
{
  return geom->NSiblingElements(id);
}

unsigned int geo::details::NSiblings(GeometryCore const* geom, WireID const& id)
{
  return geom->NSiblingElements(id);
}

unsigned int geo::details::NSiblings(GeometryCore const* geom, readout::TPCsetID const& id)
{
  return geom->NSiblingElements(id);
}

unsigned int geo::details::NSiblings(GeometryCore const* geom, readout::ROPID const& id)
{
  return geom->NSiblingElements(id);
}
