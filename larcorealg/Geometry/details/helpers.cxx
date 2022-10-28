#include "larcorealg/Geometry/details/helpers.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace geo::details {

  unsigned int NSiblings(GeometryCore const* geom, CryostatID const& id)
  {
    return geom->NSiblingElements(id);
  }

  unsigned int NSiblings(GeometryCore const* geom, TPCID const& id)
  {
    return geom->NSiblingElements(id);
  }

  unsigned int NSiblings(GeometryCore const* geom, PlaneID const& id)
  {
    return geom->NSiblingElements(id);
  }

  unsigned int NSiblings(GeometryCore const* geom, WireID const& id)
  {
    return geom->NSiblingElements(id);
  }

  unsigned int NSiblings(GeometryCore const* geom, readout::TPCsetID const& id)
  {
    return geom->NSiblingElements(id);
  }

  unsigned int NSiblings(GeometryCore const* geom, readout::ROPID const& id)
  {
    return geom->NSiblingElements(id);
  }

  CryostatGeo const* getElementPtr(GeometryCore const* geom, CryostatID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  TPCGeo const* getElementPtr(GeometryCore const* geom, TPCID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  PlaneGeo const* getElementPtr(GeometryCore const* geom, PlaneID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  WireGeo const* getElementPtr(GeometryCore const* geom, WireID const& id)
  {
    assert(geom);
    return geom->GetElementPtr(id);
  }

  bool validElement(GeometryCore const* geom, CryostatID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, TPCID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, PlaneID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

  bool validElement(GeometryCore const* geom, WireID const& id)
  {
    return geom && getElementPtr(geom, id) != nullptr;
  }

}
