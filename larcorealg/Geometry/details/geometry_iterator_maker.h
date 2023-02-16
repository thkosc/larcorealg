#ifndef LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATOR_MAKER_H
#define LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATOR_MAKER_H

#include "larcorealg/CoreUtils/span.h"
#include "larcorealg/Geometry/details/geometry_iterators.h"

#include <type_traits>

namespace geo::details {
  template <typename BaseID, typename GeoID>
  static constexpr bool is_base_of_strict{std::is_base_of<BaseID, GeoID>{} &&
                                          !std::is_same<BaseID, GeoID>{}};

  template <typename T, typename = void>
  struct IteratorMaker {
    using ID_t = typename T::ID_t;
    using begin_type = element_iterator_for<T>;
    using end_type = element_iterator_for<T>;
    using range_type = util::span<begin_type, end_type>;

    template <typename Geom>
    static begin_type create_begin(Geom const* geom)
    {
      return {geom, id_iterator<ID_t>{geom, geom->template GetBeginID<ID_t>()}};
    }

    template <typename Geom>
    static end_type create_end(Geom const* geom)
    {
      return {geom, id_iterator<ID_t>{geom, geom->template GetEndID<ID_t>()}};
    }

    template <typename Geom>
    static range_type create_range(Geom const* geom)
    {
      return {create_begin(geom), create_end(geom)};
    }

    template <typename Geom, typename BaseID>
    static begin_type create_begin(Geom const* geom, BaseID const& id)
    {
      static_assert(is_base_of_strict<BaseID, ID_t>);
      return {geom, id_iterator<ID_t>{geom, geom->template GetBeginID<ID_t>(id)}};
    }

    template <typename Geom, typename BaseID>
    static end_type create_end(Geom const* geom, BaseID const& id)
    {
      static_assert(is_base_of_strict<BaseID, ID_t>);
      return {geom, id_iterator<ID_t>{geom, geom->template GetEndID<ID_t>(id)}};
    }

    template <typename Geom, typename BaseID>
    static range_type create_range(Geom const* geom, BaseID const& id)
    {
      return {create_begin(geom, id), create_end(geom, id)};
    }
  };

  template <typename T>
  struct IteratorMaker<T, std::enable_if_t<std::is_base_of_v<CryostatID, T>>> {
    using begin_type = id_iterator<T>;
    using end_type = id_iterator<T>;
    using range_type = util::span<begin_type, end_type>;

    template <typename Geom>
    static begin_type create_begin(Geom const* geom)
    {
      return {geom, geom->template GetBeginID<T>()};
    }

    template <typename Geom>
    static end_type create_end(Geom const* geom)
    {
      return {geom, geom->template GetEndID<T>()};
    }

    template <typename Geom>
    static range_type create_range(Geom const* geom)
    {
      return {create_begin(geom), create_end(geom)};
    }

    template <typename Geom, typename BaseID>
    static begin_type create_begin(Geom const* geom, BaseID const& id)
    {
      static_assert(is_base_of_strict<BaseID, T>);
      return {geom, geom->template GetBeginID<T>(id)};
    }

    template <typename Geom, typename BaseID>
    static end_type create_end(Geom const* geom, BaseID const& id)
    {
      static_assert(is_base_of_strict<BaseID, T>);
      return {geom, geom->template GetEndID<T>(id)};
    }

    template <typename Geom, typename BaseID>
    static range_type create_range(Geom const* geom, BaseID const& id)
    {
      return {create_begin(geom, id), create_end(geom, id)};
    }
  };

  template <typename T>
  using begin_type = typename IteratorMaker<T>::begin_type;

  template <typename T>
  using end_type = typename IteratorMaker<T>::end_type;

  template <typename T>
  using range_type = typename IteratorMaker<T>::range_type;
}

#endif /* LARCOREALG_GEOMETRY_DETAILS_GEOMETRY_ITERATOR_MAKER_H */
