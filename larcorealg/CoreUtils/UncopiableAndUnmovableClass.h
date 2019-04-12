/**
 * @file   larcorealg/CoreUtils/UncopiableAndUnmovableClass.h
 * @brief  Defines classes that can't be copied nor moved.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * This library is currently a pure header.
 *
 */

#ifndef LARCORE_COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H
#define LARCORE_COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H


namespace lar {

  /** **************************************************************************
   * @brief An empty class that can't be copied (moving is allowed).
   * @see   `UnmovableClass`, `UncopiableAndUnmovableClass`
   *
   * A class derived from this one can still be copied with an explicit effort.
   * For example, to enable copy construction:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct CopiableClass: protected UncopiableAndUnmovableClass {
   *   CopiableClass(CopiableClass const& from)
   *     : UncopiableClass() // , ...
   *     {
   *       // ...
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the default constructor of the base class can be called explicitly instead
   * of the copy constructor. To provide an assignment operation,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct CopyAssignableClass: protected UncopiableClass {
   *   CopyAssignableClass& operator= (CopyAssignableClass const& from)
   *     {
   *       // ...
   *       return *this;
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  struct UncopiableClass {

    /// Default constructor
    UncopiableClass() = default;

    // @{
    /// Deleted copy and move constructors and assignments
    UncopiableClass(UncopiableClass const&) = delete;
    UncopiableClass(UncopiableClass&&) = default;

    UncopiableClass& operator= (UncopiableClass const&) = delete;
    UncopiableClass& operator= (UncopiableClass&&) = default;
    // @}

    /// Default destructor
    ~UncopiableClass() = default;

  }; // UncopiableClass


  /** **************************************************************************
   * @brief An empty class that can't be moved (copy is allowed).
   * @see   `UncopiableClass`, `UncopiableAndUnmovableClass`
   *
   * A class derived from this one can still be moved with an explicit effort.
   * For example, to enable move construction:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MoveableClass: protected UnmovableClass {
   *   MoveableClass(MoveableClass&& from)
   *     : UnmovableClass() // , ...
   *     {
   *       // ...
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the default constructor of the base class can be called explicitly instead
   * of the move constructor. To provide a move assignment operation,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MoveAssignableClass: protected UnmovableClass {
   *   MoveAssignableClass& operator= (MoveAssignableClass&& from)
   *     {
   *       // ...
   *       return *this;
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  struct UnmovableClass {

    /// Default constructor.
    UnmovableClass() = default;

    //@{
    /// Default copy constructor and assignment.
    UnmovableClass(UnmovableClass const&) = default;
    UnmovableClass& operator= (UnmovableClass const&) = default;
    //@}

    //@{
    /// Deleted move constructor and assignment.
    UnmovableClass(UnmovableClass&&) = delete;
    UnmovableClass& operator= (UnmovableClass&&) = delete;
    //@}

    /// Default destructor.
    ~UnmovableClass() = default;

  }; // UnmovableClass


  /** **************************************************************************
   * @brief An empty class that can't be copied nor moved.
   * @see `UncopiableClass`, `UnmovableClass`
   *
   * A class derived from this one can still be copied and/or moved with an
   * explicit effort. See `UncopiableClass` and `UnmovableClass` for examples.
   */
  struct UncopiableAndUnmovableClass
    : public UncopiableClass, public UnmovableClass
    {};


} // namespace lar

#endif // LARCORE_COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H
