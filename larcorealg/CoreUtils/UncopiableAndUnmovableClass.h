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

  // --- BEGIN -- Non-polymorphic classes --------------------------------------
  /**
   * @name Non-polymorphic classes
   * @anchor LArSoftCoreUtils_NonPolymorphicUncopiableUnmovable
   *
   * These base classes can be used to explicitly control movability and
   * copiability of the derived classes. It's a trade of a long class name
   * (and a long header file include statement) against five or six boilerplate
   * statements.
   *
   * The four classes in this group should be used only for classes with no
   * run-time polymorphism, for example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct UniqueData: protected lar::UncopiableClass {
   *
   *   std::vector<float> fLotsOfData{ 100'000'000U };
   *
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The class `UniqueData` in the example can be moved, but not copied.
   *
   * See
   * @ref LArSoftCoreUtils_PolymorphicUncopiableUnmovable "polymorphic classes"
   * for use with classes presenting run-time polymorphism.
   *
   */
  /// @{

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

  /// @}
  // --- END -- Non-polymorphic classes ----------------------------------------



  // --- BEGIN -- Run-time-polymorphic classes ---------------------------------
  /**
   * @name Run-time polymorphic classes
   * @anchor LArSoftCoreUtils_PolymorphicUncopiableUnmovable
   *
   * These base classes can be used to explicitly control movability and
   * copiability of the derived classes which feature run-time polymorphism.
   * It's a trade of a long class name (and a long header file include
   * statement) against five or six boilerplate statements.
   *
   * The four relevant classes in this group should be used only for classes
   * with run-time polymorphism, for example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct UniqueAlgorithmBase
   *   : protected lar::PolymorphicUncopiableAndUnmovableClass
   * {
   *
   *   virtual int compute() const = 0;
   *
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * All the classes derived from `UniqueAlgorithmBase` in the example can't be
   * copied nor moved (as `UniqueAlgorithmBase` as well, but that's moot
   * since it's pure virtual).
   * In general this may be a good thing to avoid slicing, that is an incomplete
   * copy produced by the copy or move constructor of a base class being used
   * instead of the actual underlying class.
   *
   * If run-time polymorphism is not needed, stick to
   * @ref LArSoftCoreUtils_NonPolymorphicUncopiableUnmovable "non-polymorphic classes"
   * instead.
   *
   */
  /// @{

  /**
   * @brief A simple polymorphic class, providing a virtual table.
   *
   * This class (and its derived ones) are going to be copiable and movable.
   */
  struct PolymorphicClass {

    PolymorphicClass() = default;

    virtual ~PolymorphicClass() = default;

  }; // PolymorphicClass


  /**
   * @brief A polymorphic empty class that can't be copied (moving is allowed).
   * @see   `PolymorphicClass`, `PolymorphicUnmovableClass`,
   *        `PolymorphicUncopiableAndUnmovableClass`
   *
   * A class derived from this one can still be copied with an explicit effort.
   * For example, to enable copy construction:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct CopiableClass: protected PolymorphicUncopiableClass {
   *   CopiableClass(CopiableClass const& from)
   *     : PolymorphicUncopiableClass() // , ...
   *     {
   *       // ...
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the default constructor of the base class can be called explicitly instead
   * of the copy constructor. To provide an assignment operation,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct CopyAssignableClass: protected PolymorphicUncopiableClass {
   *   CopyAssignableClass& operator= (CopyAssignableClass const& from)
   *     {
   *       // ...
   *       return *this;
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  struct PolymorphicUncopiableClass: virtual PolymorphicClass {

    PolymorphicUncopiableClass() = default;

    PolymorphicUncopiableClass(PolymorphicUncopiableClass const&) = delete;
    PolymorphicUncopiableClass(PolymorphicUncopiableClass&&) = default;

    PolymorphicUncopiableClass& operator= (PolymorphicUncopiableClass const&)
      = delete;
    PolymorphicUncopiableClass& operator= (PolymorphicUncopiableClass&&)
      = default;

  }; // PolymorphicUncopiableClass


  /**
   * @brief An empty polymorphic class that can't be moved (copy is allowed).
   * @see   `PolymorphicClass`, `PolymorphicUncopiableClass`,
   *        `PolymorphicUncopiableAndUnmovableClass`
   *
   * A class derived from this one can still be moved with an explicit effort.
   * For example, to enable move construction:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MoveableClass: protected PolymorphicUnmovableClass {
   *   MoveableClass(MoveableClass&& from)
   *     : PolymorphicUnmovableClass() // , ...
   *     {
   *       // ...
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the default constructor of the base class can be called explicitly instead
   * of the move constructor. To provide a move assignment operation,
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MoveAssignableClass: protected PolymorphicUnmovableClass {
   *   MoveAssignableClass& operator= (MoveAssignableClass&& from)
   *     {
   *       // ...
   *       return *this;
   *     }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  struct PolymorphicUnmovableClass: virtual PolymorphicClass {

    PolymorphicUnmovableClass() = default;

    PolymorphicUnmovableClass(PolymorphicUnmovableClass const&) = default;
    PolymorphicUnmovableClass(PolymorphicUnmovableClass&&) = delete;

    PolymorphicUnmovableClass& operator= (PolymorphicUnmovableClass const&)
      = default;
    PolymorphicUnmovableClass& operator= (PolymorphicUnmovableClass&&)
      = delete;

  }; // PolymorphicUnmovableClass




  /** **************************************************************************
   * @brief An empty class that can't be copied nor moved.
   * @see   `PolymorphicClass`, `PolymorphicUncopiableClass`,
   *        `PolymorphicUnmovableClass`
   *
   * A class derived from this one can still be copied and/or moved with an
   * explicit effort. See `PolymorphicUncopiableClass` and
   * `PolymorphicUnmovableClass` documentation for examples.
   */
  struct PolymorphicUncopiableAndUnmovableClass
    : PolymorphicUncopiableClass, PolymorphicUnmovableClass
  {};


  /// @}
  // --- END -- Run-time-polymorphic classes -----------------------------------

} // namespace lar

#endif // LARCORE_COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H
