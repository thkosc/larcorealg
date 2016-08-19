// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME larcorealg_GeoAlgoDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "GeoAABox.h"
#include "GeoAlgoConstants.h"
#include "GeoAlgoException.h"
#include "GeoAlgo.h"
#include "GeoAlgo-TypeDef.h"
#include "GeoCone.h"
#include "GeoCylinder.h"
#include "GeoDirectedLine.h"
#include "GeoHalfLine.h"
#include "GeoLine.h"
#include "GeoLineSegment.h"
#include "GeoObjCollection.h"
#include "GeoSphere.h"
#include "GeoTrajectory.h"
#include "GeoVector.h"

// Header files passed via #pragma extra_include

namespace geoalgo {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *geoalgo_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("geoalgo", 0 /*version*/, "GeoAlgoConstants.h", 6,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &geoalgo_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *geoalgo_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *geoalgocLcLVector_Dictionary();
   static void geoalgocLcLVector_TClassManip(TClass*);
   static void *new_geoalgocLcLVector(void *p = 0);
   static void *newArray_geoalgocLcLVector(Long_t size, void *p);
   static void delete_geoalgocLcLVector(void *p);
   static void deleteArray_geoalgocLcLVector(void *p);
   static void destruct_geoalgocLcLVector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Vector*)
   {
      ::geoalgo::Vector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Vector));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Vector", "GeoVector.h", 35,
                  typeid(::geoalgo::Vector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLVector_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Vector) );
      instance.SetNew(&new_geoalgocLcLVector);
      instance.SetNewArray(&newArray_geoalgocLcLVector);
      instance.SetDelete(&delete_geoalgocLcLVector);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLVector);
      instance.SetDestructor(&destruct_geoalgocLcLVector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Vector*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Vector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Vector*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLVector_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Vector*)0x0)->GetClass();
      geoalgocLcLVector_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLVector_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLTrajectory_Dictionary();
   static void geoalgocLcLTrajectory_TClassManip(TClass*);
   static void *new_geoalgocLcLTrajectory(void *p = 0);
   static void *newArray_geoalgocLcLTrajectory(Long_t size, void *p);
   static void delete_geoalgocLcLTrajectory(void *p);
   static void deleteArray_geoalgocLcLTrajectory(void *p);
   static void destruct_geoalgocLcLTrajectory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Trajectory*)
   {
      ::geoalgo::Trajectory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Trajectory));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Trajectory", "GeoTrajectory.h", 27,
                  typeid(::geoalgo::Trajectory), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLTrajectory_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Trajectory) );
      instance.SetNew(&new_geoalgocLcLTrajectory);
      instance.SetNewArray(&newArray_geoalgocLcLTrajectory);
      instance.SetDelete(&delete_geoalgocLcLTrajectory);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLTrajectory);
      instance.SetDestructor(&destruct_geoalgocLcLTrajectory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Trajectory*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Trajectory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Trajectory*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLTrajectory_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Trajectory*)0x0)->GetClass();
      geoalgocLcLTrajectory_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLTrajectory_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLHalfLine_Dictionary();
   static void geoalgocLcLHalfLine_TClassManip(TClass*);
   static void *new_geoalgocLcLHalfLine(void *p = 0);
   static void *newArray_geoalgocLcLHalfLine(Long_t size, void *p);
   static void delete_geoalgocLcLHalfLine(void *p);
   static void deleteArray_geoalgocLcLHalfLine(void *p);
   static void destruct_geoalgocLcLHalfLine(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::HalfLine*)
   {
      ::geoalgo::HalfLine *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::HalfLine));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::HalfLine", "GeoHalfLine.h", 26,
                  typeid(::geoalgo::HalfLine), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLHalfLine_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::HalfLine) );
      instance.SetNew(&new_geoalgocLcLHalfLine);
      instance.SetNewArray(&newArray_geoalgocLcLHalfLine);
      instance.SetDelete(&delete_geoalgocLcLHalfLine);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLHalfLine);
      instance.SetDestructor(&destruct_geoalgocLcLHalfLine);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::HalfLine*)
   {
      return GenerateInitInstanceLocal((::geoalgo::HalfLine*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::HalfLine*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLHalfLine_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::HalfLine*)0x0)->GetClass();
      geoalgocLcLHalfLine_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLHalfLine_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLLine_Dictionary();
   static void geoalgocLcLLine_TClassManip(TClass*);
   static void *new_geoalgocLcLLine(void *p = 0);
   static void *newArray_geoalgocLcLLine(Long_t size, void *p);
   static void delete_geoalgocLcLLine(void *p);
   static void deleteArray_geoalgocLcLLine(void *p);
   static void destruct_geoalgocLcLLine(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Line*)
   {
      ::geoalgo::Line *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Line));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Line", "GeoLine.h", 27,
                  typeid(::geoalgo::Line), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLLine_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Line) );
      instance.SetNew(&new_geoalgocLcLLine);
      instance.SetNewArray(&newArray_geoalgocLcLLine);
      instance.SetDelete(&delete_geoalgocLcLLine);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLLine);
      instance.SetDestructor(&destruct_geoalgocLcLLine);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Line*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Line*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Line*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLLine_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Line*)0x0)->GetClass();
      geoalgocLcLLine_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLLine_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLDirectedLine_Dictionary();
   static void geoalgocLcLDirectedLine_TClassManip(TClass*);
   static void *new_geoalgocLcLDirectedLine(void *p = 0);
   static void *newArray_geoalgocLcLDirectedLine(Long_t size, void *p);
   static void delete_geoalgocLcLDirectedLine(void *p);
   static void deleteArray_geoalgocLcLDirectedLine(void *p);
   static void destruct_geoalgocLcLDirectedLine(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::DirectedLine*)
   {
      ::geoalgo::DirectedLine *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::DirectedLine));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::DirectedLine", "GeoDirectedLine.h", 28,
                  typeid(::geoalgo::DirectedLine), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLDirectedLine_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::DirectedLine) );
      instance.SetNew(&new_geoalgocLcLDirectedLine);
      instance.SetNewArray(&newArray_geoalgocLcLDirectedLine);
      instance.SetDelete(&delete_geoalgocLcLDirectedLine);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLDirectedLine);
      instance.SetDestructor(&destruct_geoalgocLcLDirectedLine);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::DirectedLine*)
   {
      return GenerateInitInstanceLocal((::geoalgo::DirectedLine*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::DirectedLine*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLDirectedLine_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::DirectedLine*)0x0)->GetClass();
      geoalgocLcLDirectedLine_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLDirectedLine_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLLineSegment_Dictionary();
   static void geoalgocLcLLineSegment_TClassManip(TClass*);
   static void *new_geoalgocLcLLineSegment(void *p = 0);
   static void *newArray_geoalgocLcLLineSegment(Long_t size, void *p);
   static void delete_geoalgocLcLLineSegment(void *p);
   static void deleteArray_geoalgocLcLLineSegment(void *p);
   static void destruct_geoalgocLcLLineSegment(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::LineSegment*)
   {
      ::geoalgo::LineSegment *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::LineSegment));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::LineSegment", "GeoLineSegment.h", 25,
                  typeid(::geoalgo::LineSegment), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLLineSegment_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::LineSegment) );
      instance.SetNew(&new_geoalgocLcLLineSegment);
      instance.SetNewArray(&newArray_geoalgocLcLLineSegment);
      instance.SetDelete(&delete_geoalgocLcLLineSegment);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLLineSegment);
      instance.SetDestructor(&destruct_geoalgocLcLLineSegment);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::LineSegment*)
   {
      return GenerateInitInstanceLocal((::geoalgo::LineSegment*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::LineSegment*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLLineSegment_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::LineSegment*)0x0)->GetClass();
      geoalgocLcLLineSegment_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLLineSegment_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLAABox_Dictionary();
   static void geoalgocLcLAABox_TClassManip(TClass*);
   static void *new_geoalgocLcLAABox(void *p = 0);
   static void *newArray_geoalgocLcLAABox(Long_t size, void *p);
   static void delete_geoalgocLcLAABox(void *p);
   static void deleteArray_geoalgocLcLAABox(void *p);
   static void destruct_geoalgocLcLAABox(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::AABox*)
   {
      ::geoalgo::AABox *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::AABox));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::AABox", "GeoAABox.h", 34,
                  typeid(::geoalgo::AABox), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLAABox_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::AABox) );
      instance.SetNew(&new_geoalgocLcLAABox);
      instance.SetNewArray(&newArray_geoalgocLcLAABox);
      instance.SetDelete(&delete_geoalgocLcLAABox);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLAABox);
      instance.SetDestructor(&destruct_geoalgocLcLAABox);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::AABox*)
   {
      return GenerateInitInstanceLocal((::geoalgo::AABox*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::AABox*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLAABox_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::AABox*)0x0)->GetClass();
      geoalgocLcLAABox_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLAABox_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLCylinder_Dictionary();
   static void geoalgocLcLCylinder_TClassManip(TClass*);
   static void *new_geoalgocLcLCylinder(void *p = 0);
   static void *newArray_geoalgocLcLCylinder(Long_t size, void *p);
   static void delete_geoalgocLcLCylinder(void *p);
   static void deleteArray_geoalgocLcLCylinder(void *p);
   static void destruct_geoalgocLcLCylinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Cylinder*)
   {
      ::geoalgo::Cylinder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Cylinder));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Cylinder", "GeoCylinder.h", 29,
                  typeid(::geoalgo::Cylinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLCylinder_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Cylinder) );
      instance.SetNew(&new_geoalgocLcLCylinder);
      instance.SetNewArray(&newArray_geoalgocLcLCylinder);
      instance.SetDelete(&delete_geoalgocLcLCylinder);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLCylinder);
      instance.SetDestructor(&destruct_geoalgocLcLCylinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Cylinder*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Cylinder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Cylinder*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLCylinder_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Cylinder*)0x0)->GetClass();
      geoalgocLcLCylinder_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLCylinder_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLCone_Dictionary();
   static void geoalgocLcLCone_TClassManip(TClass*);
   static void *new_geoalgocLcLCone(void *p = 0);
   static void *newArray_geoalgocLcLCone(Long_t size, void *p);
   static void delete_geoalgocLcLCone(void *p);
   static void deleteArray_geoalgocLcLCone(void *p);
   static void destruct_geoalgocLcLCone(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Cone*)
   {
      ::geoalgo::Cone *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Cone));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Cone", "GeoCone.h", 27,
                  typeid(::geoalgo::Cone), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLCone_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Cone) );
      instance.SetNew(&new_geoalgocLcLCone);
      instance.SetNewArray(&newArray_geoalgocLcLCone);
      instance.SetDelete(&delete_geoalgocLcLCone);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLCone);
      instance.SetDestructor(&destruct_geoalgocLcLCone);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Cone*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Cone*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Cone*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLCone_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Cone*)0x0)->GetClass();
      geoalgocLcLCone_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLCone_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLSphere_Dictionary();
   static void geoalgocLcLSphere_TClassManip(TClass*);
   static void *new_geoalgocLcLSphere(void *p = 0);
   static void *newArray_geoalgocLcLSphere(Long_t size, void *p);
   static void delete_geoalgocLcLSphere(void *p);
   static void deleteArray_geoalgocLcLSphere(void *p);
   static void destruct_geoalgocLcLSphere(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::Sphere*)
   {
      ::geoalgo::Sphere *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::Sphere));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::Sphere", "GeoSphere.h", 25,
                  typeid(::geoalgo::Sphere), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLSphere_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::Sphere) );
      instance.SetNew(&new_geoalgocLcLSphere);
      instance.SetNewArray(&newArray_geoalgocLcLSphere);
      instance.SetDelete(&delete_geoalgocLcLSphere);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLSphere);
      instance.SetDestructor(&destruct_geoalgocLcLSphere);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::Sphere*)
   {
      return GenerateInitInstanceLocal((::geoalgo::Sphere*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::Sphere*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLSphere_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::Sphere*)0x0)->GetClass();
      geoalgocLcLSphere_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLSphere_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairlEgeoalgocLcLVectorcOstringgR_Dictionary();
   static void pairlEgeoalgocLcLVectorcOstringgR_TClassManip(TClass*);
   static void *new_pairlEgeoalgocLcLVectorcOstringgR(void *p = 0);
   static void *newArray_pairlEgeoalgocLcLVectorcOstringgR(Long_t size, void *p);
   static void delete_pairlEgeoalgocLcLVectorcOstringgR(void *p);
   static void deleteArray_pairlEgeoalgocLcLVectorcOstringgR(void *p);
   static void destruct_pairlEgeoalgocLcLVectorcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<geoalgo::Vector,string>*)
   {
      pair<geoalgo::Vector,string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<geoalgo::Vector,string>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<geoalgo::Vector,string>", "string", 96,
                  typeid(pair<geoalgo::Vector,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEgeoalgocLcLVectorcOstringgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<geoalgo::Vector,string>) );
      instance.SetNew(&new_pairlEgeoalgocLcLVectorcOstringgR);
      instance.SetNewArray(&newArray_pairlEgeoalgocLcLVectorcOstringgR);
      instance.SetDelete(&delete_pairlEgeoalgocLcLVectorcOstringgR);
      instance.SetDeleteArray(&deleteArray_pairlEgeoalgocLcLVectorcOstringgR);
      instance.SetDestructor(&destruct_pairlEgeoalgocLcLVectorcOstringgR);
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const pair<geoalgo::Vector,string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEgeoalgocLcLVectorcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<geoalgo::Vector,string>*)0x0)->GetClass();
      pairlEgeoalgocLcLVectorcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEgeoalgocLcLVectorcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLGeoAlgo_Dictionary();
   static void geoalgocLcLGeoAlgo_TClassManip(TClass*);
   static void *new_geoalgocLcLGeoAlgo(void *p = 0);
   static void *newArray_geoalgocLcLGeoAlgo(Long_t size, void *p);
   static void delete_geoalgocLcLGeoAlgo(void *p);
   static void deleteArray_geoalgocLcLGeoAlgo(void *p);
   static void destruct_geoalgocLcLGeoAlgo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::GeoAlgo*)
   {
      ::geoalgo::GeoAlgo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::GeoAlgo));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::GeoAlgo", "GeoAlgo.h", 42,
                  typeid(::geoalgo::GeoAlgo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLGeoAlgo_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::GeoAlgo) );
      instance.SetNew(&new_geoalgocLcLGeoAlgo);
      instance.SetNewArray(&newArray_geoalgocLcLGeoAlgo);
      instance.SetDelete(&delete_geoalgocLcLGeoAlgo);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLGeoAlgo);
      instance.SetDestructor(&destruct_geoalgocLcLGeoAlgo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::GeoAlgo*)
   {
      return GenerateInitInstanceLocal((::geoalgo::GeoAlgo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::GeoAlgo*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLGeoAlgo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::GeoAlgo*)0x0)->GetClass();
      geoalgocLcLGeoAlgo_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLGeoAlgo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *geoalgocLcLGeoObjCollection_Dictionary();
   static void geoalgocLcLGeoObjCollection_TClassManip(TClass*);
   static void *new_geoalgocLcLGeoObjCollection(void *p = 0);
   static void *newArray_geoalgocLcLGeoObjCollection(Long_t size, void *p);
   static void delete_geoalgocLcLGeoObjCollection(void *p);
   static void deleteArray_geoalgocLcLGeoObjCollection(void *p);
   static void destruct_geoalgocLcLGeoObjCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::geoalgo::GeoObjCollection*)
   {
      ::geoalgo::GeoObjCollection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::geoalgo::GeoObjCollection));
      static ::ROOT::TGenericClassInfo 
         instance("geoalgo::GeoObjCollection", "GeoObjCollection.h", 31,
                  typeid(::geoalgo::GeoObjCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &geoalgocLcLGeoObjCollection_Dictionary, isa_proxy, 4,
                  sizeof(::geoalgo::GeoObjCollection) );
      instance.SetNew(&new_geoalgocLcLGeoObjCollection);
      instance.SetNewArray(&newArray_geoalgocLcLGeoObjCollection);
      instance.SetDelete(&delete_geoalgocLcLGeoObjCollection);
      instance.SetDeleteArray(&deleteArray_geoalgocLcLGeoObjCollection);
      instance.SetDestructor(&destruct_geoalgocLcLGeoObjCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::geoalgo::GeoObjCollection*)
   {
      return GenerateInitInstanceLocal((::geoalgo::GeoObjCollection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::geoalgo::GeoObjCollection*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *geoalgocLcLGeoObjCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::geoalgo::GeoObjCollection*)0x0)->GetClass();
      geoalgocLcLGeoObjCollection_TClassManip(theClass);
   return theClass;
   }

   static void geoalgocLcLGeoObjCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLVector(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Vector : new ::geoalgo::Vector;
   }
   static void *newArray_geoalgocLcLVector(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Vector[nElements] : new ::geoalgo::Vector[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLVector(void *p) {
      delete ((::geoalgo::Vector*)p);
   }
   static void deleteArray_geoalgocLcLVector(void *p) {
      delete [] ((::geoalgo::Vector*)p);
   }
   static void destruct_geoalgocLcLVector(void *p) {
      typedef ::geoalgo::Vector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Vector

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLTrajectory(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Trajectory : new ::geoalgo::Trajectory;
   }
   static void *newArray_geoalgocLcLTrajectory(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Trajectory[nElements] : new ::geoalgo::Trajectory[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLTrajectory(void *p) {
      delete ((::geoalgo::Trajectory*)p);
   }
   static void deleteArray_geoalgocLcLTrajectory(void *p) {
      delete [] ((::geoalgo::Trajectory*)p);
   }
   static void destruct_geoalgocLcLTrajectory(void *p) {
      typedef ::geoalgo::Trajectory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Trajectory

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLHalfLine(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::HalfLine : new ::geoalgo::HalfLine;
   }
   static void *newArray_geoalgocLcLHalfLine(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::HalfLine[nElements] : new ::geoalgo::HalfLine[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLHalfLine(void *p) {
      delete ((::geoalgo::HalfLine*)p);
   }
   static void deleteArray_geoalgocLcLHalfLine(void *p) {
      delete [] ((::geoalgo::HalfLine*)p);
   }
   static void destruct_geoalgocLcLHalfLine(void *p) {
      typedef ::geoalgo::HalfLine current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::HalfLine

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLLine(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Line : new ::geoalgo::Line;
   }
   static void *newArray_geoalgocLcLLine(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Line[nElements] : new ::geoalgo::Line[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLLine(void *p) {
      delete ((::geoalgo::Line*)p);
   }
   static void deleteArray_geoalgocLcLLine(void *p) {
      delete [] ((::geoalgo::Line*)p);
   }
   static void destruct_geoalgocLcLLine(void *p) {
      typedef ::geoalgo::Line current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Line

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLDirectedLine(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::DirectedLine : new ::geoalgo::DirectedLine;
   }
   static void *newArray_geoalgocLcLDirectedLine(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::DirectedLine[nElements] : new ::geoalgo::DirectedLine[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLDirectedLine(void *p) {
      delete ((::geoalgo::DirectedLine*)p);
   }
   static void deleteArray_geoalgocLcLDirectedLine(void *p) {
      delete [] ((::geoalgo::DirectedLine*)p);
   }
   static void destruct_geoalgocLcLDirectedLine(void *p) {
      typedef ::geoalgo::DirectedLine current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::DirectedLine

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLLineSegment(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::LineSegment : new ::geoalgo::LineSegment;
   }
   static void *newArray_geoalgocLcLLineSegment(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::LineSegment[nElements] : new ::geoalgo::LineSegment[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLLineSegment(void *p) {
      delete ((::geoalgo::LineSegment*)p);
   }
   static void deleteArray_geoalgocLcLLineSegment(void *p) {
      delete [] ((::geoalgo::LineSegment*)p);
   }
   static void destruct_geoalgocLcLLineSegment(void *p) {
      typedef ::geoalgo::LineSegment current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::LineSegment

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLAABox(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::AABox : new ::geoalgo::AABox;
   }
   static void *newArray_geoalgocLcLAABox(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::AABox[nElements] : new ::geoalgo::AABox[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLAABox(void *p) {
      delete ((::geoalgo::AABox*)p);
   }
   static void deleteArray_geoalgocLcLAABox(void *p) {
      delete [] ((::geoalgo::AABox*)p);
   }
   static void destruct_geoalgocLcLAABox(void *p) {
      typedef ::geoalgo::AABox current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::AABox

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLCylinder(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Cylinder : new ::geoalgo::Cylinder;
   }
   static void *newArray_geoalgocLcLCylinder(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Cylinder[nElements] : new ::geoalgo::Cylinder[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLCylinder(void *p) {
      delete ((::geoalgo::Cylinder*)p);
   }
   static void deleteArray_geoalgocLcLCylinder(void *p) {
      delete [] ((::geoalgo::Cylinder*)p);
   }
   static void destruct_geoalgocLcLCylinder(void *p) {
      typedef ::geoalgo::Cylinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Cylinder

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLCone(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Cone : new ::geoalgo::Cone;
   }
   static void *newArray_geoalgocLcLCone(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Cone[nElements] : new ::geoalgo::Cone[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLCone(void *p) {
      delete ((::geoalgo::Cone*)p);
   }
   static void deleteArray_geoalgocLcLCone(void *p) {
      delete [] ((::geoalgo::Cone*)p);
   }
   static void destruct_geoalgocLcLCone(void *p) {
      typedef ::geoalgo::Cone current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Cone

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLSphere(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Sphere : new ::geoalgo::Sphere;
   }
   static void *newArray_geoalgocLcLSphere(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::Sphere[nElements] : new ::geoalgo::Sphere[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLSphere(void *p) {
      delete ((::geoalgo::Sphere*)p);
   }
   static void deleteArray_geoalgocLcLSphere(void *p) {
      delete [] ((::geoalgo::Sphere*)p);
   }
   static void destruct_geoalgocLcLSphere(void *p) {
      typedef ::geoalgo::Sphere current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::Sphere

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEgeoalgocLcLVectorcOstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<geoalgo::Vector,string> : new pair<geoalgo::Vector,string>;
   }
   static void *newArray_pairlEgeoalgocLcLVectorcOstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<geoalgo::Vector,string>[nElements] : new pair<geoalgo::Vector,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEgeoalgocLcLVectorcOstringgR(void *p) {
      delete ((pair<geoalgo::Vector,string>*)p);
   }
   static void deleteArray_pairlEgeoalgocLcLVectorcOstringgR(void *p) {
      delete [] ((pair<geoalgo::Vector,string>*)p);
   }
   static void destruct_pairlEgeoalgocLcLVectorcOstringgR(void *p) {
      typedef pair<geoalgo::Vector,string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<geoalgo::Vector,string>

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLGeoAlgo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::GeoAlgo : new ::geoalgo::GeoAlgo;
   }
   static void *newArray_geoalgocLcLGeoAlgo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::GeoAlgo[nElements] : new ::geoalgo::GeoAlgo[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLGeoAlgo(void *p) {
      delete ((::geoalgo::GeoAlgo*)p);
   }
   static void deleteArray_geoalgocLcLGeoAlgo(void *p) {
      delete [] ((::geoalgo::GeoAlgo*)p);
   }
   static void destruct_geoalgocLcLGeoAlgo(void *p) {
      typedef ::geoalgo::GeoAlgo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::GeoAlgo

namespace ROOT {
   // Wrappers around operator new
   static void *new_geoalgocLcLGeoObjCollection(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::GeoObjCollection : new ::geoalgo::GeoObjCollection;
   }
   static void *newArray_geoalgocLcLGeoObjCollection(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::geoalgo::GeoObjCollection[nElements] : new ::geoalgo::GeoObjCollection[nElements];
   }
   // Wrapper around operator delete
   static void delete_geoalgocLcLGeoObjCollection(void *p) {
      delete ((::geoalgo::GeoObjCollection*)p);
   }
   static void deleteArray_geoalgocLcLGeoObjCollection(void *p) {
      delete [] ((::geoalgo::GeoObjCollection*)p);
   }
   static void destruct_geoalgocLcLGeoObjCollection(void *p) {
      typedef ::geoalgo::GeoObjCollection current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::geoalgo::GeoObjCollection

namespace ROOT {
   static TClass *vectorlEvectorlEgeoalgocLcLVectorgRsPgR_Dictionary();
   static void vectorlEvectorlEgeoalgocLcLVectorgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p);
   static void destruct_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<geoalgo::Vector> >*)
   {
      vector<vector<geoalgo::Vector> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<geoalgo::Vector> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<geoalgo::Vector> >", -2, "vector", 214,
                  typeid(vector<vector<geoalgo::Vector> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEgeoalgocLcLVectorgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<geoalgo::Vector> >) );
      instance.SetNew(&new_vectorlEvectorlEgeoalgocLcLVectorgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEgeoalgocLcLVectorgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEgeoalgocLcLVectorgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<geoalgo::Vector> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<geoalgo::Vector> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEgeoalgocLcLVectorgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<geoalgo::Vector> >*)0x0)->GetClass();
      vectorlEvectorlEgeoalgocLcLVectorgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEgeoalgocLcLVectorgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<geoalgo::Vector> > : new vector<vector<geoalgo::Vector> >;
   }
   static void *newArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<geoalgo::Vector> >[nElements] : new vector<vector<geoalgo::Vector> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p) {
      delete ((vector<vector<geoalgo::Vector> >*)p);
   }
   static void deleteArray_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p) {
      delete [] ((vector<vector<geoalgo::Vector> >*)p);
   }
   static void destruct_vectorlEvectorlEgeoalgocLcLVectorgRsPgR(void *p) {
      typedef vector<vector<geoalgo::Vector> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<geoalgo::Vector> >

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 214,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLVectorgR_Dictionary();
   static void vectorlEgeoalgocLcLVectorgR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLVectorgR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLVectorgR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLVectorgR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLVectorgR(void *p);
   static void destruct_vectorlEgeoalgocLcLVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Vector>*)
   {
      vector<geoalgo::Vector> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Vector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Vector>", -2, "vector", 214,
                  typeid(vector<geoalgo::Vector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLVectorgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Vector>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLVectorgR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLVectorgR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLVectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLVectorgR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Vector> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Vector>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Vector>*)0x0)->GetClass();
      vectorlEgeoalgocLcLVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLVectorgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Vector> : new vector<geoalgo::Vector>;
   }
   static void *newArray_vectorlEgeoalgocLcLVectorgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Vector>[nElements] : new vector<geoalgo::Vector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLVectorgR(void *p) {
      delete ((vector<geoalgo::Vector>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLVectorgR(void *p) {
      delete [] ((vector<geoalgo::Vector>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLVectorgR(void *p) {
      typedef vector<geoalgo::Vector> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Vector>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLTrajectorygR_Dictionary();
   static void vectorlEgeoalgocLcLTrajectorygR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLTrajectorygR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLTrajectorygR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLTrajectorygR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLTrajectorygR(void *p);
   static void destruct_vectorlEgeoalgocLcLTrajectorygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Trajectory>*)
   {
      vector<geoalgo::Trajectory> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Trajectory>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Trajectory>", -2, "vector", 214,
                  typeid(vector<geoalgo::Trajectory>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLTrajectorygR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Trajectory>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLTrajectorygR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLTrajectorygR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLTrajectorygR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLTrajectorygR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLTrajectorygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Trajectory> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Trajectory>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLTrajectorygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Trajectory>*)0x0)->GetClass();
      vectorlEgeoalgocLcLTrajectorygR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLTrajectorygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLTrajectorygR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Trajectory> : new vector<geoalgo::Trajectory>;
   }
   static void *newArray_vectorlEgeoalgocLcLTrajectorygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Trajectory>[nElements] : new vector<geoalgo::Trajectory>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLTrajectorygR(void *p) {
      delete ((vector<geoalgo::Trajectory>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLTrajectorygR(void *p) {
      delete [] ((vector<geoalgo::Trajectory>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLTrajectorygR(void *p) {
      typedef vector<geoalgo::Trajectory> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Trajectory>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLSpheregR_Dictionary();
   static void vectorlEgeoalgocLcLSpheregR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLSpheregR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLSpheregR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLSpheregR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLSpheregR(void *p);
   static void destruct_vectorlEgeoalgocLcLSpheregR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Sphere>*)
   {
      vector<geoalgo::Sphere> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Sphere>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Sphere>", -2, "vector", 214,
                  typeid(vector<geoalgo::Sphere>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLSpheregR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Sphere>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLSpheregR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLSpheregR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLSpheregR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLSpheregR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLSpheregR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Sphere> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Sphere>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLSpheregR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Sphere>*)0x0)->GetClass();
      vectorlEgeoalgocLcLSpheregR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLSpheregR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLSpheregR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Sphere> : new vector<geoalgo::Sphere>;
   }
   static void *newArray_vectorlEgeoalgocLcLSpheregR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Sphere>[nElements] : new vector<geoalgo::Sphere>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLSpheregR(void *p) {
      delete ((vector<geoalgo::Sphere>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLSpheregR(void *p) {
      delete [] ((vector<geoalgo::Sphere>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLSpheregR(void *p) {
      typedef vector<geoalgo::Sphere> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Sphere>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLLineSegmentgR_Dictionary();
   static void vectorlEgeoalgocLcLLineSegmentgR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLLineSegmentgR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLLineSegmentgR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLLineSegmentgR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLLineSegmentgR(void *p);
   static void destruct_vectorlEgeoalgocLcLLineSegmentgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::LineSegment>*)
   {
      vector<geoalgo::LineSegment> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::LineSegment>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::LineSegment>", -2, "vector", 214,
                  typeid(vector<geoalgo::LineSegment>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLLineSegmentgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::LineSegment>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLLineSegmentgR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLLineSegmentgR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLLineSegmentgR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLLineSegmentgR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLLineSegmentgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::LineSegment> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::LineSegment>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLLineSegmentgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::LineSegment>*)0x0)->GetClass();
      vectorlEgeoalgocLcLLineSegmentgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLLineSegmentgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLLineSegmentgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::LineSegment> : new vector<geoalgo::LineSegment>;
   }
   static void *newArray_vectorlEgeoalgocLcLLineSegmentgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::LineSegment>[nElements] : new vector<geoalgo::LineSegment>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLLineSegmentgR(void *p) {
      delete ((vector<geoalgo::LineSegment>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLLineSegmentgR(void *p) {
      delete [] ((vector<geoalgo::LineSegment>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLLineSegmentgR(void *p) {
      typedef vector<geoalgo::LineSegment> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::LineSegment>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLLinegR_Dictionary();
   static void vectorlEgeoalgocLcLLinegR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLLinegR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLLinegR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLLinegR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLLinegR(void *p);
   static void destruct_vectorlEgeoalgocLcLLinegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Line>*)
   {
      vector<geoalgo::Line> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Line>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Line>", -2, "vector", 214,
                  typeid(vector<geoalgo::Line>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLLinegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Line>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLLinegR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLLinegR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLLinegR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLLinegR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLLinegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Line> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Line>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLLinegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Line>*)0x0)->GetClass();
      vectorlEgeoalgocLcLLinegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLLinegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLLinegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Line> : new vector<geoalgo::Line>;
   }
   static void *newArray_vectorlEgeoalgocLcLLinegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Line>[nElements] : new vector<geoalgo::Line>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLLinegR(void *p) {
      delete ((vector<geoalgo::Line>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLLinegR(void *p) {
      delete [] ((vector<geoalgo::Line>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLLinegR(void *p) {
      typedef vector<geoalgo::Line> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Line>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLHalfLinegR_Dictionary();
   static void vectorlEgeoalgocLcLHalfLinegR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLHalfLinegR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLHalfLinegR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLHalfLinegR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLHalfLinegR(void *p);
   static void destruct_vectorlEgeoalgocLcLHalfLinegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::HalfLine>*)
   {
      vector<geoalgo::HalfLine> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::HalfLine>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::HalfLine>", -2, "vector", 214,
                  typeid(vector<geoalgo::HalfLine>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLHalfLinegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::HalfLine>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLHalfLinegR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLHalfLinegR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLHalfLinegR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLHalfLinegR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLHalfLinegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::HalfLine> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::HalfLine>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLHalfLinegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::HalfLine>*)0x0)->GetClass();
      vectorlEgeoalgocLcLHalfLinegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLHalfLinegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLHalfLinegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::HalfLine> : new vector<geoalgo::HalfLine>;
   }
   static void *newArray_vectorlEgeoalgocLcLHalfLinegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::HalfLine>[nElements] : new vector<geoalgo::HalfLine>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLHalfLinegR(void *p) {
      delete ((vector<geoalgo::HalfLine>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLHalfLinegR(void *p) {
      delete [] ((vector<geoalgo::HalfLine>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLHalfLinegR(void *p) {
      typedef vector<geoalgo::HalfLine> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::HalfLine>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLDirectedLinegR_Dictionary();
   static void vectorlEgeoalgocLcLDirectedLinegR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLDirectedLinegR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLDirectedLinegR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLDirectedLinegR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLDirectedLinegR(void *p);
   static void destruct_vectorlEgeoalgocLcLDirectedLinegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::DirectedLine>*)
   {
      vector<geoalgo::DirectedLine> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::DirectedLine>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::DirectedLine>", -2, "vector", 214,
                  typeid(vector<geoalgo::DirectedLine>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLDirectedLinegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::DirectedLine>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLDirectedLinegR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLDirectedLinegR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLDirectedLinegR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLDirectedLinegR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLDirectedLinegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::DirectedLine> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::DirectedLine>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLDirectedLinegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::DirectedLine>*)0x0)->GetClass();
      vectorlEgeoalgocLcLDirectedLinegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLDirectedLinegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLDirectedLinegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::DirectedLine> : new vector<geoalgo::DirectedLine>;
   }
   static void *newArray_vectorlEgeoalgocLcLDirectedLinegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::DirectedLine>[nElements] : new vector<geoalgo::DirectedLine>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLDirectedLinegR(void *p) {
      delete ((vector<geoalgo::DirectedLine>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLDirectedLinegR(void *p) {
      delete [] ((vector<geoalgo::DirectedLine>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLDirectedLinegR(void *p) {
      typedef vector<geoalgo::DirectedLine> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::DirectedLine>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLCylindergR_Dictionary();
   static void vectorlEgeoalgocLcLCylindergR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLCylindergR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLCylindergR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLCylindergR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLCylindergR(void *p);
   static void destruct_vectorlEgeoalgocLcLCylindergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Cylinder>*)
   {
      vector<geoalgo::Cylinder> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Cylinder>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Cylinder>", -2, "vector", 214,
                  typeid(vector<geoalgo::Cylinder>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLCylindergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Cylinder>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLCylindergR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLCylindergR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLCylindergR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLCylindergR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLCylindergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Cylinder> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Cylinder>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLCylindergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Cylinder>*)0x0)->GetClass();
      vectorlEgeoalgocLcLCylindergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLCylindergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLCylindergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Cylinder> : new vector<geoalgo::Cylinder>;
   }
   static void *newArray_vectorlEgeoalgocLcLCylindergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Cylinder>[nElements] : new vector<geoalgo::Cylinder>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLCylindergR(void *p) {
      delete ((vector<geoalgo::Cylinder>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLCylindergR(void *p) {
      delete [] ((vector<geoalgo::Cylinder>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLCylindergR(void *p) {
      typedef vector<geoalgo::Cylinder> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Cylinder>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLConegR_Dictionary();
   static void vectorlEgeoalgocLcLConegR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLConegR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLConegR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLConegR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLConegR(void *p);
   static void destruct_vectorlEgeoalgocLcLConegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::Cone>*)
   {
      vector<geoalgo::Cone> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::Cone>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::Cone>", -2, "vector", 214,
                  typeid(vector<geoalgo::Cone>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLConegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::Cone>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLConegR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLConegR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLConegR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLConegR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLConegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::Cone> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::Cone>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLConegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::Cone>*)0x0)->GetClass();
      vectorlEgeoalgocLcLConegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLConegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLConegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Cone> : new vector<geoalgo::Cone>;
   }
   static void *newArray_vectorlEgeoalgocLcLConegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::Cone>[nElements] : new vector<geoalgo::Cone>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLConegR(void *p) {
      delete ((vector<geoalgo::Cone>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLConegR(void *p) {
      delete [] ((vector<geoalgo::Cone>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLConegR(void *p) {
      typedef vector<geoalgo::Cone> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::Cone>

namespace ROOT {
   static TClass *vectorlEgeoalgocLcLAABoxgR_Dictionary();
   static void vectorlEgeoalgocLcLAABoxgR_TClassManip(TClass*);
   static void *new_vectorlEgeoalgocLcLAABoxgR(void *p = 0);
   static void *newArray_vectorlEgeoalgocLcLAABoxgR(Long_t size, void *p);
   static void delete_vectorlEgeoalgocLcLAABoxgR(void *p);
   static void deleteArray_vectorlEgeoalgocLcLAABoxgR(void *p);
   static void destruct_vectorlEgeoalgocLcLAABoxgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<geoalgo::AABox>*)
   {
      vector<geoalgo::AABox> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<geoalgo::AABox>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<geoalgo::AABox>", -2, "vector", 214,
                  typeid(vector<geoalgo::AABox>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEgeoalgocLcLAABoxgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<geoalgo::AABox>) );
      instance.SetNew(&new_vectorlEgeoalgocLcLAABoxgR);
      instance.SetNewArray(&newArray_vectorlEgeoalgocLcLAABoxgR);
      instance.SetDelete(&delete_vectorlEgeoalgocLcLAABoxgR);
      instance.SetDeleteArray(&deleteArray_vectorlEgeoalgocLcLAABoxgR);
      instance.SetDestructor(&destruct_vectorlEgeoalgocLcLAABoxgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<geoalgo::AABox> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<geoalgo::AABox>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEgeoalgocLcLAABoxgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<geoalgo::AABox>*)0x0)->GetClass();
      vectorlEgeoalgocLcLAABoxgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEgeoalgocLcLAABoxgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEgeoalgocLcLAABoxgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::AABox> : new vector<geoalgo::AABox>;
   }
   static void *newArray_vectorlEgeoalgocLcLAABoxgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<geoalgo::AABox>[nElements] : new vector<geoalgo::AABox>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEgeoalgocLcLAABoxgR(void *p) {
      delete ((vector<geoalgo::AABox>*)p);
   }
   static void deleteArray_vectorlEgeoalgocLcLAABoxgR(void *p) {
      delete [] ((vector<geoalgo::AABox>*)p);
   }
   static void destruct_vectorlEgeoalgocLcLAABoxgR(void *p) {
      typedef vector<geoalgo::AABox> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<geoalgo::AABox>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 214,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *maplEgeoalgocLcLVectorcOstringgR_Dictionary();
   static void maplEgeoalgocLcLVectorcOstringgR_TClassManip(TClass*);
   static void *new_maplEgeoalgocLcLVectorcOstringgR(void *p = 0);
   static void *newArray_maplEgeoalgocLcLVectorcOstringgR(Long_t size, void *p);
   static void delete_maplEgeoalgocLcLVectorcOstringgR(void *p);
   static void deleteArray_maplEgeoalgocLcLVectorcOstringgR(void *p);
   static void destruct_maplEgeoalgocLcLVectorcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<geoalgo::Vector,string>*)
   {
      map<geoalgo::Vector,string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<geoalgo::Vector,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<geoalgo::Vector,string>", -2, "map", 96,
                  typeid(map<geoalgo::Vector,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEgeoalgocLcLVectorcOstringgR_Dictionary, isa_proxy, 4,
                  sizeof(map<geoalgo::Vector,string>) );
      instance.SetNew(&new_maplEgeoalgocLcLVectorcOstringgR);
      instance.SetNewArray(&newArray_maplEgeoalgocLcLVectorcOstringgR);
      instance.SetDelete(&delete_maplEgeoalgocLcLVectorcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEgeoalgocLcLVectorcOstringgR);
      instance.SetDestructor(&destruct_maplEgeoalgocLcLVectorcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<geoalgo::Vector,string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<geoalgo::Vector,string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEgeoalgocLcLVectorcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<geoalgo::Vector,string>*)0x0)->GetClass();
      maplEgeoalgocLcLVectorcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEgeoalgocLcLVectorcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEgeoalgocLcLVectorcOstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<geoalgo::Vector,string> : new map<geoalgo::Vector,string>;
   }
   static void *newArray_maplEgeoalgocLcLVectorcOstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<geoalgo::Vector,string>[nElements] : new map<geoalgo::Vector,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEgeoalgocLcLVectorcOstringgR(void *p) {
      delete ((map<geoalgo::Vector,string>*)p);
   }
   static void deleteArray_maplEgeoalgocLcLVectorcOstringgR(void *p) {
      delete [] ((map<geoalgo::Vector,string>*)p);
   }
   static void destruct_maplEgeoalgocLcLVectorcOstringgR(void *p) {
      typedef map<geoalgo::Vector,string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<geoalgo::Vector,string>

namespace {
  void TriggerDictionaryInitialization_liblarcorealg_GeoAlgo_Impl() {
    static const char* headers[] = {
"GeoAABox.h",
"GeoAlgoConstants.h",
"GeoAlgoException.h",
"GeoAlgo.h",
"GeoAlgo-TypeDef.h",
"GeoCone.h",
"GeoCylinder.h",
"GeoDirectedLine.h",
"GeoHalfLine.h",
"GeoLine.h",
"GeoLineSegment.h",
"GeoObjCollection.h",
"GeoSphere.h",
"GeoTrajectory.h",
"GeoVector.h",
0
    };
    static const char* includePaths[] = {
"/products/root/v6_06_04b/Linux64bit+2.6-2.12-e10-prof/include",
"/scratch/lsimons/LArLite_LArSoft_Int/larlite/UserDev/larcorealg/larcorealg/GeoAlgo/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "liblarcorealg_GeoAlgo dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAABox.h")))  Vector;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  Trajectory;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAABox.h")))  HalfLine;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  Line;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoDirectedLine.h")))  DirectedLine;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  LineSegment;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAABox.h")))  AABox;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoCylinder.h")))  Cylinder;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  Cone;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  Sphere;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoAlgo.h")))  GeoAlgo;}
namespace geoalgo{class __attribute__((annotate("$clingAutoload$GeoObjCollection.h")))  GeoObjCollection;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "liblarcorealg_GeoAlgo dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "GeoAABox.h"
#include "GeoAlgoConstants.h"
#include "GeoAlgoException.h"
#include "GeoAlgo.h"
#include "GeoAlgo-TypeDef.h"
#include "GeoCone.h"
#include "GeoCylinder.h"
#include "GeoDirectedLine.h"
#include "GeoHalfLine.h"
#include "GeoLine.h"
#include "GeoLineSegment.h"
#include "GeoObjCollection.h"
#include "GeoSphere.h"
#include "GeoTrajectory.h"
#include "GeoVector.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"geoalgo::AABox", payloadCode, "@",
"geoalgo::Cone", payloadCode, "@",
"geoalgo::Cylinder", payloadCode, "@",
"geoalgo::DirectedLine", payloadCode, "@",
"geoalgo::GeoAlgo", payloadCode, "@",
"geoalgo::GeoObjCollection", payloadCode, "@",
"geoalgo::HalfLine", payloadCode, "@",
"geoalgo::Line", payloadCode, "@",
"geoalgo::LineSegment", payloadCode, "@",
"geoalgo::Sphere", payloadCode, "@",
"geoalgo::Trajectory", payloadCode, "@",
"geoalgo::Vector", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("liblarcorealg_GeoAlgo",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_liblarcorealg_GeoAlgo_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_liblarcorealg_GeoAlgo_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_liblarcorealg_GeoAlgo() {
  TriggerDictionaryInitialization_liblarcorealg_GeoAlgo_Impl();
}
