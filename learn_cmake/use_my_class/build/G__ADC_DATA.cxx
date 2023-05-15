// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__ADC_DATA
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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

// Header files passed as explicit arguments
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/ADC_DATA.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *ADC_DATA_Dictionary();
   static void ADC_DATA_TClassManip(TClass*);
   static void *new_ADC_DATA(void *p = nullptr);
   static void *newArray_ADC_DATA(Long_t size, void *p);
   static void delete_ADC_DATA(void *p);
   static void deleteArray_ADC_DATA(void *p);
   static void destruct_ADC_DATA(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ADC_DATA*)
   {
      ::ADC_DATA *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ADC_DATA));
      static ::ROOT::TGenericClassInfo 
         instance("ADC_DATA", "", 14,
                  typeid(::ADC_DATA), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ADC_DATA_Dictionary, isa_proxy, 4,
                  sizeof(::ADC_DATA) );
      instance.SetNew(&new_ADC_DATA);
      instance.SetNewArray(&newArray_ADC_DATA);
      instance.SetDelete(&delete_ADC_DATA);
      instance.SetDeleteArray(&deleteArray_ADC_DATA);
      instance.SetDestructor(&destruct_ADC_DATA);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ADC_DATA*)
   {
      return GenerateInitInstanceLocal((::ADC_DATA*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ADC_DATA*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ADC_DATA_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ADC_DATA*)nullptr)->GetClass();
      ADC_DATA_TClassManip(theClass);
   return theClass;
   }

   static void ADC_DATA_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_ADC_DATA(void *p) {
      return  p ? new(p) ::ADC_DATA : new ::ADC_DATA;
   }
   static void *newArray_ADC_DATA(Long_t nElements, void *p) {
      return p ? new(p) ::ADC_DATA[nElements] : new ::ADC_DATA[nElements];
   }
   // Wrapper around operator delete
   static void delete_ADC_DATA(void *p) {
      delete ((::ADC_DATA*)p);
   }
   static void deleteArray_ADC_DATA(void *p) {
      delete [] ((::ADC_DATA*)p);
   }
   static void destruct_ADC_DATA(void *p) {
      typedef ::ADC_DATA current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ADC_DATA

namespace {
  void TriggerDictionaryInitialization_libADC_DATA_Impl() {
    static const char* headers[] = {
"/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/ADC_DATA.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/henrique/Documents/root_6.26.04/root_install/include",
"/home/henrique/Dropbox/APC_Paris/Root/learn_cmake/use_my_class",
"/home/henrique/Documents/root_6.26.04/root_install/include/",
"/home/henrique/Dropbox/APC_Paris/Root/learn_cmake/use_my_class/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libADC_DATA dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/ADC_DATA.h")))  ADC_DATA;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libADC_DATA dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/ADC_DATA.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ADC_DATA", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libADC_DATA",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libADC_DATA_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libADC_DATA_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libADC_DATA() {
  TriggerDictionaryInitialization_libADC_DATA_Impl();
}
