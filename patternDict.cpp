//
// File generated by rootcint at Sat Apr  5 13:35:53 2014

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME patternDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "patternDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void input_sample_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_input_sample(void *p = 0);
   static void *newArray_input_sample(Long_t size, void *p);
   static void delete_input_sample(void *p);
   static void deleteArray_input_sample(void *p);
   static void destruct_input_sample(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::input_sample*)
   {
      ::input_sample *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::input_sample >(0);
      static ::ROOT::TGenericClassInfo 
         instance("input_sample", ::input_sample::Class_Version(), "./pattern.h", 3,
                  typeid(::input_sample), DefineBehavior(ptr, ptr),
                  &::input_sample::Dictionary, isa_proxy, 4,
                  sizeof(::input_sample) );
      instance.SetNew(&new_input_sample);
      instance.SetNewArray(&newArray_input_sample);
      instance.SetDelete(&delete_input_sample);
      instance.SetDeleteArray(&deleteArray_input_sample);
      instance.SetDestructor(&destruct_input_sample);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::input_sample*)
   {
      return GenerateInitInstanceLocal((::input_sample*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::input_sample*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void solve_data_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_solve_data(void *p = 0);
   static void *newArray_solve_data(Long_t size, void *p);
   static void delete_solve_data(void *p);
   static void deleteArray_solve_data(void *p);
   static void destruct_solve_data(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::solve_data*)
   {
      ::solve_data *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::solve_data >(0);
      static ::ROOT::TGenericClassInfo 
         instance("solve_data", ::solve_data::Class_Version(), "./pattern.h", 16,
                  typeid(::solve_data), DefineBehavior(ptr, ptr),
                  &::solve_data::Dictionary, isa_proxy, 4,
                  sizeof(::solve_data) );
      instance.SetNew(&new_solve_data);
      instance.SetNewArray(&newArray_solve_data);
      instance.SetDelete(&delete_solve_data);
      instance.SetDeleteArray(&deleteArray_solve_data);
      instance.SetDestructor(&destruct_solve_data);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::solve_data*)
   {
      return GenerateInitInstanceLocal((::solve_data*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::solve_data*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *input_sample::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *input_sample::Class_Name()
{
   return "input_sample";
}

//______________________________________________________________________________
const char *input_sample::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::input_sample*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int input_sample::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::input_sample*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void input_sample::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::input_sample*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *input_sample::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::input_sample*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *solve_data::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *solve_data::Class_Name()
{
   return "solve_data";
}

//______________________________________________________________________________
const char *solve_data::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::solve_data*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int solve_data::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::solve_data*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void solve_data::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::solve_data*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *solve_data::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::solve_data*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void solve_data::Streamer(TBuffer &R__b)
{
   // Stream an object of class solve_data.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(solve_data::Class(),this);
   } else {
      R__b.WriteClassBuffer(solve_data::Class(),this);
   }
}

//______________________________________________________________________________
void solve_data::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class solve_data.
      TClass *R__cl = ::solve_data::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "inet", &inet);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "params", &params);
      R__insp.InspectMember(params, "params.");
      R__insp.Inspect(R__cl, R__insp.GetParent(), "data[800]", data);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_solve_data(void *p) {
      return  p ? new(p) ::solve_data : new ::solve_data;
   }
   static void *newArray_solve_data(Long_t nElements, void *p) {
      return p ? new(p) ::solve_data[nElements] : new ::solve_data[nElements];
   }
   // Wrapper around operator delete
   static void delete_solve_data(void *p) {
      delete ((::solve_data*)p);
   }
   static void deleteArray_solve_data(void *p) {
      delete [] ((::solve_data*)p);
   }
   static void destruct_solve_data(void *p) {
      typedef ::solve_data current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::solve_data

//______________________________________________________________________________
void input_sample::Streamer(TBuffer &R__b)
{
   // Stream an object of class input_sample.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(input_sample::Class(),this);
   } else {
      R__b.WriteClassBuffer(input_sample::Class(),this);
   }
}

//______________________________________________________________________________
void input_sample::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class input_sample.
      TClass *R__cl = ::input_sample::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "t_a", &t_a);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "t_b", &t_b);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "k[8]", k);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "n[8]", n);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "alpha", &alpha);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "beta", &beta);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_input_sample(void *p) {
      return  p ? new(p) ::input_sample : new ::input_sample;
   }
   static void *newArray_input_sample(Long_t nElements, void *p) {
      return p ? new(p) ::input_sample[nElements] : new ::input_sample[nElements];
   }
   // Wrapper around operator delete
   static void delete_input_sample(void *p) {
      delete ((::input_sample*)p);
   }
   static void deleteArray_input_sample(void *p) {
      delete [] ((::input_sample*)p);
   }
   static void destruct_input_sample(void *p) {
      typedef ::input_sample current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::input_sample

/********************************************************
* patternDict.cpp
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtablepatternDict();

extern "C" void G__set_cpp_environmentpatternDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("pattern.h");
  G__cpp_reset_tagtablepatternDict();
}
#include <new>
extern "C" int G__cpp_dllrevpatternDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* input_sample */
static int G__patternDict_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   input_sample* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new input_sample[n];
     } else {
       p = new((void*) gvp) input_sample[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new input_sample;
     } else {
       p = new((void*) gvp) input_sample;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__patternDictLN_input_sample));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) input_sample::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) input_sample::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) input_sample::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      input_sample::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((input_sample*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) input_sample::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) input_sample::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) input_sample::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_168_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) input_sample::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__patternDict_168_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   input_sample* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new input_sample(*(input_sample*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__patternDictLN_input_sample));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef input_sample G__Tinput_sample;
static int G__patternDict_168_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (input_sample*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((input_sample*) (soff+(sizeof(input_sample)*i)))->~G__Tinput_sample();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (input_sample*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((input_sample*) (soff))->~G__Tinput_sample();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__patternDict_168_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   input_sample* dest = (input_sample*) G__getstructoffset();
   *dest = *(input_sample*) libp->para[0].ref;
   const input_sample& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* solve_data */
static int G__patternDict_169_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   solve_data* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new solve_data[n];
     } else {
       p = new((void*) gvp) solve_data[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new solve_data;
     } else {
       p = new((void*) gvp) solve_data;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__patternDictLN_solve_data));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) solve_data::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) solve_data::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) solve_data::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      solve_data::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((solve_data*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) solve_data::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) solve_data::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) solve_data::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__patternDict_169_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) solve_data::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__patternDict_169_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   solve_data* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new solve_data(*(solve_data*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__patternDictLN_solve_data));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef solve_data G__Tsolve_data;
static int G__patternDict_169_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (solve_data*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((solve_data*) (soff+(sizeof(solve_data)*i)))->~G__Tsolve_data();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (solve_data*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((solve_data*) (soff))->~G__Tsolve_data();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__patternDict_169_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   solve_data* dest = (solve_data*) G__getstructoffset();
   *dest = *(solve_data*) libp->para[0].ref;
   const solve_data& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* input_sample */

/* solve_data */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncpatternDict {
 public:
  G__Sizep2memfuncpatternDict(): p(&G__Sizep2memfuncpatternDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncpatternDict::*p)();
};

size_t G__get_sizep2memfuncpatternDict()
{
  G__Sizep2memfuncpatternDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritancepatternDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__patternDictLN_input_sample))) {
     input_sample *G__Lderived;
     G__Lderived=(input_sample*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__patternDictLN_input_sample),G__get_linked_tagnum(&G__patternDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__patternDictLN_solve_data))) {
     solve_data *G__Lderived;
     G__Lderived=(solve_data*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__patternDictLN_solve_data),G__get_linked_tagnum(&G__patternDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetablepatternDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__patternDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__patternDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__patternDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__patternDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* input_sample */
static void G__setup_memvarinput_sample(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__patternDictLN_input_sample));
   { input_sample *p; p=(input_sample*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->t_a)-(long)(p)),100,0,0,-1,-1,-1,1,"t_a=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->t_b)-(long)(p)),100,0,0,-1,-1,-1,1,"t_b=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->k)-(long)(p)),100,0,0,-1,-1,-1,1,"k[8]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->n)-(long)(p)),100,0,0,-1,-1,-1,1,"n[8]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->alpha)-(long)(p)),100,0,0,-1,-1,-1,1,"alpha=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->beta)-(long)(p)),100,0,0,-1,-1,-1,1,"beta=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__patternDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* solve_data */
static void G__setup_memvarsolve_data(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__patternDictLN_solve_data));
   { solve_data *p; p=(solve_data*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->inet)-(long)(p)),105,0,0,-1,-1,-1,1,"inet=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->params)-(long)(p)),117,0,0,G__get_linked_tagnum(&G__patternDictLN_input_sample),-1,-1,1,"params=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->data)-(long)(p)),100,0,0,-1,-1,-1,1,"data[800]=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__patternDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarpatternDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncinput_sample(void) {
   /* input_sample */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__patternDictLN_input_sample));
   G__memfunc_setup("input_sample",1297,G__patternDict_168_0_1, 105, G__get_linked_tagnum(&G__patternDictLN_input_sample), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__patternDict_168_0_2, 85, G__get_linked_tagnum(&G__patternDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&input_sample::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__patternDict_168_0_3, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&input_sample::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__patternDict_168_0_4, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&input_sample::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__patternDict_168_0_5, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&input_sample::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__patternDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__patternDict_168_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__patternDict_168_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&input_sample::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__patternDict_168_0_11, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&input_sample::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__patternDict_168_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&input_sample::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__patternDict_168_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&input_sample::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("input_sample", 1297, G__patternDict_168_0_14, (int) ('i'), G__get_linked_tagnum(&G__patternDictLN_input_sample), -1, 0, 1, 1, 1, 0, "u 'input_sample' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~input_sample", 1423, G__patternDict_168_0_15, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__patternDict_168_0_16, (int) ('u'), G__get_linked_tagnum(&G__patternDictLN_input_sample), -1, 1, 1, 1, 1, 0, "u 'input_sample' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncsolve_data(void) {
   /* solve_data */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__patternDictLN_solve_data));
   G__memfunc_setup("solve_data",1058,G__patternDict_169_0_1, 105, G__get_linked_tagnum(&G__patternDictLN_solve_data), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__patternDict_169_0_2, 85, G__get_linked_tagnum(&G__patternDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&solve_data::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__patternDict_169_0_3, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&solve_data::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__patternDict_169_0_4, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&solve_data::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__patternDict_169_0_5, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&solve_data::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__patternDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__patternDict_169_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__patternDict_169_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&solve_data::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__patternDict_169_0_11, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&solve_data::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__patternDict_169_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&solve_data::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__patternDict_169_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&solve_data::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("solve_data", 1058, G__patternDict_169_0_14, (int) ('i'), G__get_linked_tagnum(&G__patternDictLN_solve_data), -1, 0, 1, 1, 1, 0, "u 'solve_data' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~solve_data", 1184, G__patternDict_169_0_15, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__patternDict_169_0_16, (int) ('u'), G__get_linked_tagnum(&G__patternDictLN_solve_data), -1, 1, 1, 1, 1, 0, "u 'solve_data' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncpatternDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalpatternDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcpatternDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__patternDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__patternDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__patternDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__patternDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__patternDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__patternDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__patternDictLN_input_sample = { "input_sample" , 99 , -1 };
G__linked_taginfo G__patternDictLN_solve_data = { "solve_data" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtablepatternDict() {
  G__patternDictLN_TClass.tagnum = -1 ;
  G__patternDictLN_TBuffer.tagnum = -1 ;
  G__patternDictLN_TMemberInspector.tagnum = -1 ;
  G__patternDictLN_TObject.tagnum = -1 ;
  G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__patternDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__patternDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__patternDictLN_input_sample.tagnum = -1 ;
  G__patternDictLN_solve_data.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtablepatternDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__patternDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__patternDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__patternDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__patternDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__patternDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__patternDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__patternDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__patternDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__patternDictLN_input_sample),sizeof(input_sample),-1,291072,(char*)NULL,G__setup_memvarinput_sample,G__setup_memfuncinput_sample);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__patternDictLN_solve_data),sizeof(solve_data),-1,291072,(char*)NULL,G__setup_memvarsolve_data,G__setup_memfuncsolve_data);
}
extern "C" void G__cpp_setuppatternDict(void) {
  G__check_setup_version(30051515,"G__cpp_setuppatternDict()");
  G__set_cpp_environmentpatternDict();
  G__cpp_setup_tagtablepatternDict();

  G__cpp_setup_inheritancepatternDict();

  G__cpp_setup_typetablepatternDict();

  G__cpp_setup_memvarpatternDict();

  G__cpp_setup_memfuncpatternDict();
  G__cpp_setup_globalpatternDict();
  G__cpp_setup_funcpatternDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncpatternDict();
  return;
}
class G__cpp_setup_initpatternDict {
  public:
    G__cpp_setup_initpatternDict() { G__add_setup_func("patternDict",(G__incsetup)(&G__cpp_setuppatternDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initpatternDict() { G__remove_setup_func("patternDict"); }
};
G__cpp_setup_initpatternDict G__cpp_setup_initializerpatternDict;

