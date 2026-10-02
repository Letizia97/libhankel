#ifndef LIBHANKEL_EXPORT_H
#define LIBHANKEL_EXPORT_H

/* On Windows the compiler only exports DLL symbols that are explicitly marked.
 * LIBHANKEL_BUILDING_DLL is defined by the build system when compiling the
 * library itself; consumer code that includes this header gets dllimport. */
#ifdef _WIN32
  #ifdef LIBHANKEL_BUILDING_DLL
    #define LIBHANKEL_API __declspec(dllexport)
  #else
    #define LIBHANKEL_API __declspec(dllimport)
  #endif
#else
  #define LIBHANKEL_API
#endif

#endif /* LIBHANKEL_EXPORT_H */
