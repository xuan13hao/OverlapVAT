#ifndef __SORT_UTIL_H__
#define __SORT_UTIL_H__



#include "Utils.h"
#include "../../common.h"

#ifdef AVX2

namespace avx2{
  // Regular
  template <typename InType, typename RegType>
  void SortBlock64(InType *&arr, size_t offset);
  template <typename InType, typename RegType>
  void SortBlock16(InType *&arr, size_t offset);

  // Masked
  template <typename InType, typename RegType>
  void MaskedSortBlock4x8(InType *&arr, size_t offset);
  template <typename InType, typename RegType>
  void MaskedSortBlock2x4(InType *&arr, size_t offset);

};

#endif
#endif // __SORT_UTIL_H__