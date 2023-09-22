#include <iostream>
#include <immintrin.h>  // SSE intrinsics header

// Provided implementation of popcount_3 function
unsigned popcount_3(uint64_t x)
{
    const uint64_t m1  = 0x5555555555555555;
    const uint64_t m2  = 0x3333333333333333;
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
    const uint64_t h01 = 0x0101010101010101;

    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    return (x * h01) >> 56;
}

template<typename _val>
unsigned match_block(const _val *x, const _val *y)
{
    static const __m128i mask = _mm_set1_epi8(0x7F);
    __m128i r1 = _mm_loadu_si128 ((__m128i const*)(x));
    __m128i r2 = _mm_loadu_si128 ((__m128i const*)(y));
    r2 = _mm_and_si128(r2, mask);
    return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
}

template<typename _val>
unsigned fast_match(const _val *q, const _val *s)
{ return popcount_3(match_block(q-8, s-8)<<16 | match_block(q+8, s+8)); }

int main()
{
    // Example data: integer values
    int q_values[] = {12, 34, 56, 78, 90, 12, 34, 56};
    int s_values[] = {12, 34, 56, 78, 91, 12, 34, 56};

    // Call the fast_match function and store the result
    unsigned similarity = fast_match(q_values, s_values);

    // Expected result based on the example data and logic
    unsigned expected_similarity = 4;

    // Check if the calculated similarity matches the expected result
    if (similarity == expected_similarity) {
        std::cout << "Test case passed. Similarity: " << similarity << std::endl;
    } else {
        std::cout << "Test case failed. Expected: " << expected_similarity << ", Actual: " << similarity << std::endl;
    }

    return 0;
}
