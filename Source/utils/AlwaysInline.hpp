#ifndef ALWAYS_INLINE_HPP_
#define ALWAYS_INLINE_HPP_

#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((always_inline)) __inline__
#else
#define ALWAYS_INLINE inline
#endif

#endif /* ALWAYS_INLINE_HPP_ */
