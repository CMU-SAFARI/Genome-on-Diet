#if !defined(SIMDE_TESTS_ARM_INTERNAL_H)
#define SIMDE_TESTS_CURRENT_ARCH arm

#include "../run-tests.h"

HEDLEY_BEGIN_C_DECLS

#define SIMDE_TESTS_ARM_GENERATE_SUITE_GETTERS(isax) \
  MunitSuite* SIMDE_TESTS_GENERATE_SYMBOL_FULL(suite, SIMDE_TESTS_CURRENT_ARCH, isax, native, c)(void); \
  MunitSuite* SIMDE_TESTS_GENERATE_SYMBOL_FULL(suite, SIMDE_TESTS_CURRENT_ARCH, isax, native, cpp)(void); \
  MunitSuite* SIMDE_TESTS_GENERATE_SYMBOL_FULL(suite, SIMDE_TESTS_CURRENT_ARCH, isax, emul,   c)(void); \
  MunitSuite* SIMDE_TESTS_GENERATE_SYMBOL_FULL(suite, SIMDE_TESTS_CURRENT_ARCH, isax, emul,   cpp)(void)

SIMDE_TESTS_ARM_GENERATE_SUITE_GETTERS(neon);

HEDLEY_END_C_DECLS

#endif /* !defined(SIMDE_TESTS_ARM_INTERNAL_H) */
