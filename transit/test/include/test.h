/****************************** START LICENSE ******************************
Transit, a code to solve for the radiative-transifer equation for
planetary atmospheres.

This project was completed with the support of the NASA Planetary
Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
Joseph Harrington. Principal developers included graduate students
Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
undergraduate Andrew S. D. Foster.  The included
'transit' radiative transfer code is based on an earlier program of
the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
he was a graduate student at Cornell University under Joseph
Harrington.

Copyright (C) 2015 University of Central Florida.  All rights reserved.

This is a test version only, and may not be redistributed to any third
party.  Please refer such requests to us.  This program is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

Our intent is to release this software under an open-source,
reproducible-research license, once the code is mature and the first
research paper describing the code has been accepted for publication
in a peer-reviewed journal.  We are committed to development in the
open, and have posted this code on github.com so that others can test
it and give us feedback.  However, until its first publication and
first stable release, we do not permit others to redistribute the code
in either original or modified form, nor to publish work based in
whole or in part on the output of this code.  By downloading, running,
or modifying this code, you agree to these conditions.  We do
encourage sharing any modifications with us and discussing them
openly.

We welcome your feedback, but do not guarantee support.  Please send
feedback or inquiries to:

Joseph Harrington <jh@physics.ucf.edu>
Patricio Cubillos <pcubillos@fulbrightmail.org>
Jasmina Blecic <jasmina@physics.ucf.edu>

or alternatively,

Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
UCF PSB 441
4111 Libra Drive
Orlando, FL 32816-2385
USA

Thank you for using transit!
******************************* END LICENSE ******************************/


/* Test Framework for Transit
 * Created by AJ Foster <aj.foster@knights.ucf.edu>
 *
 * For more information, see the comments in test.c and the Wiki page located
 * at https://github.com/exosports/transit/wiki/Testing-Transit.
 *
 */

#ifndef _TEST_TRANSIT_H
#define _TEST_TRANSIT_H

/* Required for compilation with shared memory libraries.                   */
#define _XOPEN_SOURCE 700

#include <transit.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* The following are the return types for tests and test batches.
 *
 */
#define TR_TEST static char *
#define TR_BATCH int *


/* tr_setup_tests()
 *
 * Setup variables for the test suite. We would like to track the total number
 * of tests run as well as the number of failures. In the future, we may be
 * interested in the number of assertions made.
 *
 */
#define tr_setup_tests() int tr_num_tests = 0, tr_num_fails = 0


/* tr_run_test( test_name )
 *
 * Run a test function defined by test_name. This macro runs the function given
 * (without parens), which is expected to return a static character array or
 * NULL. Note that the test name is printed and flushed immediately in case the
 * test causes an error. Only one failing assertion is displayed at a time.
 *
 */
#define tr_run_test(test_name) do { \
  printf(#test_name ": "); \
  fflush(stdout); \
  char *msg = test_name(); \
  tr_num_tests++; \
  if (msg) { \
    printf("Fail\n  %s\n", msg); \
    tr_num_fails++; \
  } else { \
    printf("Pass\n"); \
  } \
} while(0)


/* tr_setup_batch()
 *
 * Setup variables for the test batch. This sets up variables necessary to help
 * track the number of tests and failures within the batch.
 *
 */
#define tr_setup_batch() int tr_num_tests = 0, tr_num_fails = 0


/* tr_finish_batch()
 *
 * Clean up and return the tracking numbers for this batch of tests.
 *
 */
#define tr_finish_batch() do { \
  int *_tr_return = malloc(sizeof(int) * 2); \
  _tr_return[0] = tr_num_tests; \
  _tr_return[1] = tr_num_fails; \
  return _tr_return; \
} while(0)


/* tr_run_batch( batch_name )
 *
 * Run a batch of tests with the given name. The batch is expected to be a
 * funcion with the TR_BATCH return type. Test batches should begin with the
 * tr_setup_batch() function and end with tr_finish_batch().
 *
 */
#define tr_run_batch(batch_name) do { \
 int *tr_batch_stats = batch_name(); \
 tr_num_tests += tr_batch_stats[0]; \
 tr_num_fails += tr_batch_stats[1]; \
 free(tr_batch_stats); \
} while(0)


/* tr_run_main( argc, argv )
 *
 * Run transit's main function with the given arguments. These arguments are
 * expected to be of the correct format. This avoids the issue of having two
 * main() functions.
 *
 */
#define tr_run_main(_argc, _argv) _tr_main(_argc, _argv)


/* tr_assert( assertion, message )
 *
 * Tests the truth value of the given assertion, which is expected to be a
 * boolean statement or value. If the test returns false, then the given message
 * is returned as a failed assertion.
 *
 */
#define tr_assert(assertion, msg) do { \
  if (!(assertion)) return msg; \
} while(0)


/* tr_assert_not( assertion, message )
 *
 * Tests the truth value of the given assertion, which is expected to be a
 * boolean statement or value. If the test returns true, then the given message
 * is returned as a failed assertion.
 *
 */
#define tr_assert_not(assertion, msg) do { \
  if (assertion) return msg; \
} while(0)


/* tr_assert_equal( value 1, value 2, message )
 *
 * Tests for integer equality between the given values. This is not intended to
 * be used for floating point values (see tr_assert_close) or strings (see
 * tr_assert_equal_str). Note that this test can be relevant for other boolean-
 * like values such as 0, 1, NULL, etc.
 *
 */
#define tr_assert_equal(val1, val2, msg) do { \
  if (!((val1) == (val2))) return msg; \
} while(0)


/* tr_assert_not_equal( value 1, value 2, message )
 *
 * Tests for integer inequality between the given values. This is not intended
 * to be used for floating point values or strings (see
 * tr_assert_not_equal_str). Note that this test can be relevant for other
 * boolean-like values such as 0, 1, NULL, etc.
 *
 */
#define tr_assert_not_equal(val1, val2, msg) do { \
  if ((val1) == (val2)) return msg; \
} while(0)


/* tr_assert_close( value 1, value 2, margin, message )
 *
 * Tests that the two values - assumed to be floating point or promotable - are
 * within the given margin (in absolute value).
 *
 */
#define tr_assert_close(val1, val2, margin, msg) do { \
  if (fabs((double)(val1) - (double)(val2)) > (double)(margin)) return msg; \
} while(0)


/* tr_assert_equal_str( string 1, string 2, message )
 *
 * Tests that the given strings (expected to be given as (char *)) are equal,
 * using strcmp(). If unequal, the given message is returned as a failed
 * assertion.
 *
 */
#define tr_assert_equal_str(str1, str2, msg) do { \
  if (strcmp((char *)(str1), (char *)(str2))) return msg; \
} while(0)


/* tr_assert_not_equal_str( string 1, string 2, message )
 *
 * Tests that the given strings (expected to be given as (char *)) are not
 * equal, using strcmp(). If equal, the given message is returned as a failed
 * assertion.
 *
 */
#define tr_assert_not_equal_str(str1, str2, msg) do { \
  if (!strcmp((char *)(str1), (char *)(str2))) return msg; \
} while(0)


/* tr_skip()
 *
 * Skips the remaining assertions in the test. Useful for debugging.
 *
 */
#define tr_skip() return NULL


/* tr_finish_tests()
 *
 * Wraps up the test process and prints summary information.
 *
 */
#define tr_finish_tests() do { \
  printf("\nSummary | Tests: %d, Failures: %d\n", \
    tr_num_tests, tr_num_fails); \
} while(0)


#ifndef TEST_TRANSIT

// Define a placeholder for the renamed main();
int _tr_main(int argc, char **argv) { return 0; }

#else

// Define transit's main function in its new name.
int _tr_main(int argc, char **argv);

#endif

#endif // _TEST_TRANSIT_H