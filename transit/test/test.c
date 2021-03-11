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
 * The following, along with include/test.h, creates a framework for testing
 * Transit. Below you can find a brief introduction to its usage. For more
 * information and examples, see:
 *   https://github.com/exosports/transit/wiki/Testing-Transit
 *
 * This test framework has several layers:
 *
 * - Assertions are individual checks that test, for example, the equality of
 *   two values. Assertions have messages attached to describe the failure.
 * - Tests are named collections of assertions. Think of tests as functions
 *   that run several assertions for a specific function or use case. Only one
 *   failing assertion will occur per test (the rest are skipped).
 * - Batches are collections of tests. Think of batches as functions that run
 *   multiple tests relating to the same .c file.
 *
 *
 * Format
 * ------
 *
 * The main() function of test.c must begin with tr_setup_tests() and end with
 * tr_finish_tests(). Inside of these, you can call test batches using
 * tr_run_batch(batch_name) or individual tests with tr_run_test(test_name).
 *
 * Test batches are functions, and must have the TR_BATCH return type. They
 * begin with tr_setup_batch() and end with tr_finish_batch(). Inside of a
 * batch, you can call individual tests using tr_run_test(test_name). Tests
 * should return NULL.
 *
 * Tests are also functions and should have the TR_TEST return type. Within a
 * test you can use any number of assertions. For a full list of assertions, see
 * the Wiki page or include/test.h.
 *
 *
 * Example
 * -------
 *
 * TR_TEST test_the_truth () {
 *   tr_assert(1, "One does not evaluate as true.");
 *   return NULL;
 * }
 *
 * TR_BATCH example_test_batch () {
 *   tr_setup_batch();
 *   tr_run_test(test_the_truth);
 *   tr_finish_batch();
 * }
 *
 * int main (...) {
 *   tr_setup_tests();
 *   tr_run_batch(example_test_batch);
 *   tr_finish_tests();
 *   return 0;
 * }
 *
 */

#include <test.h>

#ifdef TEST_TRANSIT
int main(int argc, char **argv) {
#else
int _tr_test(int argc, char **argv) {
#endif
  tr_setup_tests();

  // Define tests and batches to run here

  tr_finish_tests();
  return 0;
}