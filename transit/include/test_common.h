// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _TEST_COMMON_H
#define _TEST_COMMON_H

#define test_result(...)        do{                 \
         fprintf(stdout,__VA_ARGS__);  \
                                       }while(0)

#endif
