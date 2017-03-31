#ifndef _DEBUGGINGTOOLS_HPP
#define _DEBUGGINGTOOLS_HPP

#include <string.h>
#include "parstream.H" //Gives us pout()

//This file contains a collection of helpful #defines and other definitions
//that make debugging quicker.

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define DEBUG_SHOW(VAR) pout() << #VAR << ": " << VAR << " "
#define DEBUG_HEADER pout () << __FILENAME__ << ": "
#define DEBUG_END pout () << std::endl

///The macros DEBUG_OUT make debugging quicker and allow easy printing of a variable.
#define DEBUG_OUT(VAR) DEBUG_HEADER; DEBUG_SHOW(VAR); DEBUG_END
#define DEBUG_OUT2(VAR1,VAR2) DEBUG_HEADER; DEBUG_SHOW(VAR1); DEBUG_SHOW(VAR2); DEBUG_END
#define DEBUG_OUT3(VAR1,VAR2,VAR3) DEBUG_HEADER; DEBUG_SHOW(VAR1); DEBUG_SHOW(VAR2); DEBUG_SHOW(VAR3); DEBUG_END
#define DEBUG_OUT4(VAR1,VAR2,VAR3,VAR4) DEBUG_HEADER; DEBUG_SHOW(VAR1); DEBUG_SHOW(VAR2); DEBUG_SHOW(VAR3); DEBUG_SHOW(VAR4); DEBUG_END

#endif /* _DEBUGGINGTOOLS_HPP */
