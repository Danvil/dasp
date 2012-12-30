/*
 * SolveDenseTemplate.cpp
 *
 *  Created on: Okt 20, 2012
 *      Author: david
 */

#include <arpack++/arlssym.h>
#include <iostream>

void MemoryOverflow()
{
	std::cerr << "ArpackError: MEMORY_OVERFLOW" << std::endl;
	throw ArpackError(ArpackError::MEMORY_OVERFLOW);
}

void ArpackError::Set(ArpackError::ErrorCode code, char const* msg)
{
//	ArpackError::code = code;
	std::cerr << "ArpackError: code=" << code << ", msg=" << std::string(msg) << std::endl;
}

