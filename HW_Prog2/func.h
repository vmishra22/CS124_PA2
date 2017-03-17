#pragma once

#ifndef _FUNC_
#define _FUNC_

#include <vector>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cassert>

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	assert(a.size() == b.size());

	std::vector<T> result;
	result.reserve(a.size());

	std::transform(a.begin(), a.end(), b.begin(),
		std::back_inserter(result), std::plus<T>());
	return result;
}

#endif
