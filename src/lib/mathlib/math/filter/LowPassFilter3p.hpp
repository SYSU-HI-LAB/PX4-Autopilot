/****************************************************************************
 *
 *   Copyright (C) 2012-2021 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/// @file	LowPassFilter3p.hpp
/// @brief	A class to implement a third order butterworth low pass filter
/// Author: LYU Ximin <lvxm6@mail.sysu.edu.cn>
///

#pragma once

#include <mathlib/math/Functions.hpp>
#include <float.h>
#include <matrix/math.hpp>

namespace math
{

template<typename T>
class LowPassFilter3p
{
public:
	LowPassFilter3p() = default;

	LowPassFilter3p(float sample_freq, float cutoff_freq)
	{
		// set initial parameters
		set_cutoff_frequency(sample_freq, cutoff_freq);
	}

	// Change filter parameters
	void set_cutoff_frequency(float sample_freq, float cutoff_freq)
	{
		if ((sample_freq <= 0.f) || (cutoff_freq <= 0.f) || (cutoff_freq >= sample_freq / 2)
		    || !isFinite(sample_freq) || !isFinite(cutoff_freq)) {

			disable();
			return;
		}

		// reset delay elements on filter change
		_delay_element_1 = {};
		_delay_element_2 = {};
		_delay_element_3 = {};

		_cutoff_freq = math::max(cutoff_freq, sample_freq * 0.001f);
		_sample_freq = sample_freq;

		const float wch = M_PI_F*_cutoff_freq/_sample_freq;
		const float tn3 = tan(wch)*tan(wch)*tan(wch);
		const float tn2 = tan(wch)*tan(wch);
		const float tn 	= tan(wch);
		const float a0 = 2*tn2+2*tn+tn3+1;

		_a1 = (2*tn2-2*tn+3*tn3-3)/a0;
		_a2 = (-2*tn2-2*tn+3*tn3+3)/a0;
		_a3 = (-2*tn2+2*tn+tn3-1)/a0;

		_b0 = tn3/a0;
		_b3 = _b0;
		_b1 = 3*tn3/a0;
		_b2 = _b1;

		if (!isFinite(_b0) || !isFinite(_b1) || !isFinite(_b2) || !isFinite(_b3) || !isFinite(_a1) || !isFinite(_a2) || !isFinite(_a3)) {
			disable();
		}
	}

	/**
	 * Add a new raw value to the filter
	 *
	 * @return retrieve the filtered result
	 */
	inline T apply(const T &sample)
	{
		// Direct Form II implementation
		T delay_element_0{sample - _delay_element_1 *_a1 - _delay_element_2 * _a2 - _delay_element_3 * _a3};

		const T output{delay_element_0 *_b0 + _delay_element_1 *_b1 + _delay_element_2 * _b2 + _delay_element_3 * _b3};

		_delay_element_3 = _delay_element_2;
		_delay_element_2 = _delay_element_1;
		_delay_element_1 = delay_element_0;

		return output;
	}

	// Filter array of samples in place using the Direct form II.
	inline void applyArray(T samples[], int num_samples)
	{
		for (int n = 0; n < num_samples; n++) {
			samples[n] = apply(samples[n]);
		}
	}

	// Return the cutoff frequency
	float get_cutoff_freq() const { return _cutoff_freq; }

	// Return the sample frequency
	float get_sample_freq() const { return _sample_freq; }

	float getMagnitudeResponse(float frequency) const;

	// Reset the filter state to this value
	T reset(const T &sample)
	{
		const T input = isFinite(sample) ? sample : T{};

		if (fabsf(1 + _a1 + _a2 + _a3) > FLT_EPSILON) {
			_delay_element_1 = _delay_element_2 = _delay_element_3 = input / (1 + _a1 + _a2 + _a3);

			if (!isFinite(_delay_element_1) || !isFinite(_delay_element_2) || !isFinite(_delay_element_3)) {
				_delay_element_1 = _delay_element_2 = _delay_element_3 = input;
			}

		} else {
			_delay_element_1 = _delay_element_2 = _delay_element_3 = input;
		}

		return apply(input);
	}

	void disable()
	{
		// no filtering
		_sample_freq = 0.f;
		_cutoff_freq = 0.f;

		_delay_element_1 = {};
		_delay_element_2 = {};
		_delay_element_3 = {};

		_b0 = 1.f;
		_b1 = 0.f;
		_b2 = 0.f;
		_b3 = 0.f;

		_a1 = 0.f;
		_a2 = 0.f;
		_a3 = 0.f;
	}

protected:
	T _delay_element_1{}; // buffered sample -1
	T _delay_element_2{}; // buffered sample -2
	T _delay_element_3{}; // buffered sample -3

	// All the coefficients are normalized by a0, so a0 becomes 1 here
	float _a1{0.f};
	float _a2{0.f};
	float _a3{0.f};

	float _b0{1.f};
	float _b1{0.f};
	float _b2{0.f};
	float _b3{0.f};

	float _cutoff_freq{0.f};
	float _sample_freq{0.f};
};

} // namespace math
