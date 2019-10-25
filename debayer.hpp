/*
  Copyright (C) 2012 Hoyoung Lee

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <cassert>
#include <limits>
#include <iostream>
#include <cstring>
#include <float.h>

#include <cpixmap.hpp>
#include <cchunk.hpp>

/*
  template <typename T>
  void debayerBilinearGRBG2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  template <typename T>
  void debayerHiQLinearGRBG2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  template <typename T>
  void debayerEdgingGRBG2RGB(cpixmap<T>& bayer, int r_lower_limit, int r_higher_limit, int g_lower_limit, int g_higher_limit, int b_lower_limit, int b_higher_limit, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  template <typename T>
  void debayerBilinearRGGB2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  template <typename T>
  void debayerHiQLinearRGGB2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  template <typename T>
  void debayerEdgingRGGB2RGB(cpixmap<T>& bayer, int r_lower_limit, int r_higher_limit, int g_lower_limit, int g_higher_limit, int b_lower_limit, int b_higher_limit, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b);

  // http://www.analog.com/media/en/technical-documentation/application-notes/EE358.pdf
  template <typename T>
  void debayerRCCC2GrayR(cpixmap<T>& rccc, cpixmap<T>& gray, cpixmap<T>& r);
*/

template <typename T>
void debayerBilinearGRBG2RGB(const cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));
  
  cslice<T> linebuf(bayer, 2/*lines*/, 1/*hpadding*/, 1/*vpadding*/);
  linebuf.draftSlice(bayer);
  
  for (size_t y = 0; y < bayer.getHeight(); y += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	size_t i = y;

	//** GRGR...
	T *r_evenline = r.getLine(i, 0);
	T *g_evenline = g.getLine(i, 0);
	T *b_evenline = b.getLine(i, 0);
    
	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  // R on rGr
	  r_evenline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  // G on rGr
	  g_evenline[j] = linebuf(i, j);
	  // B on rGr
	  b_evenline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
      
	  j = x+1;
	  // R on gRg
	  r_evenline[j] = linebuf(i, j);
	  // G on gRg
#if 1
	  T hdiff = std::abs((int)linebuf(i, j-1) - (int)linebuf(i, j+1));
	  T vdiff = std::abs((int)linebuf(i-1, j) - (int)linebuf(i+1, j));
	  if (hdiff < vdiff) g_evenline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  else g_evenline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
#else
	  g_evenline[j] = (linebuf(i-1, j)>>2) + (linebuf(i+1, j)>>2) + (linebuf(i, j-1)>>2) + (linebuf(i, j+1)>>2);
#endif
	  // B on gRg
	  b_evenline[j] = (linebuf(i-1, j-1)>>2) + (linebuf(i-1, j+1)>>2) + (linebuf(i+1, j-1)>>2) + (linebuf(i+1, j+1)>>2);
	}
      }
# pragma omp section
      {
	size_t i = y+1;
    
	//** BGBG...
	T *r_oddline = r.getLine(i, 0);
	T *g_oddline = g.getLine(i, 0);
	T *b_oddline = b.getLine(i, 0);

	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  // R on gBg
	  r_oddline[j] = (linebuf(i-1, j-1) >> 2) + (linebuf(i-1, j+1) >> 2) + (linebuf(i+1, j-1) >> 2) + (linebuf(i+1, j+1) >> 2);
	  // G on gBg
#if 1
	  T hdiff = std::abs((int)linebuf(i, j-1) - (int)linebuf(i, j+1));
	  T vdiff = std::abs((int)linebuf(i-1, j) - (int)linebuf(i+1, j));
	  if (hdiff < vdiff) g_oddline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  else g_oddline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
#else
	  g_oddline[j] = (linebuf(i, j-1)>>2) + (linebuf(i, j+1)>>2) + (linebuf(i-1, j)>>2) + (linebuf(i+1, j)>>2);
#endif
	  // B on gBg
	  b_oddline[j] = linebuf(i, j);

	  j = x+1;
	  // R on bGb
	  r_oddline[j] = (linebuf(i-1, j) >> 1) + (linebuf(i+1, j) >> 1);
	  // G on bGb
	  g_oddline[j] = linebuf(i, j);
	  // B on bGb
	  b_oddline[j] = (linebuf(i, j-1) >> 1) + (linebuf(i, j+1) >> 1);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
}

/* Referenced from
 * "High-quality linear interpolation for demosaicing of bayer-patterned color images"
 * - Henrique S. Malvar, Li-wei He, and Ross Cutler
 */
template <typename T>
void debayerHiQLinearGRBG2RGB(const cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));

  cslice<T> linebuf(bayer, 2/*lines*/, 2/*hpadding*/, 2/*vpadding*/);
  linebuf.draftSlice(bayer);

  for (size_t y = 0; y < bayer.getHeight(); y += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	size_t i = y;
	
	// GRGR...
	T *r_evenline = r.getLine(i, 0);
	T *g_evenline = g.getLine(i, 0);
	T *b_evenline = b.getLine(i, 0);

	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  int rval, gval, bval;
	  // R on rGr
	  rval =
	    + ((int)linebuf(i-2, j)>>4)
	    - ((int)linebuf(i-1, j-1)>>3) - ((int)linebuf(i-1, j+1)>>3)
	    - ((int)linebuf(i, j-2)>>3) + ((int)linebuf(i, j-1)>>1) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>3) + ((int)linebuf(i, j+1)>>1) - ((int)linebuf(i, j+2)>>3)
	    - ((int)linebuf(i+1, j-1)>>3) - ((int)linebuf(i+1, j+1)>>3)
	    + ((int)linebuf(i+2, j)>>4);
	  r_evenline[j] = (T)std::max(rval, 0);
	  // G on rGr
	  g_evenline[j] = linebuf(i, j);
	  // B on rGr
	  bval =
	    - ((int)linebuf(i-2, j)>>3)
	    - ((int)linebuf(i-1, j-1)>>3) + ((int)linebuf(i-1, j)>>1) - ((int)linebuf(i-1, j+1)>>3)
	    + ((int)linebuf(i, j-2)>>4) + ((int)linebuf(i,j)>>1) + ((int)linebuf(i, j)>>3) + ((int)linebuf(i, j+2)>>4)
	    - ((int)linebuf(i+1, j-1)>>3) + ((int)linebuf(i+1, j)>>1) - ((int)linebuf(i+1, j+1)>>3)
	    - ((int)linebuf(i+2, j)>>3);
	  b_evenline[j] = (T)std::max(bval, 0);
      
	  j = x+1;
      
	  // R on gRg
	  r_evenline[j] = linebuf(i, j);
	  // G on gRg
	  gval = 
	    - ((int)linebuf(i-2, j)>>3)
	    + ((int)linebuf(i-1, j)>>2)
	    - ((int)linebuf(i, j-2)>>3) + ((int)linebuf(i, j-1)>>2) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j+1)>>2) - ((int)linebuf(i, j+2)>>3)
	    + ((int)linebuf(i+1, j)>>2)
	    - ((int)linebuf(i+2, j)>>3);
	  g_evenline[j] = (T)std::max(gval, 0);
	  // B on gRg
	  bval =
	    - ((int)linebuf(i-2, j)>>3) - ((int)linebuf(i-2, j)>>4)
	    + ((int)linebuf(i-1, j-1)>>2) + ((int)linebuf(i-1, j+1)>>2)
	    - ((int)linebuf(i, j-2)>>3) - ((int)linebuf(i, j-2)>>4) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>2) - ((int)linebuf(i, j+2)>>3) - ((int)linebuf(i, j+2)>>4)
	    + ((int)linebuf(i+1, j-1)>>2) + ((int)linebuf(i+1, j+1)>>2)
	    - ((int)linebuf(i+2, j)>>3) - ((int)linebuf(i+2, j)>>4);
	  b_evenline[j] = (T)std::max(bval, 0);
	}
      }
# pragma omp section
      {
	size_t i = y+1;
    
	// BGBG...
	T *r_oddline = r.getLine(i, 0);
	T *g_oddline = g.getLine(i, 0);
	T *b_oddline = b.getLine(i, 0);
    
	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  // R on gBg
	  int rval, gval, bval;
	  rval =
	    - ((int)linebuf(i-2, j)>>3) - ((int)linebuf(i-2, j)>>4)
	    + ((int)linebuf(i-1, j-1)>>2) + ((int)linebuf(i-1, j+1)>>2)
	    - ((int)linebuf(i, j-2)>>3) - ((int)linebuf(i, j-2)>>4) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>2) - ((int)linebuf(i, j+2)>>3) - ((int)linebuf(i, j+2)>>4)
	    + ((int)linebuf(i+1, j-1)>>2) + ((int)linebuf(i+1, j+1)>>2)
	    - ((int)linebuf(i+2, j)>>3) - ((int)linebuf(i+2, j)>>4);
	  r_oddline[j] = (T)std::max(rval, 0);
	  // G on gBg
	  gval =
	    - ((int)linebuf(i-2, j)>>3)
	    + ((int)linebuf(i-1, j)>>2)
	    - ((int)linebuf(i, j-2)>>3) + ((int)linebuf(i, j-1)>>2) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j+1)>>2) - ((int)linebuf(i, j+2)>>3)
	    + ((int)linebuf(i+1, j)>>2)
	    - ((int)linebuf(i+2, j)>>3);
	  g_oddline[j] = (T)std::max(gval, 0);
	  // B on gBg
	  b_oddline[j] = linebuf(i, j);
      
	  j = x+1;
      
	  // R on bGb
	  rval =
	    - ((int)linebuf(i-2, j)>>3)
	    - ((int)linebuf(i-1, j-1)>>3) + ((int)linebuf(i-1, j)>>1) - ((int)linebuf(i-1, j+1)>>3)
	    + ((int)linebuf(i, j-2)>>4) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>3) + ((int)linebuf(i, j+2)>>4)
	    - ((int)linebuf(i+1, j-1)>>3) + ((int)linebuf(i+1, j)>>1) - ((int)linebuf(i+1, j+1)>>3)
	    - ((int)linebuf(i+2, j)>>3);
	  r_oddline[j] = (T)std::max(rval, 0);
	  // G on bGb
	  g_oddline[j] = linebuf(i, j);
	  // B on bGb
	  bval =
	    + ((int)linebuf(i-2, j)>>4)
	    - ((int)linebuf(i-1, j-1)>>3) - ((int)linebuf(i-1, j+1)>>3)
	    - ((int)linebuf(i, j-2)>>3) + ((int)linebuf(i, j-1)>>1) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>3) + ((int)linebuf(i, j+1)>>1) - ((int)linebuf(i, j+2)>>3)
	    - ((int)linebuf(i+1, j-1)>>3) - ((int)linebuf(i+1, j+1)>>3)
	    + ((int)linebuf(i+2, j)>>4);
	  b_oddline[j] = (T)std::max(bval, 0);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
}

template <typename T>
void debayerEdgingGRBG2RGB(const cpixmap<T>& bayer, const int r_lower_limit, const int r_higher_limit, const int g_lower_limit, const int g_higher_limit, const int b_lower_limit, const int b_higher_limit, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits < std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));
  
  cslice<T> linebuf(bayer, 2/*lines*/, 2/*hpadding*/, 2/*vpadding*/);
  linebuf.draftSlice(bayer);

  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	/****** Even line ******/
	//** GRGR... (Even line)
	T *g_evenline = g.getLine(i, 0);
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  //// Even pixel
	  // G on rGr
	  g_evenline[j] = linebuf(i, j);
	  //// Odd pixel
	  size_t jj = j+1;
	  // G on gRg
	  int hdiff =
	    std::abs((int)(linebuf(i, jj-1)>>2) - (int)(linebuf(i, jj+1)>>2)) +
	    std::abs(-(int)(linebuf(i, jj-2)>>2) + (int)(linebuf(i, jj)>>1) - (int)(linebuf(i, jj+2)>>2));
	  int vdiff =
	    std::abs((int)(linebuf(i-1, jj)>>2) - (int)(linebuf(i+1, jj)>>2)) +
	    std::abs(-(int)(linebuf(i-2, jj)>>2) + (int)(linebuf(i, jj)>>1) - (int)(linebuf(i+2, jj)>>2));
	  int maxdiff = std::max(hdiff, vdiff);

	  int hweight =
	    -(int)(linebuf(i, jj-2)>>3)
	    +(int)(linebuf(i, jj-1)>>2)
	    +(int)(linebuf(i,   jj)>>2)
	    +(int)(linebuf(i, jj+1)>>2)
	    -(int)(linebuf(i, jj+2)>>3);
	  int vweight =
	    -(int)(linebuf(i-2, jj)>>3)
	    +(int)(linebuf(i-1, jj)>>2)
	    +(int)(linebuf(  i, jj)>>2)
	    +(int)(linebuf(i+1, jj)>>2)
	    -(int)(linebuf(i+2, jj)>>3);
	  int meanweight = (hweight>>1) + (vweight>>1);

	  int normalized;
	  if (maxdiff < g_lower_limit) normalized = 0;
	  else if (maxdiff > g_higher_limit) normalized = 256;
	  else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);
      
	  if (hdiff < vdiff) {
	    //g_evenline[jj] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	    int value = (((hweight-meanweight) * normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) g_evenline[jj] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) g_evenline[jj] = std::numeric_limits<T>::min();
	    else g_evenline[jj] = (T)value;
	  } else {
	    //g_evenline[jj] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	    int value = (((vweight-meanweight) * normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) g_evenline[jj] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) g_evenline[jj] = std::numeric_limits<T>::min();
	    else g_evenline[jj] = (T)value;
	  }
	}
      }
# pragma omp section
      {
	//// Odd line
	size_t ii = i+1;
	// BGBG...
	T *g_oddline = g.getLine(ii, 0);
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  //// Even pixel
	  // G on gBg
	  int hdiff =
	    std::abs((int)(linebuf(ii, j-1)>>2) - (int)(linebuf(ii, j+1)>>2)) +
	    std::abs(-(int)(linebuf(ii, j-2)>>2) + (int)(linebuf(ii, j)>>1) - (int)(linebuf(ii, j+2)>>2));
	  int vdiff =
	    std::abs((int)(linebuf(ii-1, j)>>2) - (int)(linebuf(ii+1, j)>>2)) +
	    std::abs(-(int)(linebuf(ii-2, j)>>2) + (int)(linebuf(ii, j)>>1) - (int)(linebuf(ii+2, j)>>2));
	  int maxdiff = std::max(hdiff, vdiff);

	  int hweight =
	    -(int)(linebuf(ii, j-2)>>3)
	    +(int)(linebuf(ii, j-1)>>2)
	    +(int)(linebuf(ii,   j)>>2)
	    +(int)(linebuf(ii, j+1)>>2)
	    -(int)(linebuf(ii, j+2)>>3);
	  int vweight =
	    -(int)(linebuf(ii-2, j)>>3)
	    +(int)(linebuf(ii-1, j)>>2)
	    +(int)(linebuf(  ii, j)>>2)
	    +(int)(linebuf(ii+1, j)>>2)
	    -(int)(linebuf(ii+2, j)>>3);
	  int meanweight = (hweight>>1) + (vweight>>1);

	  int normalized;
	  if (maxdiff < g_lower_limit) normalized = 0;
	  else if (maxdiff > g_higher_limit) normalized = 256;
	  else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);
      
	  if (hdiff < vdiff) {
	    //g_oddline[j] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	    int value = (((hweight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) g_oddline[j] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) g_oddline[j] = std::numeric_limits<T>::min();
	    else g_oddline[j] = (T)value;
	  } else {
	    //g_oddline[j] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	    int value = (((vweight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) g_oddline[j] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) g_oddline[j] = std::numeric_limits<T>::min();
	    else g_oddline[j] = (T)value;
	  }
	  //// Odd pixel
	  size_t jj = j+1;
	  // G on bGb
	  g_oddline[jj] = linebuf(ii, jj);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
  
  // processing for remained R, and B
  linebuf.setSlice(bayer, 2, 1, 1);
  linebuf.draftSlice(bayer);

  cslice<T> glinebuf(g, 2, 1, 1);
  glinebuf.draftSlice(g);
  
  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	/****** Even line ******/
	//** GRGR...
	T *r_evenline = r.getLine(i, 0);
	T *b_evenline = b.getLine(i, 0);
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  //// Even pixel
	  // R on rGr
	  int value =
	    +(int)(linebuf(i, j-1)>>1)
	    +(int)(linebuf(i, j+1)>>1)
	    -(int)(glinebuf(i, j-1)>>2)
	    +(int)(glinebuf(i, j)>>1)
	    -(int)(glinebuf(i, j+1)>>2);
	  if (value > std::numeric_limits<T>::max()) r_evenline[j] = std::numeric_limits<T>::max();
	  else if (value < std::numeric_limits<T>::min()) r_evenline[j] = std::numeric_limits<T>::min();
	  else r_evenline[j] = (T)value;
	  // B on rGr
	  value =
	    +(int)(linebuf(i-1, j)>>1)
	    +(int)(linebuf(i+1, j)>>1)
	    -(int)(glinebuf(i-1, j)>>2)
	    +(int)(glinebuf(i, j)>>1)
	    -(int)(glinebuf(i+1, j)>>2);
	  if (value > std::numeric_limits<T>::max()) b_evenline[j] = std::numeric_limits<T>::max();
	  else if (value < std::numeric_limits<T>::min()) b_evenline[j] = std::numeric_limits<T>::min();
	  else b_evenline[j] = (T)value;

	  //// Odd pixel
	  size_t jj = j+1;
	  // R on gRg
	  r_evenline[jj] = linebuf(i, jj);
	  // B on gRg
	  int d1diff =
	    std::abs((int)(linebuf(i-1, jj+1)>>2) - (int)(linebuf(i+1, jj-1)>>2)) +
	    std::abs(-(int)(glinebuf(i-1, jj+1)>>2) + (int)(glinebuf(i, jj)>>1) - (int)(glinebuf(i+1, jj-1)>>2));
	  int d2diff =
	    std::abs((int)(linebuf(i-1, jj-1)>>2) - (int)(linebuf(i+1, jj+1)>>2)) +
	    std::abs(-(int)(glinebuf(i-1, jj-1)>>2) + (int)(glinebuf(i, jj)>>1) - (int)(glinebuf(i+1, jj+1)>>2));
	  int maxdiff = std::max(d1diff, d2diff);

	  int d1weight =
	    +(int)(linebuf(i-1, jj+1)>>2)
	    +(int)(linebuf(i+1, jj-1)>>2)
	    -(int)(glinebuf(i-1, jj+1)>>3)
	    +(int)(glinebuf(i, jj)>>2)
	    -(int)(glinebuf(i+1, jj-1)>>3);
	  int d2weight =
	    +(int)(linebuf(i-1, jj-1)>>2)
	    +(int)(linebuf(i+1, jj+1)>>2)
	    -(int)(glinebuf(i-1, jj-1)>>3)
	    +(int)(glinebuf(i, jj)>>2)
	    -(int)(glinebuf(i+1, jj+1)>>3);
	  int meanweight = (d1weight>>1) + (d2weight>>1);

	  int normalized;
	  if (maxdiff < b_lower_limit) normalized = 0;
	  else if (maxdiff > b_higher_limit) normalized = 256;
	  else normalized = ((maxdiff - b_lower_limit)<<8) / (b_higher_limit - b_lower_limit);

	  if (d1diff < d2diff) {
	    //b_evenline[jj] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	    //b_evenline[jj] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	    value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) b_evenline[jj] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) b_evenline[jj] = std::numeric_limits<T>::min();
	    else b_evenline[jj] = (T)value;
	  } else {
	    //b_evenline[jj] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	    //b_evenline[jj] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	    value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) b_evenline[jj] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) b_evenline[jj] = std::numeric_limits<T>::min();
	    else b_evenline[jj] = (T)value;
	  }
	}
      }
# pragma omp section
      {
	/****** Odd line ******/
	size_t ii = i+1;
	//** BGBG...
	T *r_oddline = r.getLine(ii, 0);
	T *b_oddline = b.getLine(ii, 0);
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  //// Even pixel
	  // R on gBg
	  int d1diff =
	    std::abs((int)(linebuf(ii-1, j+1)>>2) - (int)(linebuf(ii+1, j-1)>>2)) +
	    std::abs(-(int)(glinebuf(ii-1, j+1)>>2) + (int)(glinebuf(ii, j)>>1) - (int)(glinebuf(ii+1, j-1)>>2));
	  int d2diff =
	    std::abs((int)(linebuf(ii-1, j-1)>>2) - (int)(linebuf(ii+1, j+1)>>2)) +
	    std::abs(-(int)(glinebuf(ii-1, j-1)>>2) + (int)(glinebuf(ii, j)>>1) - (int)(glinebuf(ii+1, j+1)>>2));
	  int maxdiff = std::max(d1diff, d2diff);

	  int d1weight =
	    +(int)(linebuf(ii-1, j+1)>>2)
	    +(int)(linebuf(ii+1, j-1)>>2)
	    -(int)(glinebuf(ii-1, j+1)>>3)
	    +(int)(glinebuf(ii, j)>>2)
	    -(int)(glinebuf(ii+1, j-1)>>3);
	  int d2weight =
	    +(int)(linebuf(ii-1, j-1)>>2)
	    +(int)(linebuf(ii+1, j+1)>>2)
	    -(int)(glinebuf(ii-1, j-1)>>3)
	    +(int)(glinebuf(ii, j)>>2)
	    -(int)(glinebuf(ii+1, j+1)>>3);
	  int meanweight = (d1weight>>1) + (d2weight>>1);

	  int normalized;
	  if (maxdiff < r_lower_limit) normalized = 0;
	  else if (maxdiff > r_higher_limit) normalized = 256;
	  else normalized = ((maxdiff - r_lower_limit)<<8) / (r_higher_limit - r_lower_limit);
      
	  if (d1diff < d2diff) {
	    //r_oddline[j] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	    //r_oddline[j] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	    int value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) r_oddline[j] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) r_oddline[j] = std::numeric_limits<T>::min();
	    else r_oddline[j] = (T)value;
	  } else {
	    //r_oddline[j] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	    //r_oddline[j] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	    int value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	    if (value > std::numeric_limits<T>::max()) r_oddline[j] = std::numeric_limits<T>::max();
	    else if (value < std::numeric_limits<T>::min()) r_oddline[j] = std::numeric_limits<T>::min();
	    else r_oddline[j] = (T)value;
	  }
	  // B on gBg
	  b_oddline[j] = linebuf(ii, j);
      
	  size_t jj = j+1;
	  //// Odd pixel
	  // R on bGb
	  int value =
	    +(int)(linebuf(ii-1, jj)>>1)
	    +(int)(linebuf(ii+1, jj)>>1)
	    -(int)(glinebuf(ii-1, jj)>>2)
	    +(int)(glinebuf(ii, jj)>>1)
	    -(int)(glinebuf(ii+1, jj)>>2);
	  if (value > std::numeric_limits<T>::max()) r_oddline[jj] = std::numeric_limits<T>::max();
	  else if (value < std::numeric_limits<T>::min()) r_oddline[jj] = std::numeric_limits<T>::min();
	  else r_oddline[jj] = (T)value;
	  // B on bGb 
	  value =
	    +(int)(linebuf(ii, jj-1)>>1)
	    +(int)(linebuf(ii, jj+1)>>1)
	    -(int)(glinebuf(ii, jj-1)>>2)
	    +(int)(glinebuf(ii, jj)>>1)
	    -(int)(glinebuf(ii, jj+1)>>2);
	  if (value > std::numeric_limits<T>::max()) b_oddline[jj] = std::numeric_limits<T>::max();
	  else if (value < std::numeric_limits<T>::min()) b_oddline[jj] = std::numeric_limits<T>::min();
	  else b_oddline[jj] = (T)value;
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
    glinebuf.shiftSlice(2, g);
  }
}

template <typename T>
void debayerBilinearRGGB2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));
  
  cslice<T> linebuf(bayer, 2/*lines*/, 1/*hpadding*/, 1/*vpadding*/);
  linebuf.draftSlice(bayer);

  for (size_t y = 0; y < bayer.getHeight(); y += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	//** RGRG...
	T *r_evenline = r.getLine(y, 0);
	T *g_evenline = g.getLine(y, 0);
	T *b_evenline = b.getLine(y, 0);

	size_t i = y;	
	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  // R on gRg
	  r_evenline[j] = linebuf(i, j);
	  // G on gRg
#if 1
	  T hdiff = std::abs((int)linebuf(i, j-1) - (int)linebuf(i, j+1));
	  T vdiff = std::abs((int)linebuf(i-1, j) - (int)linebuf(i+1, j));
	  if (hdiff < vdiff) g_evenline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  else g_evenline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
#else
	  g_evenline[j] = (linebuf(i-1, j)>>2) + (linebuf(i+1, j)>>2) + (linebuf(i, j-1)>>2) + (linebuf(i, j+1)>>2);
#endif
	  // B on gRg
	  b_evenline[j] = (linebuf(i-1, j-1)>>2) + (linebuf(i-1, j+1)>>2) + (linebuf(i+1, j-1)>>2) + (linebuf(i+1, j+1)>>2);

	  j = x+1;
	  // R on rGr
	  r_evenline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  // G on rGr
	  g_evenline[j] = linebuf(i, j);
	  // B on rGr
	  b_evenline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
	}
      }
# pragma omp section
      {
	//** GBGB...
	T *r_oddline = r.getLine(y+1, 0);
	T *g_oddline = g.getLine(y+1, 0);
	T *b_oddline = b.getLine(y+1, 0);

	size_t i = y+1;
	for (size_t x = 0; x < bayer.getWidth(); x += 2) {
	  size_t j = x;
	  // R on bGb
	  r_oddline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
	  // G on bGb
	  g_oddline[j] = linebuf(i, j);
	  // B on bGb
	  b_oddline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);

	  j = x+1;
	  // R on gBg
	  r_oddline[j] = (linebuf(i-1, j-1)>>2) + (linebuf(i-1, j+1)>>2) + (linebuf(i+1, j-1)>>2) + (linebuf(i+1, j+1)>>2);
	  // G on gBg
#if 1
	  T hdiff = std::abs((int)linebuf(i, j-1) - (int)linebuf(i, j+1));
	  T vdiff = std::abs((int)linebuf(i-1, j) - (int)linebuf(i+1, j));
	  if (hdiff < vdiff) g_oddline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  else g_oddline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
#else
	  g_oddline[j] = (linebuf(i, j-1)>>2) + (linebuf(i, j+1)>>2) + (linebuf(i-1, j)>>2) + (linebuf(i+1, j)>>2);
#endif
	  // B on gBg
	  b_oddline[j] = linebuf(i, j);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
}

template <typename T>
void debayerBilinearRGGB2RGB(cpixmap<T>& bayer, cpixmap<T>& rgb)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits);
  assert(rgb.getWidth() == bayer.getWidth());
  assert(rgb.getHeight() == bayer.getHeight());
  assert(rgb.getBands() >= 3);

  cslice<T> linebuf(bayer, 2/*lines*/, 1/*hpadding*/, 1/*vpadding*/);
  linebuf.draftSlice(bayer);
  
  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	//** RGRG...
	T *r_evenline = rgb.getLine(i, 2);
	T *g_evenline = rgb.getLine(i, 1);
	T *b_evenline = rgb.getLine(i, 0);
    
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  // R on gRg
	  r_evenline[j] = linebuf(i, j);
	  // G on gRg
#if 1
	  T hdiff = std::abs((int)linebuf(i, j-1) - (int)linebuf(i, j+1));
	  T vdiff = std::abs((int)linebuf(i-1, j) - (int)linebuf(i+1, j));
	  if (hdiff < vdiff) g_evenline[j] = (linebuf(i, j-1)>>1) + (linebuf(i, j+1)>>1);
	  else g_evenline[j] = (linebuf(i-1, j)>>1) + (linebuf(i+1, j)>>1);
#else
	  g_evenline[j] = (linebuf(i-1, j)>>2) + (linebuf(i+1, j)>>2) + (linebuf(i, j-1)>>2) + (linebuf(i, j+1)>>2);
#endif
	  // B on gRg
	  b_evenline[j] = (linebuf(i-1, j-1)>>2) + (linebuf(i-1, j+1)>>2) + (linebuf(i+1, j-1)>>2) + (linebuf(i+1, j+1)>>2);

	  size_t jj = j+1;
	  // R on rGr
	  r_evenline[jj] = (linebuf(i, jj-1)>>1) + (linebuf(i, jj+1)>>1);
	  // G on rGr
	  g_evenline[jj] = linebuf(i, jj);
	  // B on rGr
	  b_evenline[jj] = (linebuf(i-1, jj)>>1) + (linebuf(i+1, jj)>>1);
	}
      }
# pragma omp section
      {
	//** GBGB...
	T *r_oddline = rgb.getLine(i+1, 2);
	T *g_oddline = rgb.getLine(i+1, 1);
	T *b_oddline = rgb.getLine(i+1, 0);

	size_t ii = i+1;
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  // R on bGb
	  r_oddline[j] = (linebuf(ii-1, j) >> 1) + (linebuf(ii+1, j) >> 1);
	  // G on bGb
	  g_oddline[j] = linebuf(ii, j);
	  // B on bGb
	  b_oddline[j] = (linebuf(ii, j-1) >> 1) + (linebuf(ii, j+1) >> 1);

	  size_t jj = j+1;
	  // R on gBg
	  r_oddline[jj] = (linebuf(ii-1, jj-1) >> 2) + (linebuf(ii-1, jj+1) >> 2) + (linebuf(ii+1, jj-1) >> 2) + (linebuf(ii+1, jj+1) >> 2);
	  // G on gBg
#if 1
	  T hdiff = std::abs((int)linebuf(ii, jj-1) - (int)linebuf(ii, jj+1));
	  T vdiff = std::abs((int)linebuf(ii-1, jj) - (int)linebuf(ii+1, jj));
	  if (hdiff < vdiff) g_oddline[jj] = (linebuf(ii, jj-1)>>1) + (linebuf(ii, jj+1)>>1);
	  else g_oddline[jj] = (linebuf(ii-1, jj)>>1) + (linebuf(ii+1, jj)>>1);
#else
	  g_oddline[jj] = (linebuf(ii, jj-1)>>2) + (linebuf(ii, jj+1)>>2) + (linebuf(ii-1, jj)>>2) + (linebuf(ii+1, jj)>>2);
#endif
	  // B on gBg
	  b_oddline[jj] = linebuf(ii, jj);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
}


/* Referenced from
 * "High-quality linear interpolation for demosaicing of bayer-patterned color images"
 * - Henrique S. Malvar, Li-wei He, and Ross Cutler
 */
template <typename T>
void debayerHiQLinearRGGB2RGB(cpixmap<T>& bayer, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits <= std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));

  cslice<T> linebuf(bayer, 2/*lines*/, 2/*hpadding*/, 2/*vpadding*/);
  linebuf.draftSlice(bayer);

  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
    {
# pragma omp section
      {
	// RGRG...
	T *r_evenline = r.getLine(i, 0);
	T *g_evenline = g.getLine(i, 0);
	T *b_evenline = b.getLine(i, 0);

	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  int rval, gval, bval;
      
	  // R on gRg
	  r_evenline[j] = linebuf(i, j);
	  // G on gRg
	  gval = 
	    - ((int)linebuf(i-2, j)>>3)
	    + ((int)linebuf(i-1, j)>>2)
	    - ((int)linebuf(i, j-2)>>3) + ((int)linebuf(i, j-1)>>2) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j+1)>>2) - ((int)linebuf(i, j+2)>>3)
	    + ((int)linebuf(i+1, j)>>2)
	    - ((int)linebuf(i+2, j)>>3);
	  g_evenline[j] = (T)std::max(gval, 0);
	  // B on gRg
	  bval =
	    - ((int)linebuf(i-2, j)>>3) - ((int)linebuf(i-2, j)>>4)
	    + ((int)linebuf(i-1, j-1)>>2) + ((int)linebuf(i-1, j+1)>>2)
	    - ((int)linebuf(i, j-2)>>3) - ((int)linebuf(i, j-2)>>4) + ((int)linebuf(i, j)>>1) + ((int)linebuf(i, j)>>2) - ((int)linebuf(i, j+2)>>3) - ((int)linebuf(i, j+2)>>4)
	    + ((int)linebuf(i+1, j-1)>>2) + ((int)linebuf(i+1, j+1)>>2)
	    - ((int)linebuf(i+2, j)>>3) - ((int)linebuf(i+2, j)>>4);
	  b_evenline[j] = (T)std::max(bval, 0);

	  size_t jj = j+1;

	  // R on rGr
	  rval =
	    + ((int)linebuf(i-2, jj)>>4)
	    - ((int)linebuf(i-1, jj-1)>>3) - ((int)linebuf(i-1, jj+1)>>3)
	    - ((int)linebuf(i, jj-2)>>3) + ((int)linebuf(i, jj-1)>>1) + ((int)linebuf(i, jj)>>1) + ((int)linebuf(i, jj)>>3) + ((int)linebuf(i, jj+1)>>1) - ((int)linebuf(i, jj+2)>>3)
	    - ((int)linebuf(i+1, jj-1)>>3) - ((int)linebuf(i+1, jj+1)>>3)
	    + ((int)linebuf(i+2, jj)>>4);
	  r_evenline[jj] = (T)std::max(rval, 0);
	  // G on rGr
	  g_evenline[jj] = linebuf(i, jj);
	  // B on rGr
	  bval =
	    - ((int)linebuf(i-2, jj)>>3)
	    - ((int)linebuf(i-1, jj-1)>>3) + ((int)linebuf(i-1, jj)>>1) - ((int)linebuf(i-1, jj+1)>>3)
	    + ((int)linebuf(i, jj-2)>>4) + ((int)linebuf(i,jj)>>1) + ((int)linebuf(i, jj)>>3) + ((int)linebuf(i, jj+2)>>4)
	    - ((int)linebuf(i+1, jj-1)>>3) + ((int)linebuf(i+1, jj)>>1) - ((int)linebuf(i+1, jj+1)>>3)
	    - ((int)linebuf(i+2, jj)>>3);
	  b_evenline[jj] = (T)std::max(bval, 0);
	}
      }
# pragma omp section
      {
	size_t ii = i+1;
    
	// GBGB...
	T *r_oddline = r.getLine(ii, 0);
	T *g_oddline = g.getLine(ii, 0);
	T *b_oddline = b.getLine(ii, 0);
    
	for (size_t j = 0; j < bayer.getWidth(); j += 2) {
	  int rval, gval, bval;

	  // R on bGb
	  rval =
	    - ((int)linebuf(ii-2, j)>>3)
	    - ((int)linebuf(ii-1, j-1)>>3) + ((int)linebuf(ii-1, j)>>1) - ((int)linebuf(ii-1, j+1)>>3)
	    + ((int)linebuf(ii, j-2)>>4) + ((int)linebuf(ii, j)>>1) + ((int)linebuf(ii, j)>>3) + ((int)linebuf(ii, j+2)>>4)
	    - ((int)linebuf(ii+1, j-1)>>3) + ((int)linebuf(ii+1, j)>>1) - ((int)linebuf(ii+1, j+1)>>3)
	    - ((int)linebuf(ii+2, j)>>3);
	  r_oddline[j] = (T)std::max(rval, 0);
	  // G on bGb
	  g_oddline[j] = linebuf(ii, j);
	  // B on bGb
	  bval =
	    + ((int)linebuf(ii-2, j)>>4)
	    - ((int)linebuf(ii-1, j-1)>>3) - ((int)linebuf(ii-1, j+1)>>3)
	    - ((int)linebuf(ii, j-2)>>3) + ((int)linebuf(ii, j-1)>>1) + ((int)linebuf(ii, j)>>1) + ((int)linebuf(ii, j)>>3) + ((int)linebuf(ii, j+1)>>1) - ((int)linebuf(ii, j+2)>>3)
	    - ((int)linebuf(ii+1, j-1)>>3) - ((int)linebuf(ii+1, j+1)>>3)
	    + ((int)linebuf(ii+2, j)>>4);
	  b_oddline[j] = (T)std::max(bval, 0);
      
	  size_t jj = j+1;

	  // R on gBg
	  rval =
	    - ((int)linebuf(ii-2, jj)>>3) - ((int)linebuf(ii-2, jj)>>4)
	    + ((int)linebuf(ii-1, jj-1)>>2) + ((int)linebuf(ii-1, jj+1)>>2)
	    - ((int)linebuf(ii, jj-2)>>3) - ((int)linebuf(ii, jj-2)>>4) + ((int)linebuf(ii, jj)>>1) + ((int)linebuf(ii, jj)>>2) - ((int)linebuf(ii, jj+2)>>3) - ((int)linebuf(ii, jj+2)>>4)
	    + ((int)linebuf(ii+1, jj-1)>>2) + ((int)linebuf(ii+1, jj+1)>>2)
	    - ((int)linebuf(ii+2, jj)>>3) - ((int)linebuf(ii+2, jj)>>4);
	  r_oddline[jj] = (T)std::max(rval, 0);
	  // G on gBg
	  gval =
	    - ((int)linebuf(ii-2, jj)>>3)
	    + ((int)linebuf(ii-1, jj)>>2)
	    - ((int)linebuf(ii, jj-2)>>3) + ((int)linebuf(ii, jj-1)>>2) + ((int)linebuf(ii, jj)>>1) + ((int)linebuf(ii, jj+1)>>2) - ((int)linebuf(ii, jj+2)>>3)
	    + ((int)linebuf(ii+1, jj)>>2)
	    - ((int)linebuf(ii+2, jj)>>3);
	  g_oddline[jj] = (T)std::max(gval, 0);
	  // B on gBg
	  b_oddline[jj] = linebuf(ii, jj);
	}
      }
    }
    linebuf.shiftSlice(2, bayer);
  }
}

template <typename T>
void debayerEdgingRGGB2RGB(cpixmap<T>& bayer, int r_lower_limit, int r_higher_limit, int g_lower_limit, int g_higher_limit, int b_lower_limit, int b_higher_limit, cpixmap<T>& r, cpixmap<T>& g, cpixmap<T>& b)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits < std::numeric_limits<int>::digits);

  assert(bayer.isMatched(r));
  assert(bayer.isMatched(g));
  assert(bayer.isMatched(b));
  
  cslice<T> linebuf(bayer, 2, 2, 2);
  linebuf.draftSlice(bayer);

  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
  {
# pragma omp section
   {
    /****** Even line ******/
    //** RGRG... (Even line)
     T *g_evenline = g.getLine(i, 0);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // G on gRg
      int hdiff =
	std::abs((int)(linebuf(i, j-1)>>2) - (int)(linebuf(i, j+1)>>2)) +
	std::abs(-(int)(linebuf(i, j-2)>>2) + (int)(linebuf(i, j)>>1) - (int)(linebuf(i, j+2)>>2));
      int vdiff =
	std::abs((int)(linebuf(i-1, j)>>2) - (int)(linebuf(i+1, j)>>2)) +
	std::abs(-(int)(linebuf(i-2, j)>>2) + (int)(linebuf(i, j)>>1) - (int)(linebuf(i+2, j)>>2));
      int maxdiff = std::max(hdiff, vdiff);

      int hweight =
	-(int)(linebuf(i, j-2)>>3)
	+(int)(linebuf(i, j-1)>>2)
	+(int)(linebuf(i,   j)>>2)
	+(int)(linebuf(i, j+1)>>2)
	-(int)(linebuf(i, j+2)>>3);
      int vweight =
	-(int)(linebuf(i-2, j)>>3)
	+(int)(linebuf(i-1, j)>>2)
	+(int)(linebuf(  i, j)>>2)
	+(int)(linebuf(i+1, j)>>2)
	-(int)(linebuf(i+2, j)>>3);
      int meanweight = (hweight>>1) + (vweight>>1);

      int normalized;
      if (maxdiff < g_lower_limit) normalized = 0;
      else if (maxdiff > g_higher_limit) normalized = 256;
      else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);

      if (hdiff < vdiff) {
        //g_evenline[j] = (hweight * normalized + meanweight * (256 - normalized)) >> 8;
        //g_evenline[j] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((hweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_evenline[j] = std::numeric_limits<T>::min();
	else g_evenline[j] = (T)value;
      } else {
        //g_evenline[j] = (vweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_evenline[j] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((vweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_evenline[j] = std::numeric_limits<T>::min();
	else g_evenline[j] = (T)value;
      }
      //// Odd pixel
      size_t jj = j+1;
      // G on rGr
      g_evenline[jj] = linebuf(i, jj);
    }
   }
# pragma omp section
   {
    //// Odd line
    size_t ii = i+1;
    // GBGB...
    T *g_oddline = g.getLine(ii, 0);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // G on bGb
      g_oddline[j] = linebuf(ii, j);
      
      //// Odd pixel
      size_t jj = j+1;
      // G on gBg
      int hdiff =
	std::abs((int)(linebuf(ii, jj-1)>>2) - (int)(linebuf(ii, jj+1)>>2)) +
	std::abs(-(int)(linebuf(ii, jj-2)>>2) + (int)(linebuf(ii, jj)>>1) - (int)(linebuf(ii, jj+2)>>2));
      int vdiff =
	std::abs((int)(linebuf(ii-1, jj)>>2) - (int)(linebuf(ii+1, jj)>>2)) +
	std::abs(-(int)(linebuf(ii-2, jj)>>2) + (int)(linebuf(ii, jj)>>1) - (int)(linebuf(ii+2, jj)>>2));
      int maxdiff = std::max(hdiff, vdiff);

      int hweight =
	-(int)(linebuf(ii, jj-2)>>3)
	+(int)(linebuf(ii, jj-1)>>2)
	+(int)(linebuf(ii,   jj)>>2)
	+(int)(linebuf(ii, jj+1)>>2)
	-(int)(linebuf(ii, jj+2)>>3);
      int vweight =
	-(int)(linebuf(ii-2, jj)>>3)
	+(int)(linebuf(ii-1, jj)>>2)
	+(int)(linebuf(  ii, jj)>>2)
	+(int)(linebuf(ii+1, jj)>>2)
	-(int)(linebuf(ii+2, jj)>>3);
      int meanweight = (hweight>>1) + (vweight>>1);

      int normalized;
      if (maxdiff < g_lower_limit) normalized = 0;
      else if (maxdiff > g_higher_limit) normalized = 256;
      else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);
      
      if (hdiff < vdiff) {
	//g_oddline[jj] = (hweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_oddline[jj] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((hweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_oddline[jj] = std::numeric_limits<T>::min();
	else g_oddline[jj] = (T)value;
      } else {
	//g_oddline[jj] = (vweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_oddline[jj] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((vweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_oddline[jj] = std::numeric_limits<T>::min();
	else g_oddline[jj] = (T)value;
      }
    }
   }
  }
    linebuf.shiftSlice(2, bayer);
  }
  
  // processing for remained R, and B
  linebuf.setSlice(bayer, 2, 1, 1);
  linebuf.draftSlice(bayer);

  cslice<T> glinebuf(g, 2, 1, 1);
  glinebuf.draftSlice(g);
  
  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
  {
# pragma omp section
   {
    /****** Even line ******/
    //** RGRG...
    T *r_evenline = r.getLine(i, 0);
    T *b_evenline = b.getLine(i, 0);
    
    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // R on gRg
      r_evenline[j] = linebuf(i, j);
      // G on gRg
      // B on gRg
      int d1diff =
	std::abs((int)(linebuf(i-1, j+1)>>2) - (int)(linebuf(i+1, j-1)>>2)) +
	std::abs(-(int)(glinebuf(i-1, j+1)>>2) + (int)(glinebuf(i, j)>>1) - (int)(glinebuf(i+1, j-1)>>2));
      int d2diff =
	std::abs((int)(linebuf(i-1, j-1)>>2) - (int)(linebuf(i+1, j+1)>>2)) +
	std::abs(-(int)(glinebuf(i-1, j-1)>>2) + (int)(glinebuf(i, j)>>1) - (int)(glinebuf(i+1, j+1)>>2));
      int maxdiff = std::max(d1diff, d2diff);

      int d1weight =
	+(int)(linebuf(i-1, j+1)>>2)
	+(int)(linebuf(i+1, j-1)>>2)
	-(int)(glinebuf(i-1, j+1)>>3)
	+(int)(glinebuf(i, j)>>2)
	-(int)(glinebuf(i+1, j-1)>>3);
      int d2weight =
	+(int)(linebuf(i-1, j-1)>>2)
	+(int)(linebuf(i+1, j+1)>>2)
	-(int)(glinebuf(i-1, j-1)>>3)
	+(int)(glinebuf(i, j)>>2)
	-(int)(glinebuf(i+1, j+1)>>3);
      int meanweight = (d1weight>>1) + (d2weight>>1);

      int normalized;
      if (maxdiff < b_lower_limit) normalized = 0;
      else if (maxdiff > b_higher_limit) normalized = 256;
      else normalized = ((maxdiff - b_lower_limit)<<8) / (b_higher_limit - b_lower_limit);

      if (d1diff < d2diff) {
	//b_evenline[j] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	//b_evenline[j] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	int value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) b_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) b_evenline[j] = std::numeric_limits<T>::min();
	else b_evenline[j] = (T)value;
      } else {
	//b_evenline[j] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	//b_evenline[j] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	int value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) b_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) b_evenline[j] = std::numeric_limits<T>::min();
	else b_evenline[j] = (T)value;
      }

      //// Odd pixel
      size_t jj = j+1;
      // G on rGr
      // R on rGr
      int value =
	+(int)(linebuf(i, jj-1)>>1)
	+(int)(linebuf(i, jj+1)>>1)
	-(int)(glinebuf(i, jj-1)>>2)
	+(int)(glinebuf(i, jj)>>1)
	-(int)(glinebuf(i, jj+1)>>2);
      if (value > std::numeric_limits<T>::max()) r_evenline[jj] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) r_evenline[jj] = std::numeric_limits<T>::min();
      else r_evenline[jj] = (T)value;
      // B on rGr
      value =
	+(int)(linebuf(i-1, jj)>>1)
	+(int)(linebuf(i+1, jj)>>1)
	-(int)(glinebuf(i-1, jj)>>2)
	+(int)(glinebuf(i, jj)>>1)
	-(int)(glinebuf(i+1, jj)>>2);
      if (value > std::numeric_limits<T>::max()) b_evenline[jj] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) b_evenline[jj] = std::numeric_limits<T>::min();
      else b_evenline[jj] = (T)value;
    }
   }
# pragma omp section
   {
    /****** Odd line ******/
    size_t ii = i+1;
    //** GBGB...
    T *r_oddline = r.getLine(ii, 0);
    T *b_oddline = b.getLine(ii, 0);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // R on bGb
      int value =
	+(int)(linebuf(ii-1, j)>>1)
	+(int)(linebuf(ii+1, j)>>1)
	-(int)(glinebuf(ii-1, j)>>2)
	+(int)(glinebuf(ii, j)>>1)
	-(int)(glinebuf(ii+1, j)>>2);
      if (value > std::numeric_limits<T>::max()) r_oddline[j] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) r_oddline[j] = std::numeric_limits<T>::min();
      else r_oddline[j] = (T)value;
      // B on bGb 
      value =
	+(int)(linebuf(ii, j-1)>>1)
	+(int)(linebuf(ii, j+1)>>1)
	-(int)(glinebuf(ii, j-1)>>2)
	+(int)(glinebuf(ii, j)>>1)
	-(int)(glinebuf(ii, j+1)>>2);
      if (value > std::numeric_limits<T>::max()) b_oddline[j] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) b_oddline[j] = std::numeric_limits<T>::min();
      else b_oddline[j] = (T)value;
      
      size_t jj = j+1;
      //// Odd pixel
      // R on gBg
      int d1diff =
	std::abs((int)(linebuf(ii-1, jj+1)>>2) - (int)(linebuf(ii+1, jj-1)>>2)) +
	std::abs(-(int)(glinebuf(ii-1, jj+1)>>2) + (int)(glinebuf(ii, jj)>>1) - (int)(glinebuf(ii+1, jj-1)>>2));
      int d2diff =
	std::abs((int)(linebuf(ii-1, jj-1)>>2) - (int)(linebuf(ii+1, jj+1)>>2)) +
	std::abs(-(int)(glinebuf(ii-1, jj-1)>>2) + (int)(glinebuf(ii, jj)>>1) - (int)(glinebuf(ii+1, jj+1)>>2));
      int maxdiff = std::max(d1diff, d2diff);

      int d1weight =
	+(int)(linebuf(ii-1, jj+1)>>2)
	+(int)(linebuf(ii+1, jj-1)>>2)
	-(int)(glinebuf(ii-1, jj+1)>>3)
	+(int)(glinebuf(ii, jj)>>2)
	-(int)(glinebuf(ii+1, jj-1)>>3);
      int d2weight =
	+(int)(linebuf(ii-1, jj-1)>>2)
	+(int)(linebuf(ii+1, jj+1)>>2)
	-(int)(glinebuf(ii-1, jj-1)>>3)
	+(int)(glinebuf(ii, jj)>>2)
	-(int)(glinebuf(ii+1, jj+1)>>3);
      int meanweight = (d1weight>>1) + (d2weight>>1);

      int normalized;
      if (maxdiff < r_lower_limit) normalized = 0;
      else if (maxdiff > r_higher_limit) normalized = 256;
      else normalized = ((maxdiff - r_lower_limit)<<8) / (r_higher_limit - r_lower_limit);
      
      if (d1diff < d2diff) {
	//r_oddline[jj] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	//r_oddline[jj] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) r_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) r_oddline[jj] = std::numeric_limits<T>::min();
	else r_oddline[jj] = (T)value;
      } else {
	//r_oddline[jj] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	//r_oddline[jj] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) r_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) r_oddline[jj] = std::numeric_limits<T>::min();
	else r_oddline[jj] = (T)value;
      }
      // B on gBg
      b_oddline[jj] = linebuf(ii, jj);
    }
   }
  }
    linebuf.shiftSlice(2, bayer);
    glinebuf.shiftSlice(2, g);
  }
}

template <typename T>
void debayerEdgingRGGB2RGB(cpixmap<T>& bayer, int r_lower_limit, int r_higher_limit, int g_lower_limit, int g_higher_limit, int b_lower_limit, int b_higher_limit, cpixmap<T>& rgb)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits < std::numeric_limits<int>::digits);
  assert(rgb.getWidth() == bayer.getWidth());
  assert(rgb.getHeight() == bayer.getHeight());
  assert(rgb.getBands() >= 3);

  cslice<T> linebuf(bayer, 2, 2, 2);
  linebuf.draftSlice(bayer);

  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
  {
# pragma omp section
   {
    /****** Even line ******/
    //** RGRG... (Even line)
    T *g_evenline = rgb.getLine(i, 1);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // G on gRg
      int hdiff =
	std::abs((int)(linebuf(i, j-1)>>2) - (int)(linebuf(i, j+1)>>2)) +
	std::abs(-(int)(linebuf(i, j-2)>>2) + (int)(linebuf(i, j)>>1) - (int)(linebuf(i, j+2)>>2));
      int vdiff =
	std::abs((int)(linebuf(i-1, j)>>2) - (int)(linebuf(i+1, j)>>2)) +
	std::abs(-(int)(linebuf(i-2, j)>>2) + (int)(linebuf(i, j)>>1) - (int)(linebuf(i+2, j)>>2));
      int maxdiff = std::max(hdiff, vdiff);

      int hweight =
	-(int)(linebuf(i, j-2)>>3)
	+(int)(linebuf(i, j-1)>>2)
	+(int)(linebuf(i,   j)>>2)
	+(int)(linebuf(i, j+1)>>2)
	-(int)(linebuf(i, j+2)>>3);
      int vweight =
	-(int)(linebuf(i-2, j)>>3)
	+(int)(linebuf(i-1, j)>>2)
	+(int)(linebuf(  i, j)>>2)
	+(int)(linebuf(i+1, j)>>2)
	-(int)(linebuf(i+2, j)>>3);
      int meanweight = (hweight>>1) + (vweight>>1);

      int normalized;
      if (maxdiff < g_lower_limit) normalized = 0;
      else if (maxdiff > g_higher_limit) normalized = 256;
      else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);

      if (hdiff < vdiff) {
        //g_evenline[j] = (hweight * normalized + meanweight * (256 - normalized)) >> 8;
        //g_evenline[j] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = ((hweight-meanweight)>>7) * normalized + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_evenline[j] = std::numeric_limits<T>::min();
	else g_evenline[j] = (T)value;
      } else {
        //g_evenline[j] = (vweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_evenline[j] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = ((vweight-meanweight)>>7) * normalized + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_evenline[j] = std::numeric_limits<T>::min();
	else g_evenline[j] = (T)value;
      }
      //// Odd pixel
      size_t jj = j+1;
      // G on rGr
      g_evenline[jj] = linebuf(i, jj);
    }
   }
# pragma omp section
   {
    //// Odd line
    size_t ii = i+1;
    // GBGB...
    T *g_oddline = rgb.getLine(ii, 1);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // G on bGb
      g_oddline[j] = linebuf(ii, j);
      
      //// Odd pixel
      size_t jj = j+1;
      // G on gBg
      int hdiff =
	std::abs((int)(linebuf(ii, jj-1)>>2) - (int)(linebuf(ii, jj+1)>>2)) +
	std::abs(-(int)(linebuf(ii, jj-2)>>2) + (int)(linebuf(ii, jj)>>1) - (int)(linebuf(ii, jj+2)>>2));
      int vdiff =
	std::abs((int)(linebuf(ii-1, jj)>>2) - (int)(linebuf(ii+1, jj)>>2)) +
	std::abs(-(int)(linebuf(ii-2, jj)>>2) + (int)(linebuf(ii, jj)>>1) - (int)(linebuf(ii+2, jj)>>2));
      int maxdiff = std::max(hdiff, vdiff);

      int hweight =
	-(int)(linebuf(ii, jj-2)>>3)
	+(int)(linebuf(ii, jj-1)>>2)
	+(int)(linebuf(ii,   jj)>>2)
	+(int)(linebuf(ii, jj+1)>>2)
	-(int)(linebuf(ii, jj+2)>>3);
      int vweight =
	-(int)(linebuf(ii-2, jj)>>3)
	+(int)(linebuf(ii-1, jj)>>2)
	+(int)(linebuf(  ii, jj)>>2)
	+(int)(linebuf(ii+1, jj)>>2)
	-(int)(linebuf(ii+2, jj)>>3);
      int meanweight = (hweight>>1) + (vweight>>1);

      int normalized;
      if (maxdiff < g_lower_limit) normalized = 0;
      else if (maxdiff > g_higher_limit) normalized = 256;
      else normalized = ((maxdiff - g_lower_limit)<<8) / (g_higher_limit - g_lower_limit);
      
      if (hdiff < vdiff) {
	//g_oddline[jj] = (hweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_oddline[jj] = (hweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((hweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_oddline[jj] = std::numeric_limits<T>::min();
	else g_oddline[jj] = (T)value;
      } else {
	//g_oddline[jj] = (vweight * normalized + meanweight * (256 - normalized)) >> 8;
	//g_oddline[jj] = (vweight * normalized + meanweight * (256 - normalized)) >> 7;
	int value = (((vweight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) g_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) g_oddline[jj] = std::numeric_limits<T>::min();
	else g_oddline[jj] = (T)value;
      }
      // B on gBg
    }
   }
  }
    linebuf.shiftSlice(2, bayer);
  }
  
  // processing for remained R, and B
  linebuf.setSlice(bayer, 2, 1, 1);
  linebuf.draftSlice(bayer);

  cslice<T> glinebuf(rgb, 2, 1, 1);
  glinebuf.draftSlice(rgb, 1);
  
  for (size_t i = 0; i < bayer.getHeight(); i += 2) {
#pragma omp parallel sections
  {
# pragma omp section
   {
    /****** Even line ******/
    //** RGRG...
    T *r_evenline = rgb.getLine(i, 2);
    T *b_evenline = rgb.getLine(i, 0);
    
    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // R on gRg
      r_evenline[j] = linebuf(i, j);
      // B on gRg
      int d1diff =
	std::abs((int)(linebuf(i-1, j+1)>>2) - (int)(linebuf(i+1, j-1)>>2)) +
	std::abs(-(int)(glinebuf(i-1, j+1)>>2) + (int)(glinebuf(i, j)>>1) - (int)(glinebuf(i+1, j-1)>>2));
      int d2diff =
	std::abs((int)(linebuf(i-1, j-1)>>2) - (int)(linebuf(i+1, j+1)>>2)) +
	std::abs(-(int)(glinebuf(i-1, j-1)>>2) + (int)(glinebuf(i, j)>>1) - (int)(glinebuf(i+1, j+1)>>2));
      int maxdiff = std::max(d1diff, d2diff);

      int d1weight =
	+(int)(linebuf(i-1, j+1)>>2)
	+(int)(linebuf(i+1, j-1)>>2)
	-(int)(glinebuf(i-1, j+1)>>3)
	+(int)(glinebuf(i, j)>>2)
	-(int)(glinebuf(i+1, j-1)>>3);
      int d2weight =
	+(int)(linebuf(i-1, j-1)>>2)
	+(int)(linebuf(i+1, j+1)>>2)
	-(int)(glinebuf(i-1, j-1)>>3)
	+(int)(glinebuf(i, j)>>2)
	-(int)(glinebuf(i+1, j+1)>>3);
      int meanweight = (d1weight>>1) + (d2weight>>1);

      int normalized;
      if (maxdiff < b_lower_limit) normalized = 0;
      else if (maxdiff > b_higher_limit) normalized = 256;
      else normalized = ((maxdiff - b_lower_limit)<<8) / (b_higher_limit - b_lower_limit);

      if (d1diff < d2diff) {
	//b_evenline[j] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	//b_evenline[j] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	int value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) b_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) b_evenline[j] = std::numeric_limits<T>::min();
	else b_evenline[j] = (T)value;
      } else {
	//b_evenline[j] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	//b_evenline[j] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	int value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) b_evenline[j] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) b_evenline[j] = std::numeric_limits<T>::min();
	else b_evenline[j] = (T)value;
      }

      //// Odd pixel
      size_t jj = j+1;
      // R on rGr
      int value =
	+(int)(linebuf(i, jj-1)>>1)
	+(int)(linebuf(i, jj+1)>>1)
	-(int)(glinebuf(i, jj-1)>>2)
	+(int)(glinebuf(i, jj)>>1)
	-(int)(glinebuf(i, jj+1)>>2);
      if (value > std::numeric_limits<T>::max()) r_evenline[jj] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) r_evenline[jj] = std::numeric_limits<T>::min();
      else r_evenline[jj] = (T)value;
      // B on rGr
      value =
	+(int)(linebuf(i-1, jj)>>1)
	+(int)(linebuf(i+1, jj)>>1)
	-(int)(glinebuf(i-1, jj)>>2)
	+(int)(glinebuf(i, jj)>>1)
	-(int)(glinebuf(i+1, jj)>>2);
      if (value > std::numeric_limits<T>::max()) b_evenline[jj] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) b_evenline[jj] = std::numeric_limits<T>::min();
      else b_evenline[jj] = (T)value;
    }
   }
# pragma omp section
   {
    /****** Odd line ******/
    size_t ii = i+1;
    //** GBGB...
    T *r_oddline = rgb.getLine(ii, 2);
    T *b_oddline = rgb.getLine(ii, 0);

    for (size_t j = 0; j < bayer.getWidth(); j += 2) {
      //// Even pixel
      // R on bGb
      int value =
	+(int)(linebuf(ii-1, j)>>1)
	+(int)(linebuf(ii+1, j)>>1)
	-(int)(glinebuf(ii-1, j)>>2)
	+(int)(glinebuf(ii, j)>>1)
	-(int)(glinebuf(ii+1, j)>>2);
      if (value > std::numeric_limits<T>::max()) r_oddline[j] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) r_oddline[j] = std::numeric_limits<T>::min();
      else r_oddline[j] = (T)value;
      // G on bGb
      // B on bGb 
      value =
	+(int)(linebuf(ii, j-1)>>1)
	+(int)(linebuf(ii, j+1)>>1)
	-(int)(glinebuf(ii, j-1)>>2)
	+(int)(glinebuf(ii, j)>>1)
	-(int)(glinebuf(ii, j+1)>>2);
      if (value > std::numeric_limits<T>::max()) b_oddline[j] = std::numeric_limits<T>::max();
      else if (value < std::numeric_limits<T>::min()) b_oddline[j] = std::numeric_limits<T>::min();
      else b_oddline[j] = (T)value;
      
      size_t jj = j+1;
      //// Odd pixel
      // R on gBg
      int d1diff =
	std::abs((int)(linebuf(ii-1, jj+1)>>2) - (int)(linebuf(ii+1, jj-1)>>2)) +
	std::abs(-(int)(glinebuf(ii-1, jj+1)>>2) + (int)(glinebuf(ii, jj)>>1) - (int)(glinebuf(ii+1, jj-1)>>2));
      int d2diff =
	std::abs((int)(linebuf(ii-1, jj-1)>>2) - (int)(linebuf(ii+1, jj+1)>>2)) +
	std::abs(-(int)(glinebuf(ii-1, jj-1)>>2) + (int)(glinebuf(ii, jj)>>1) - (int)(glinebuf(ii+1, jj+1)>>2));
      int maxdiff = std::max(d1diff, d2diff);

      int d1weight =
	+(int)(linebuf(ii-1, jj+1)>>2)
	+(int)(linebuf(ii+1, jj-1)>>2)
	-(int)(glinebuf(ii-1, jj+1)>>3)
	+(int)(glinebuf(ii, jj)>>2)
	-(int)(glinebuf(ii+1, jj-1)>>3);
      int d2weight =
	+(int)(linebuf(ii-1, jj-1)>>2)
	+(int)(linebuf(ii+1, jj+1)>>2)
	-(int)(glinebuf(ii-1, jj-1)>>3)
	+(int)(glinebuf(ii, jj)>>2)
	-(int)(glinebuf(ii+1, jj+1)>>3);
      int meanweight = (d1weight>>1) + (d2weight>>1);

      int normalized;
      if (maxdiff < r_lower_limit) normalized = 0;
      else if (maxdiff > r_higher_limit) normalized = 256;
      else normalized = ((maxdiff - r_lower_limit)<<8) / (r_higher_limit - r_lower_limit);
      
      if (d1diff < d2diff) {
	//r_oddline[jj] = (d1weight * normalized + meanweight * (256 - normalized)) >> 8;
	//r_oddline[jj] = (d1weight * normalized - meanweight * (normalized - 256)) >> 7;
	value = (((d1weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) r_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) r_oddline[jj] = std::numeric_limits<T>::min();
	else r_oddline[jj] = (T)value;
      } else {
	//r_oddline[jj] = (d2weight * normalized + meanweight * (256 - normalized)) >> 8;
	//r_oddline[jj] = (d2weight * normalized - meanweight * (normalized - 256)) >> 7;
	value = (((d2weight-meanweight)*normalized)>>7) + (meanweight<<1);
	if (value > std::numeric_limits<T>::max()) r_oddline[jj] = std::numeric_limits<T>::max();
	else if (value < std::numeric_limits<T>::min()) r_oddline[jj] = std::numeric_limits<T>::min();
	else r_oddline[jj] = (T)value;
      }
      // B on gBg
      b_oddline[jj] = linebuf(ii, jj);
    }
   }
  }
    linebuf.shiftSlice(2, bayer);
    glinebuf.shiftSlice(2, rgb, 1);
  }
}

/* source(Gray): http://www.analog.com/media/en/technical-documentation/application-notes/EE358.pdf
   source(R): https://courses.cs.washington.edu/courses/cse467/08au/pdfs/lectures/09-Demosaicing.pdf
   constant Hue based interpolation: Hue does not change abruptly within a small neighborhood.
    1. interpolate green channel(gray in the case) first.
    2. interpolate hue (defined as either color differences or color ratios).
    3. Estimate the missing (red/blue) from the interpolated hue.
*/
template <typename T>
void debayerRCCC2GrayR(cpixmap<T>& rccc, cpixmap<T>& gray, cpixmap<T>& r)
{
  assert(std::numeric_limits<T>::is_integer);
  assert(!std::numeric_limits<T>::is_signed);
  assert(std::numeric_limits<T>::digits < std::numeric_limits<int>::digits);
  
  assert(rccc.getWidth() == gray.getWidth() && rccc.getHeight() == gray.getHeight());
  assert(rccc.getWidth() == r.getWidth() && rccc.getHeight() == r.getHeight());

  cslice<T> linebuf(rccc, 2/*lines*/, 2/*hpadding*/, 2/*vpadding*/);
  linebuf.draftSlice(rccc);

  for (size_t y = 0; y < rccc.getHeight(); y += 2) {
    T *grayline = gray.getLine(y);
    for (size_t x = 0; x < rccc.getWidth(); x += 2) {
      int value =
	-(linebuf(y-2, x)>>3) // -1/8
	+(linebuf(y-1, x)>>2) // +1/4
	-(linebuf(y, x-2)>>3) + (linebuf(y, x-1)>>2) + (linebuf(y, x)>>1) + (linebuf(y, x+1)>>2) - (linebuf(y, x+2)>>3)
	+(linebuf(y+1, x)>>2)
	-(linebuf(y+2, x)>>3);
      
      if (value < std::numeric_limits<T>::min()) grayline[x] = std::numeric_limits<T>::min();
      else if (value > std::numeric_limits<T>::max()) grayline[x] = std::numeric_limits<T>::max();
      else grayline[x] = (T)value;
    }
    linebuf.shiftSlice(2, rccc);
  }
  
}





