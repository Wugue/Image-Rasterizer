#include "texture.h"
#include "CGL/color.h"

namespace CGL {

Color Texture::sample(const SampleParams &sp) {
  // Part 5: Fill this in.
  if (sp.lsm == L_ZERO){
    if (sp.psm == P_NEAREST){
      return sample_nearest(sp.p_uv,0);
    } else {
      return sample_bilinear(sp.p_uv,0);
    }
  } else if (sp.lsm == L_NEAREST){
    if (sp.psm == P_NEAREST){
      return sample_nearest(sp.p_uv, get_level(sp));
    } else {
      return sample_bilinear(sp.p_uv, get_level(sp));
    }
  } else {
    return sample_trilinear(sp.p_uv, sp.p_dx_uv, sp.p_dy_uv);
  }
}

float Texture::get_level(const SampleParams &sp) {
  // Part 6: Fill this in.
  float l = fmax(sqrt(pow((sp.p_dx_uv[0]-sp.p_uv[0])*width,2) + pow((sp.p_dx_uv[1]-sp.p_uv[1])*height, 2)), sqrt(pow((sp.p_dy_uv[0]-sp.p_uv[0])*width, 2) + pow((sp.p_dy_uv[1]-sp.p_uv[1])*height, 2)));
  return log2(l);
}

Color Texture::sample_nearest(Vector2D uv, int level) {
  // Part 5: Fill this in.
  float width = mipmap[level].width;
  float height = mipmap[level].height;
  int x = floor(width * uv[0] + 0.5);
  int y = floor(height * uv[1] + 0.5);

  unsigned char killme [3] = {mipmap[level].texels[y*width*4 + x*4], mipmap[level].texels[y*width*4 + x*4 +1], mipmap[level].texels[y*width*4 + x*4 + 2]};
  return Color(killme);
}

Color Texture::sample_bilinear(Vector2D uv, int level) {
  // Part 5: Fill this in.
  // Used formula from https://helloacm.com/cc-function-to-compute-the-bilinear-interpolation/
  float width = mipmap[level].width;
  float height = mipmap[level].height;
  int x = floor(width * uv[0]);
  int y = floor(height * uv[1]);
  int x2 = floor(width * uv[0]) + 1;
  int y2 = floor(height * uv[1]) + 1;

  float x2x1 = x2 - x;
  float y2y1 = y2 - y;
  float x2x = x2 - uv[0]*width;
  float y2y = y2 - uv[1]*height;
  float yy1 = uv[1]*height - y;
  float xx1 = uv[0]*width - x;

  float p0 = x2x * y2y;
  float p1 = xx1 * y2y;
  float p2 = x2x * yy1;
  float p3 = xx1 * yy1;

  float r0 = mipmap[level].texels[y*width*4 + x*4];
  float r1 = mipmap[level].texels[y*width*4 + x2*4];
  float r2 = mipmap[level].texels[y2*width*4 + x*4];
  float r3 = mipmap[level].texels[y2*width*4 + x2*4];
  float g0 = mipmap[level].texels[y*width*4 + x*4+1];
  float g1 = mipmap[level].texels[y*width*4 + x2*4+1];
  float g2 = mipmap[level].texels[y2*width*4 + x*4+1];
  float g3 = mipmap[level].texels[y2*width*4 + x2*4+1];
  float b0 = mipmap[level].texels[y*width*4 + x*4+2];
  float b1 = mipmap[level].texels[y*width*4 + x2*4+2];
  float b2 = mipmap[level].texels[y2*width*4 + x*4+2];
  float b3 = mipmap[level].texels[y2*width*4 + x2*4+2];
  
  unsigned char killme [3] = {1.0*(p0*r0 + p1*r1 + p2*r2 + p3*r3)/(x2x1*y2y1), (1.0*p0*g0 + p1*g1 + p2*g2 + p3*g3)/(x2x1*y2y1), (1.0*p0*b0 + p1*b1 + p2*b2 + p3*b3)/(x2x1*y2y1)};
  return Color(killme);
}

Color Texture::sample_trilinear(Vector2D uv, Vector2D du, Vector2D dv) {
  // Part 6: Fill this in.
  float l = log2(fmax(sqrt(pow((du[0]-uv[0])*width,2) + pow((du[1]-uv[1])*height, 2)), sqrt(pow((dv[0]-uv[0])*width, 2) + pow((dv[1]-uv[1])*height, 2))));
  float high = ceil(l);
  float low = floor(l);
  float chigh = l - low;
  float clow = high - l;
  return chigh * sample_nearest(uv, high) + clow * sample_nearest(uv, low);
}




/****************************************************************************/



inline void uint8_to_float(float dst[4], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[4]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  dst_uint8[3] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[3])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 4; // 32 bit RGBA
    int currLevelPitch = currLevel.width * 4; // 32 bit RGBA

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[4];
    float input[4];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 4 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
          result[3] += wWeight[ii] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (4 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
          result[3] += hWeight[jj] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = result[3] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        4 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
              result[3] += weight * input[3];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 4 * i, result);
        }
      }
    }
  }
}

}
