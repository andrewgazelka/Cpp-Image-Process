//CSCI 5607 HW 2 - Image Conversion Instructor: S. J. Guy <sjguy@umn.edu>
//In this assignment you will load and convert between various image formats.
//Additionally, you will manipulate the stored image data by quantizing, cropping, and supressing channels

#include "image.h"
#include "Vec3D.h"
#include "FloatImg.h"
#include "Vec2D.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <iostream>

#include <fstream>

using namespace std;

/**
 * Image
 **/
Image::Image(int width_, int height_) {

    assert(width_ > 0);
    assert(height_ > 0);

    width = width_;
    height = height_;
    num_pixels = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels * 4];
    int b = 0; //which byte to write to
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
            data.raw[b++] = 0;
        }
    }

    assert(data.raw != NULL);
}

Image::Image(const Image &src) {
    width = src.width;
    height = src.height;
    num_pixels = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels * sizeof(Pixel)];

    memcpy(data.raw, src.data.raw, num_pixels * sizeof(Pixel));
}

Image::Image(char *fname) {

    int lastc = strlen(fname);
    int numComponents; //(e.g., Y, YA, RGB, or RGBA)
    data.raw = stbi_load(fname, &width, &height, &numComponents, 4);

    if (data.raw == NULL) {
        printf("Error loading image: %s", fname);
        exit(-1);
    }

    num_pixels = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

}

Image::~Image() {
    delete data.raw;
    data.raw = nullptr;
}

void Image::Write(char *fname) {

    int lastc = strlen(fname);

    switch (fname[lastc - 1]) {
        case 'g': //jpeg (or jpg) or png
            if (fname[lastc - 2] == 'p' || fname[lastc - 2] == 'e') //jpeg or jpg
                stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
            else //png
                stbi_write_png(fname, width, height, 4, data.raw, width * 4);
            break;
        case 'a': //tga (targa)
            stbi_write_tga(fname, width, height, 4, data.raw);
            break;
        case 'p': //bmp
        default:
            stbi_write_bmp(fname, width, height, 4, data.raw);
    }
}


void Image::Brighten(double factor) {
    int x, y;
    for (x = 0; x < Width(); x++) {
        for (y = 0; y < Height(); y++) {
            Pixel p = GetPixel(x, y);
            Pixel scaled_p = p * factor;
            GetPixel(x, y) = scaled_p;
        }
    }
}

void Image::ExtractChannel(int channel) {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto &pixel = GetPixel(x, y);
            pixel = pixel.GetChannel((char) channel);
        }
    }
}


void Image::Quantize(int nbits) {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto &pixel = GetPixel(x, y);
            pixel = PixelQuant(pixel, nbits);
        }
    }
}

Image *Image::Crop(int x, int y, int w, int h) const {

    int dx = w - x;
    int dy = h - y;

    auto newImage = new Image(dx, dy);

    for (int ix = 0; ix < dx; ++ix) {
        for (int iy = 0; iy < dy; ++iy) {
            auto pixel = GetPixel(ix + x, iy + y);
            newImage->SetPixel(ix, iy, pixel);
        }
    }
    return newImage;
}


void Image::AddNoise(double factor) {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto noisePixel = PixelNoise(factor);
            auto &originPixel = GetPixel(x, y);
            originPixel = originPixel + noisePixel;
        }
    }
}

void Image::ChangeContrast(double factor) {
    double meanLuminance = MeanLuminance();

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto &pixel = GetPixel(x, y);
            auto dLuminance = pixel.Luminance() - meanLuminance;
            auto scaleFactor = dLuminance * factor;
            pixel = pixel * scaleFactor;
        }
    }
}


void Image::ChangeSaturation(double factor) {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto &pixel = GetPixel(x, y);

            auto luminance = pixel.Luminance();

            auto dRed = (pixel.r - luminance) * factor;
            auto dGreen = (pixel.g - luminance) * factor;
            auto dBlue = (pixel.b - luminance) * factor;

            auto dPixel = Pixel(dRed, dGreen, dBlue);

            pixel = pixel + dPixel;
        }
    }
}


//For full credit, check that your dithers aren't making the pictures systematically brighter or darker
void Image::RandomDither(int nbits) {
    AddNoise(0.5); // TODO: ummm yea
    Quantize(nbits);
}

//This bayer method gives the quantization thresholds for an ordered dither.
//This is a 4x4 dither pattern, assumes the values are quantized to 16 levels.
//You can either expand this to a larger bayer pattern. Or (more likely), scale
//the threshold based on the target quantization levels.
static int Bayer4[4][4] = {
        {15, 7,  13, 5},
        {3,  11, 1,  9},
        {12, 4,  14, 6},
        {0,  8,  2,  10}
};


void Image::OrderedDither(int nbits) {
    // TODO
//    int n = 2;
//    for (int x = 0; x < width; ++x) {
//        for (int y = 0; y < height; ++y) {
//
//            auto &pixel = GetPixel(x, y);
//
//            int i = x % n;
//            int j = y % n;
//
//            int e = t
//
//        }
//    }
}

/* Error-diffusion parameters */
const double
        ALPHA = 7.0 / 16.0,
        BETA = 3.0 / 16.0,
        GAMMA = 5.0 / 16.0,
        DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits) {

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto oldPixel = GetPixel(x, y);

            auto newPixel = PixelQuant(oldPixel, nbits);
            GetPixel(x, y) = newPixel;

            auto errorVec = Vec3D::fromPixel(oldPixel) - Vec3D::fromPixel(newPixel);


            auto *alphaPixel = GetPixelOrNull(x + 1, y);
            if (alphaPixel != nullptr) {
                *alphaPixel = (Vec3D::fromPixel(*alphaPixel) + errorVec * ALPHA).toPixel();
            }

            auto *betaPixel = GetPixelOrNull(x - 1, y + 1);
            if (betaPixel != nullptr) {
                *betaPixel = (Vec3D::fromPixel(*betaPixel) + errorVec * BETA).toPixel();
            }


            auto *gammaPixel = GetPixelOrNull(x, y + 1);
            if (gammaPixel != nullptr) {
                *gammaPixel = (Vec3D::fromPixel(*gammaPixel) + errorVec * GAMMA).toPixel();
            }

            auto *deltaPixel = GetPixelOrNull(x + 1, y + 1);
            if (deltaPixel != nullptr) {
                *deltaPixel = (Vec3D::fromPixel(*deltaPixel) + errorVec * DELTA).toPixel();
            }
        }
    }

}

double gauss(double x, double y, double sigma) {
    double sigma2 = pow(sigma, 2);
    double r = -(x * x + y * y) / (2 * sigma2);
    return 1.0 / (2.0 * M_PI * sigma2) * exp(r);
}

double *gaussFilter(int radius, double sigma) {
    int width = radius * 2 + 1;
    int length = width * width;
    auto filter = new double[length];
    for (int x = -radius; x <= radius; ++x) {
        for (int y = -radius; y <= radius; ++y) {
            filter[(x + radius) + (y + radius) * width] = gauss(x, y, sigma);
        }
    }
    return filter;
}

void normalizeFilter(double *filter, int radius) {
    int width = radius * 2 + 1;
    int length = width * width;
    double sum = 0.0;
    for (int i = 0; i < length; ++i) {
        sum += filter[i];
    }
    double scale = 1.0 / sum; // how much to scale the filter so that no dark/brighten
    for (int i = 0; i < length; ++i) {
        filter[i] *= scale;
    }
}

void Image::Blur(int n) {
    int radius = 5;
    auto filter = gaussFilter(radius, n);
    normalizeFilter(filter, radius);
    ApplyFilter(filter, radius);
}

void Image::Sharpen(int n) {
    Image blurredCopy = *this;
    blurredCopy.Blur(n);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto &actualPixel = GetPixel(x, y);
            auto blurredPixel = blurredCopy.GetPixel(x, y);

            double dRed = actualPixel.r - blurredPixel.r;
            double dGreen = actualPixel.g - blurredPixel.g;
            double dBlue = actualPixel.b - blurredPixel.b;
            actualPixel.SetClamp(actualPixel.r + dRed, actualPixel.g + dGreen, actualPixel.b + dBlue);
        }
    }
}

void Image::EdgeDetect() {
    auto filter = new double[]{
            -1, -1, -1,
            -1, 8, -1,
            -1, -1, -1
    };
    ApplyFilter(filter, 1);
}

Image *Image::Scale(double sx, double sy) {
    int newWidth = (int) (sx * width);
    int newHeight = (int) (sy * height);

    auto *newImage = new Image(newWidth, newHeight);

    for (int x = 0; x < newWidth; ++x) {
        for (int y = 0; y < newHeight; ++y) {
            double prevX = x / sx;
            double prevY = y / sy;

            newImage->GetPixel(x, y) = Sample(prevX, prevY);
        }
    }
    return newImage;
}

Image *Image::Rotate(double degrees) {

    const auto angle = degrees * M_PI / 180.0;

    const auto center = Vec2D(width, height) / 2.0;

    // the four points of the frame centered around (0,0) and rotated and angle
    const auto p1 = (Vec2D(width, height) / 2.0).rotate(angle);
    const auto p2 = (Vec2D(-width, height) / 2.0).rotate(angle);
    const auto p3 = (Vec2D(-width, -height) / 2.0).rotate(angle);
    const auto p4 = (Vec2D(width, -height) / 2.0).rotate(angle);

    // after rotation the min/max X
    const double maxCenteredX = max({p1.x, p2.x, p3.x, p4.x});

    // after rotation the min/max Y
    const double maxCenteredY = max({p1.y, p2.y, p3.y, p4.y});

    // the new width, height of the frame
    const int newWidth = (int) (maxCenteredX * 2);
    const int newHeight = (int) (maxCenteredY * 2);

    const auto newImage = new Image(newWidth, newHeight);

    for (int finishX = 0; finishX < newWidth; ++finishX) {
        for (int finishY = 0; finishY < newHeight; ++finishY) {
            const Vec2D finish = Vec2D(finishX, finishY);
            const Vec2D relFinish = finish - Vec2D(maxCenteredX, maxCenteredY);
            const Vec2D relStart = relFinish.rotate(angle); // TODO: inverse rotate???
            const Vec2D start = relStart + center;

            newImage->GetPixel(finishX, finishY) = Sample(start.x, start.y);
        }
    }
    return newImage;
}

void Image::Fun() {
    /* WORK HERE */
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method) {
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}

double Image::MeanLuminance() {
    long summedLuminance = 0;
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            auto pixel = GetPixel(x, y);
            summedLuminance += pixel.Luminance();
        }
    }
    return (double) summedLuminance / NumPixels();
}

void Image::ApplyFilter(const double *filter, int radius) const {

    // because we need to make a copy to apply the filter
    auto newImage = Image(width, height);
    int filterWidth = 2 * radius + 1;

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {

            double rSum = 0;
            double gSum = 0;
            double bSum = 0;

            for (int dx = -radius; dx <= radius; ++dx) {
                for (int dy = -radius; dy <= radius; ++dy) {
                    double filterOn = filter[(dx + radius) + (dy + radius) * filterWidth];

                    int xOn = x + dx;
                    if (xOn < 0) xOn = 0;
                    if (xOn >= width) xOn = width - 1;

                    int yOn = y + dy;
                    if (yOn < 0) yOn = 0;
                    if (yOn >= height) yOn = height - 1;

                    auto pixelOn = GetPixel(xOn, yOn);

                    rSum += pixelOn.r * filterOn;
                    gSum += pixelOn.g * filterOn;
                    bSum += pixelOn.b * filterOn;
                }
            }

            Pixel newPixel;
            newPixel.SetClamp(rSum, gSum, bSum);

            newImage.GetPixel(x, y) = newPixel;
        }
    }

    // setting the old image to the same pixel values as the new one
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            GetPixel(x, y) = newImage.GetPixel(x, y);
        }
    }
}

Pixel Image::Sample(double x, double y) {

    double stdDev = 0.7;
    int r = 3;

    int iX = (int) x;
    int iY = (int) y;

    // allow for a little bit of leniency
    if (x < -1.1 || y <= -1.1 || x > width + 0.1 || y > height + 0.1) {
        return {0, 0, 0, 0}; // transparent pixel
    }

    double sumWeight = 0;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int dx = -r; dx <= r; ++dx) {
        for (int dy = -r; dy <= r; ++dy) {
            int onX = iX + dx;
            int onY = iY + dy;

            Pixel onPixel = GetPixelUnbounded(onX, onY);

            double trueDx = onX - x;
            double trueDy = onY - y;

            double weight;
            double dist2 = trueDy * trueDy + trueDx * trueDx;
            double manhattenDist = abs(trueDy) + abs(trueDx);
            switch (sampling_method) {
                case IMAGE_SAMPLING_POINT:
                    weight = manhattenDist <= 1.0 ? 1.0 : 0.0;
                    break;
                case IMAGE_SAMPLING_BILINEAR:
                    weight = dist2 <= 1.0 ? dist2 : 0.0;
                    break;
                case IMAGE_SAMPLING_GAUSSIAN:
                    weight = gauss(trueDx, trueDy, stdDev);
                    break;
            }
            sumWeight += weight;

            red += onPixel.r * weight;
            green += onPixel.g * weight;
            blue += onPixel.b * weight;

        }
    }
    Pixel pixel;
    pixel.SetClamp(red / sumWeight, green / sumWeight, blue / sumWeight);
    return pixel;
}


