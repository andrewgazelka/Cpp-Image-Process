# CSCI 5607 Project 2

## Commands
Where `jf.jpg` is the jelly fish image

![Jellyfish](.github/JellyFish.jpg)


|Category|Command|Image|
|---|---|---|
|**Noise**|`image -input jf.jpg -noise 0.1 -output img.jpg` | ![Jellyfish noise](.github/noise0.1.jpg) |
||`image -input jf.jpg -noise 0.5 -output img.jpg` | ![Jellyfish noise](.github/noise0.5.jpg) |
|**Brightness**|`image -input jf.jpg -brightness 0.5 -output img.jpg` | ![Jellyfish dark](.github/bright0.5.jpg) |
||`image -input jf.jpg -brightness 2 -output img.jpg` | ![Jellyfish light](.github/bright2.jpg) |
|**Contrast**|`image -input jf.jpg -contrast 0.5 -output img.jpg` | ![Jellyfish light](.github/contrast0.5.jpg) |
||`image -input jf.jpg -contrast 1.5 -output img.jpg` | ![Jellyfish light](.github/contrast1.5.jpg) | 
|**Saturation**|`image -input jf.jpg -saturation 0.2 -output img.jpg` | ![Jellyfish light](.github/sat0.2.jpg) |
||`image -input jf.jpg -saturation 0.8 -output img.jpg` | ![Jellyfish light](.github/sat0.8.jpg) |
|**Crop**|`image -input jf.jpg -crop 600 300 100 200  -output img.jpg` | ![Jellyfish cropped](.github/crop.jpg) |
|**Extract Channel**|`image -input jf.jpg -extractChannel 0 -output img.jpg` | ![Jellyfish red](.github/ec0.jpg) |
|**Quantize**|`image -input jf.jpg -quantize 2 -output img.jpg` | ![Jellyfish red](.github/q2.jpg) |
|**Random Dither**|`image -input jf.jpg -randomDither 2 -output img.jpg` | ![Jellyfish red](.github/rd2.jpg) |
|**Blur**|`image -input jf.jpg -blur 2 -output img.jpg` | ![Jellyfish red](.github/bl2.jpg) |
|**Sharpen**|`image -input jf.jpg -sharpen 2 -output img.jpg` | ![Jellyfish red](.github/s2.jpg) |
|**Edge Detect**|`image -input jf.jpg -edgeDetect -output img.jpg` | ![Jellyfish red](.github/ed.jpg) |
|**Floyd-Steinberg Dither**|`image -input jf.jpg -FloydSteinbergDither 2 -output img.jpg` | ![Jellyfish red](.github/fs.jpg) |
|**Scale**|`image -input jf.jpg -scale 1.5 2.3 -output img.jpg` | ![Jellyfish red](.github/scale.jpg) |
|**Rotate**|`image -input jf.jpg -rotate 45 -output img.jpg` | ![Jellyfish red](.github/rot.jpg) |