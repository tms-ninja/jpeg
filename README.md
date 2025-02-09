# jpeg
A simple jpeg encoder for Python written in C++. Features include:
- encoding greyscale and colour images via the functions `encode_greyscale()` and `encode_colour()`
- Adjusting encoding quality by a quality factor using the [Independent JPEG Group](https://github.com/libjpeg-turbo/ijg) quality factor algorithm for generating quantization tables
- Ability to generate image specific Huffman tables
- Chromiance subsampling including 4:4:4, 4:2:2 and 4:2:0 modes

It can be built using `pip install .` and requires `numpy` as a build dependency.
