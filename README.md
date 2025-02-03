# jpeg
A simple jpeg encoder for Python written in C++. It supports encoding greyscale and colour images via the functions `encode_greyscale()` and `encode_colour()`. It further supports adjusting encoding quality by a quality factor using the [Independent JPEG Group](https://github.com/libjpeg-turbo/ijg) quality factor algorithm for generating quantization tables.

It can be built using `pip install .` and requires `numpy` as a build dependency.
