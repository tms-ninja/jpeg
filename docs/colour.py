# Example of encoding a colour image
import pathlib

import numpy as np
import PIL.Image

import jpeg

docs_path = pathlib.Path(__file__).parent
asset_path = docs_path.parent / "test/assets"
input_path = asset_path / "mandrill.png"
output_path = docs_path / "mandrill.jpg"

with PIL.Image.open(input_path) as im:
    data = np.array(im)

# Unpack the image data
r, g, b = [data[:, :, i] for i in range(3)]

encoded_image = jpeg.encode_colour(r, g, b)

with open(output_path, 'wb') as f:
    f.write(encoded_image)

