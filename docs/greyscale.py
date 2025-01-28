import pathlib

import numpy as np
import PIL.Image

import jpeg

docs_path = pathlib.Path(__file__).parent
asset_path = docs_path.parent / "test/assets"
input_path = asset_path / "pirate.png"
output_path = docs_path / "pirate.jpg"

with PIL.Image.open(input_path) as im:
    data = np.array(im)

encoded_image = jpeg.encode_greyscale(data)

with open(output_path, 'wb') as f:
    f.write(encoded_image)

