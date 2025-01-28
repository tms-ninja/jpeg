import io
import pathlib
import unittest

import numpy as np
import PIL.Image

import jpeg

class Test_encode_greyscale(unittest.TestCase):
    """Tests encode_greyscale()"""

    def setUp(self):
        """Loads the pirate.png image"""

        test_dir = pathlib.Path(__file__).parent.parent

        pirate_path = test_dir / 'assets/pirate.png'

        with PIL.Image.open(pirate_path) as im:
            self.pirate_img = np.array(im)

    def test_mean_diff(self):
        """Tests the mean difference between the original and encoded images"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_greyscale(self.pirate_img)

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        # Now check the decoded image is similar to the original
        mean_diff = (decoded_image.astype(np.float64) - self.pirate_img.astype(np.float64)).mean()

        # Hard to say what a reasonable value for the expected mean difference is
        # Seems to be about 0.006 for the spec luminance table
        self.assertLess(mean_diff, 0.01)
    
    def test_reject_non_numpy_array(self):
        """Tests objects other than numpy arrays are rejected"""

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(1)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale([1, 2, 3])

    def test_reject_incorrect_element_type(self):
        """Tests numpy arrays with types other than np.uint8 are rejected"""

        test_array = np.array([[1, 2, 3, 4]])

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(test_array.astype(np.float64))

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(test_array.astype(np.int64))

    def test_reject_incorrect_shape(self):
        """Tests numpy arrays that are not 2d are rejected"""

        # 1d array
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(np.array(
                [1, 2, 3, 4],
                dtype=np.uint8
            ))

        # 3d array
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(np.array(
                [[[1, 2, 3, 4]]],
                dtype=np.uint8
            ))
