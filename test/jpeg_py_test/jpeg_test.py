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
        mean_diff = abs((decoded_image.astype(np.float64) - self.pirate_img.astype(np.float64)).mean())

        # Hard to say what a reasonable value for the expected mean difference is
        # Seems to be about 0.006 for the spec luminance table
        self.assertLess(mean_diff, 0.01)

    def test_mean_diff_qf_90(self):
        """Tests the mean difference between the original and encoded images using a quality factor of 90"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_greyscale(self.pirate_img, qf=90)

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        # Now check the decoded image is similar to the original
        mean_diff = abs((decoded_image.astype(np.float64) - self.pirate_img.astype(np.float64)).mean())

        # Hard to say what a reasonable value for the expected mean difference is
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

    def test_reject_incorrect_quality_factor(self):
        """Tests quality factors less than 0 or greater than 100 are rejected"""

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(self.pirate_img, qf=-1)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_greyscale(self.pirate_img, qf=101)

        

class Test_encode_colour(unittest.TestCase):
    """Tests encode_colour()"""

    def setUp(self):
        """Loads the mandrill.png image"""

        test_dir = pathlib.Path(__file__).parent.parent

        pirate_path = test_dir / 'assets/mandrill.png'

        with PIL.Image.open(pirate_path) as im:
            img_data = np.array(im)

        self.mandrill_img = tuple(img_data[:, :, i] for i in range(3))

    def test_mean_diff(self):
        """Tests the mean difference between the original and encoded images"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_colour(*self.mandrill_img)

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        decoded_image = tuple(decoded_image[:, :, i] for i in range(3))

        # Now check the decoded image is similar to the original, component by component
        for c_orig, c_dec in zip(self.mandrill_img, decoded_image):
                
            mean_diff = abs((c_dec.astype(np.float64) - c_orig.astype(np.float64)).mean())

            # Hard to say what a reasonable value for the expected mean difference is
            # Max seems to be about 0.13 for the blue component using spec quantization tables
            self.assertLess(mean_diff, 1.0)

    def test_mean_diff_qf_90(self):
        """Tests the mean difference between the original and encoded images using a quality factor of 90"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_colour(*self.mandrill_img, qf=90)

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        decoded_image = tuple(decoded_image[:, :, i] for i in range(3))

        # Now check the decoded image is similar to the original, component by component
        for c_orig, c_dec in zip(self.mandrill_img, decoded_image):
                
            mean_diff = abs((c_dec.astype(np.float64) - c_orig.astype(np.float64)).mean())

            # Hard to say what a reasonable value for the expected mean difference is
            # Max seems to be about 0.13 for the blue component using spec quantization tables
            self.assertLess(mean_diff, 1.0)

    def test_mean_diff_ss_4_2_0(self):
        """Tests the mean difference between the original and encoded images using 4:2:0 subsampling"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_colour(*self.mandrill_img, ss="4:2:0")

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        decoded_image = tuple(decoded_image[:, :, i] for i in range(3))

        # Now check the decoded image is similar to the original, component by component
        for c_orig, c_dec in zip(self.mandrill_img, decoded_image):
                
            mean_diff = abs((c_dec.astype(np.float64) - c_orig.astype(np.float64)).mean())

            # Hard to say what a reasonable value for the expected mean difference is
            self.assertLess(mean_diff, 1.0)

    def test_mean_diff_ss_4_2_2(self):
        """Tests the mean difference between the original and encoded images using 4:2:2 subsampling"""

        # Encode and then decode the pirate image
        encoded_image = jpeg.encode_colour(*self.mandrill_img, ss="4:2:2")

        with PIL.Image.open(io.BytesIO(encoded_image)) as im:
            decoded_image = np.array(im)

        decoded_image = tuple(decoded_image[:, :, i] for i in range(3))

        # Now check the decoded image is similar to the original, component by component
        for c_orig, c_dec in zip(self.mandrill_img, decoded_image):
                
            mean_diff = abs((c_dec.astype(np.float64) - c_orig.astype(np.float64)).mean())

            # Hard to say what a reasonable value for the expected mean difference is
            self.assertLess(mean_diff, 1.0)
    
    def test_reject_non_numpy_array(self):
        """Tests objects other than numpy arrays are rejected"""

        valid_component = self.mandrill_img[0]

        # Test each component in turn
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(1, valid_component, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, 1, valid_component)
        
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, 1)

    def test_reject_incorrect_element_type(self):
        """Tests numpy arrays with types other than np.uint8 are rejected"""

        valid_component = self.mandrill_img[0]
        test_array = np.array([[1, 2, 3, 4]])

        # Test with floats
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(test_array.astype(np.float64), valid_component, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, test_array.astype(np.float64), valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, test_array.astype(np.float64))

        # Test with ints
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(test_array.astype(np.int64), valid_component, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, test_array.astype(np.int64), valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, test_array.astype(np.int64))

    def test_reject_incorrect_shape(self):
        """Tests numpy arrays that are not 2d are rejected"""

        valid_component = self.mandrill_img[0]

        # 1d array
        test_array_1d = np.array(
            [1, 2, 3, 4],
            dtype=np.uint8
        )

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(test_array_1d, valid_component, valid_component)
        
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, test_array_1d, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, test_array_1d)

        # 3d array
        test_array_3d = np.array(
            [[[1, 2, 3, 4]]],
            dtype=np.uint8
        )

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(test_array_3d, valid_component, valid_component)
        
        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, test_array_3d, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, test_array_3d)

    def test_reject_inconsistent_shape(self):
        """Tests components that have different shapes are rejected"""

        mandrill_shape = self.mandrill_img[0].shape

        array_bad_shape = self.mandrill_img[0][:mandrill_shape[0]//2, :mandrill_shape[1]//2]
        valid_component = self.mandrill_img[0]

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(array_bad_shape, valid_component, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, array_bad_shape, valid_component)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(valid_component, valid_component, array_bad_shape)

    def test_reject_incorrect_quality_factor(self):
        """Tests quality factors less than 0 or greater than 100 are rejected"""

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(*self.mandrill_img, qf=-1)

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(*self.mandrill_img, qf=101)

    def test_reject_incorrect_subsamplings(self):
        """Tests invalid subsamplings are rejected"""

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(*self.mandrill_img, ss="dsfadf")

        with self.assertRaises(ValueError) as context:
            _ = jpeg.encode_colour(*self.mandrill_img, ss="42:0")
