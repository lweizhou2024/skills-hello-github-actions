from PIL import Image
import pytesseract
import os
import cv2

# Path to the Tesseract executable
pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR\tesseract.exe'
os.environ['TESSDATA_PREFIX'] = r'C:\Program Files\Tesseract-OCR\tessdata'

# Load the image
image_path = 'test.png'
image_path = 'example.png'
image = cv2.imread(image_path)

x_start, y_start = 20, 41   # Top-left corner
x_end, y_end = 196, 100      # Bottom-right corner
x_start, y_start = 171, 97  # Top-left corner
x_end, y_end = 206, 142      # Bottom-right corner
x_start, y_start = 170, 248  # Top-left corner
x_end, y_end = 206, 275      # Bottom-right corner


# Crop the image
image = image[y_start:y_end, x_start:x_end]
#blurred = cv2.GaussianBlur(image, (5, 5), 0)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]

# Convert to PIL Image
image = Image.fromarray(thresh)

image.show()
# Use pytesseract to do OCR on the image
custom_config = r'--oem 3 --psm 6 outputbase digits'
text = pytesseract.image_to_string(image,  lang='lets', config='--psm 6 --oem 3 -c tessedit_char_whitelist=0123456789')
#text=pytesseract.image_to_string(image, config='--psm 6 --oem 3 -c tessedit_char_whitelist=0123456789')
#text=pytesseract.image_to_string(image, config='--psm 6 --oem 3')

# Print the extracted text
print("Extracted Numbers:", text)