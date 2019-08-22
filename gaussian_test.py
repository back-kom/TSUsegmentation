import numpy as np
import cv2

img = cv2.imread('image/cat.jpg', cv2.IMREAD_GRAYSCALE)
img_out = img.copy()

height = img.shape[0]
width = img.shape[1]

gauss = (1.0 / 57) * np.array(
    [[0, 1, 2, 1, 1],
     [1, 3, 5, 3, 1],
     [2, 5, 9, 5, 2],
     [1, 3, 5, 3, 1],
     [0, 1, 2, 1, 0]])
sum(sum(gauss))
print(sum(sum(gauss)))
for i in np.arange(2, height - 2):
    for j in np.arange(2, width - 2):
        sum2 = 0
        for k in np.arange(-2, 3):
            for l in np.arange(-2, 3):
                a = img.item(i + k, j + l)
                p = gauss[2 + k, 2 + l]
                sum2 = sum2 + (p * a)
        b = sum2
        img_out.itemset((i, j), b)

cv2.imwrite('image/filter_gauss.jpg', img_out)

cv2.imshow('image', img_out)
cv2.waitKey(0)
