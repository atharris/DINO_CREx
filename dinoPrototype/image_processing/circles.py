import cv2
import image_processing_functions as ip


if __name__ == "__main__":
    img_folder = './bodies/'

    img = cv2.imread(img_folder + './Saturns-Crescent-Moons.jpg', cv2.IMREAD_GRAYSCALE)
    ip.hough_circles(img, center_dist=300, blur=3, show_img=True)

    img = cv2.imread(img_folder + './the-waning-crescent-moon-9.jpg', cv2.IMREAD_GRAYSCALE)
    ip.hough_circles(img, center_dist=300, blur=15, show_img=True)
