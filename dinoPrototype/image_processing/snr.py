import numpy as np
from scipy import misc
import math
import matplotlib.pyplot as plt
import image_processing_functions as imfunc
import search_location_functions as locfunc
import cv2
import time

cv2.setRNGSeed(int(time.time()))

extended = True
point = False

# function to get the power signal noise ratio of two images
# if the images are the same the result is undefined and 0 is returned
def getPSNR(img1, img2):
    i, j = img1.shape
    diff = np.empty(img1.shape)
    cv2.absdiff(img1, img2, diff)
    diff = np.square(diff)
    s = diff.sum()
    if s < 1e-10:
        return 0
    mse = s / (i*j)
    psnr = 10*np.log10(255**2 / mse)
    return psnr

# signal_threshold, noise_threshold, ROI_size (n x n pixel border), single side ROI_border_width
ROI_parameters = {}
ROI_parameters['signal_threshold'] = 1.5
ROI_parameters['noise_threshold'] = 1E-6
ROI_parameters['ROI_size'] = 100
ROI_parameters['ROI_border_width'] = 1
ROI_parameters['max_search_dist'] = 50

cam_res = (512, 512)
cam_pixel_size = (39E-6, 39E-6)  # horizontal, vertical [m]
cam_focal_length = .05  # [m]
cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0] / 2., cam_focal_length)),
           2 * math.degrees(math.atan2(cam_sensor_size[1] / 2., cam_focal_length)))

cameraParam = {}
cameraParam['resolution'] = cam_res
cameraParam['focal length'] = cam_focal_length
cameraParam['sensor size'] = cam_sensor_size
cameraParam['field of view'] = cam_fov
cameraParam['pixel size'] = cam_pixel_size


if extended:
    file_in = np.load('CDR_save_files/90_deg.npz')
    pos_beacon = np.vstack((file_in['earth_pos'], file_in['moon_pos']))
    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))
    ex_image = file_in['detector_array']
    ex_image = (ex_image / max(ex_image)) * 512
    ex_image = ex_image.reshape(512, 512)
    # plt.imshow(ex_image)
    # plt.show()

    # generate pixel line estimates for beacons in camera field of view
    pixel_line_beacon_i = locfunc.initial_beacons_estimate(
        pos_beacon, pos_sc, attde_sc, cameraParam)

    # pass in three initial estimates to test ROI generation logic
    pixel_truth = np.array([256, 394.99206801, 240])
    line_truth = np.array([256, 256, 245])
    truth_vals = np.stack((pixel_truth, line_truth))

    truth = [(256, 256), (256, 394.99206801), None]

    ROI_estimates = []
    for ind in range(len(pixel_truth)):
        ROI_estimates.append((pixel_truth[ind], line_truth[ind]))

    # crop original image to an ROI based on initial
    corner_ROI, image_ROI = imfunc.generate_point_source_ROI(ex_image, ROI_estimates, ROI_parameters)

    centers = np.empty([len(ROI_parameters)], dtype=tuple)
    for i in range(0, len(image_ROI)):

        # determine average value of region of interest border, subtract from rest of pixel map
        image_ROI[i] = imfunc.apply_ROI_border(image_ROI[i], ROI_parameters)

        # plt.savefig('saved_output/cropped_image_' + str(i) + '.png')
        # np.savez('saved_output/cropped_image_' + str(i)  + '.npz')

        # normalize ROI for center-finding function
        currentROI = np.empty(image_ROI[i].shape)
        cv2.normalize(image_ROI[i], currentROI, 255, 0, cv2.NORM_MINMAX)
        # print "current roi"
        # print currentROI
        # print np.max(currentROI)
        # currentROI = (image_ROI[i] / maxROIvalue) * 255
        # currentROI = np.uint8(currentROI)


        for run in range(0,6):
            psnr = []
            errors = []
            noiseScale = []
            noise_scale = 5
            noise_static = np.empty(currentROI.shape)
            cv2.randn(noise_static, 0, 0.2)
            while noise_scale <= 200:
                noise = noise_static * noise_scale
                # print "noise"
                # print noise

                image = np.empty(currentROI.shape)
                cv2.add(noise, currentROI, image)
                cv2.normalize(image, image, 255, 0, cv2.NORM_MINMAX)

                psnr.append(getPSNR(image, currentROI))

                # print "image"
                # print image
                # plt.imshow(currentROI)
                # plt.show()
                # if noise_scale == 200 or noise_scale == 0:
                    # plt.figure(20)
                    # plt.imshow(image)
                    # plt.show()

                center = imfunc.hough_circles(image, high_noise=True, center_dist=1, canny_thresh=175, blur=5, accum=7, show_img=False)
                if center is not None:
                    center = (center[0] + corner_ROI[i][1], center[1] + corner_ROI[i][0])
                    error = np.linalg.norm([(center[0] - pixel_truth[i]), (center[1] - line_truth[i])])
                else:
                    error = np.nan

                errors.append(error)

                noise_scale += 5

                if error is np.nan:
                    continue

                # if noise_scale == 200:
                #     crop_center = (int(round(center[0], 1)), \
                #                     int(round(center[1], 1)))
                #     crop_box = 30

                #     plt.figure(2)
                #     plt.imshow(ex_image[crop_center[1]-crop_box:crop_center[1]+crop_box,
                #                 crop_center[0]-crop_box:crop_center[0]+crop_box],
                #                 interpolation='none', cmap='viridis')

                #     # plot measured center
                #     plt.scatter(center[0]-crop_center[0]+crop_box,
                #                 center[1]-crop_center[1]+crop_box,
                #                 color='r', marker='x', s=75)

                #     # plot truth value
                #     plt.scatter(pixel_truth[i]-crop_center[0]+crop_box-.5,
                #                 line_truth[i]-crop_center[1]+crop_box-.5,
                #                 color='b', marker='o', s=75)
                #     plt.show()

            print "run: ", run
            print "psnr = ", psnr
            print "errors = ", errors
            not_found = 100 * np.count_nonzero(np.isnan(errors)) / float(len(errors))
            print "% of nans: ", not_found

            plt.figure(1)
            plt.scatter(psnr, errors, label='Centers not found: %.1f%%' % not_found)
            plt.xlabel("PSNR [dB]")
            plt.ylabel('Distance from true center [pixels]')

        plt.legend()
        plt.show()

if point:
    file_in = np.load('CDR_save_files/stars_only_cdr.npz')

    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))

    ex_image = file_in['detector_array']
    ex_image = (ex_image/np.amax(ex_image)) * 255
    ex_image = ex_image.reshape(512, 512)

    # generate pixel line estimates for stars in camera field of view
    # pixel_truth, line_star, star_catalog = locfunc.initial_stars_estimate(
    #     attde_sc, cameraParam, fname_catalog)

    # find centroid
    pixel_truth = np.array([443.71657484, 493.96318093, 95.3500384]) - 0.5
    line_truth = np.array([103.55964507, 105.8875924, 152.30559197]) - 0.5
    ROI_estimates = []
    for ind in range(len(pixel_truth)):
        ROI_estimates.append((pixel_truth[ind], line_truth[ind]))

# def find_centroid_point_source(pixel_map, pixel_line_beacon_i, ROI_parameters):
    loc_centroid = np.empty([len(ROI_estimates)], dtype=tuple)
    DN = np.empty([len(ROI_estimates)], dtype=int)

    # crop original image to an ROI based on initial
    corner_ROI, image_ROI = imfunc.generate_point_source_ROI(ex_image, ROI_estimates, ROI_parameters)

    for i in range(0, len(image_ROI)):

        # determine average value of region of interest border, subtract from rest of pixel map
        image_ROI[i] = imfunc.apply_ROI_border(image_ROI[i], ROI_parameters)

        currentROI = np.empty(image_ROI[i].shape)
        cv2.normalize(image_ROI[i], currentROI, 255, 0, cv2.NORM_MINMAX)
        for run in range(0,12):
            psnr = []
            errors = []
            noise_scale = 5
            noise_static = np.empty(currentROI.shape)
            cv2.randn(noise_static, 0, 0.1)

            while noise_scale <= 200:
                noise = noise_static * noise_scale

                image = np.empty(currentROI.shape)
                cv2.add(noise, currentROI, image)
                cv2.normalize(image, image, 255, 0, cv2.NORM_MINMAX)

                psnr.append(getPSNR(image, currentROI))

                # print "image"
                # print image
                # plt.imshow(currentROI)
                # plt.show()
                # if noise_scale == 200: # or noise_scale == 0:
                #     plt.figure(20)
                #     plt.imshow(image)
                #     plt.show()
                # calculate centroid position and ROI brightness value <-- centroid location is pixel number not index value
                # (starts with 1)
                center, _ = imfunc.find_centroid(image, corner_ROI[i], ROI_parameters)

                error = np.linalg.norm([(center[0] - pixel_truth[i]), (center[1] - line_truth[i])])
                errors.append(error)
                # print "error: ", error

                noise_scale += 5

                print(center)
                print(pixel_truth[i], line_truth[i])
                error = np.linalg.norm([(center[0] - pixel_truth[i]), (center[1] - line_truth[i])])
                print "error: ", error

                crop_center = (int(round(center[0], 1)),
                                int(round(center[1], 1)))
                print '\nSubplot Center Coordinate'
                print crop_center
                crop_box = 5
                plt.imshow(ex_image[crop_center[1]-crop_box:crop_center[1]+crop_box,
                        crop_center[0]-crop_box:crop_center[0]+crop_box],
                        interpolation='none', cmap='viridis')

                # plot measured centroid location
                plt.scatter(center[0]-crop_center[0]+crop_box,
                        center[1]-crop_center[1]+crop_box,
                                color='r', marker='x', s=75)

                # plot truth value
                plt.scatter(pixel_truth[i]-crop_center[0]+crop_box,
                                line_truth[i]-crop_center[1]+crop_box,
                                color='b', marker='+', s=75)
                # plt.scatter(0, 0, color='w', marker='+', s=75)
                plt.show()


            print "run: ", run
            print "psnr = ", psnr
            print "errors = ", errors

            plt.figure(1)
            plt.scatter(psnr, errors)
            plt.xlabel("PSNR [dB]")
            plt.ylabel('Distance from true center [pixels]')

        plt.figure(1)
        plt.legend()
        plt.show()
